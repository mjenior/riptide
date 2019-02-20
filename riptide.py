#!/usr/bin/python

# Dependencies
import copy
import time
import numpy
import cobra
import pandas
import bisect
import symengine
from cobra.util import solver
from optlang.symbolics import Zero
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis.sampling import ACHRSampler
from cobra.flux_analysis import flux_variability_analysis

# Read in transcriptomic read abundances, default is tsv with no header 
def read_transcription_file(read_abundances_file, header=False, replicates=False, sep='\t'):
    '''Generates dictionary of transcriptomic abundances from a file.
    
    Parameters
    ----------
    read_abundances_file : string
        User-provided file name which contains gene IDs and associated transcription values
    header : boolean
        Defines if read abundance file has a header that needs to be ignored
    replicates : boolean
        Defines if read abundances contains replicates and medians require calculation
    sep: string
        Defines what character separates entries on each line
    '''
    abund_dict = {}
    with open(read_abundances_file, 'r') as transcription:
        if header == True: header_line = transcription.readline()

        for line in transcription:
            line = line.split(sep)
            gene = str(line[0])
            
            if replicates == True:
                abundance = float(numpy.median([float(x) for x in line[1:]]))
            else:
                abundance = float(line[1])
            
            if gene in abund_dict.keys():
                abund_dict[gene] += abundance
            else:
                abund_dict[gene] = abundance

    return abund_dict


# Ensure that the user provided model and transcriptomic data are ready for RIPTiDe
def initialize_model(model):
    
    # Create a copy of the original model and set new id
    riptide_model = copy.deepcopy(model)
    riptide_model.id = str(riptide_model.id) + '_riptide'
    
    # Check that the model can grow
    solution = riptide_model.optimize()
    if solution.objective_value < 1e-6 or str(solution.objective_value) == 'nan':
        raise ValueError('ERROR: Provided model objective cannot carry flux! Please correct')
    
    # Calculate flux ranges and remove totally blocked reactions
    flux_span = flux_variability_analysis(riptide_model, fraction_of_optimum=0.1)
    flux_ranges = {}
    blocked_rxns = []
    for rxn_id, min_max in flux_span.iterrows():
        if max(abs(min_max)) < 1e-6:
            blocked_rxns.append(rxn_id)
        else:
            flux_ranges[rxn_id] = [min(min_max), max(min_max)]
    for rxn in blocked_rxns: 
        riptide_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
        
    # Calculate original solution space volume
    ellipsoid_vol = calculate_polytope_volume(flux_ranges)
    
    return riptide_model, ellipsoid_vol


# Converts a dictionary of transcript distribution percentiles
# Loosely based on:
# Schultz, A, & Qutub, AA (2016). Reconstruction of Tissue-Specific Metabolic Networks Using CORDA. 
#      PLoS Computational Biology. https://doi.org/10.1371/journal.pcbi.1004808
def assign_coefficients(raw_transcription_dict, model, percentiles, min_coefficients):
    
    # Screen transcriptomic abundances for genes that are included in model
    transcription_dict = {}
    for gene in model.genes:
        try:
            transcription_dict[gene.id] = raw_transcription_dict[gene.id]
        except KeyError:
            continue
    
    # Calculate transcript abundance cutoffs
    max_coefficients = min_coefficients[::-1]
    distribution = transcription_dict.values()
    abund_cutoffs = [numpy.percentile(distribution, x) for x in percentiles]
    
    # Screen transcript distribution by newly defined abundance intervals
    min_coefficient_dict = {}
    max_coefficient_dict = {}
    for gene in transcription_dict.iterkeys():
        transcription = transcription_dict[gene]
        if transcription in abund_cutoffs:
            index = abund_cutoffs.index(transcription)
            min_coefficient = min_coefficients[index]
            max_coefficient = max_coefficients[index]
        else:
            index = bisect.bisect_right(abund_cutoffs, transcription) - 1
            min_coefficient = min_coefficients[index]
            max_coefficient = max_coefficients[index]
            
        # Assign corresponding coefficients to reactions associated with each gene
        gene_rxn_dict = {}
        for rxn in list(model.genes.get_by_any(gene)[0].reactions):
            gene_rxn_dict[rxn.id] = gene
            
            if rxn.id in min_coefficient_dict.keys():
                min_coefficient_dict[rxn.id].append(min_coefficient)
                max_coefficient_dict[rxn.id].append(max_coefficient)
            else:
                min_coefficient_dict[rxn.id] = [min_coefficient]
                max_coefficient_dict[rxn.id] = [max_coefficient]
    
    # Assign final coefficients
    for rxn in model.reactions:
        try:
            # Take smallest value for reactions assigned multiple coefficients
            min_coefficient_dict[rxn.id] = min(min_coefficient_dict[rxn.id])
            max_coefficient_dict[rxn.id] = max(max_coefficient_dict[rxn.id])
        # Assign coefficients to reactions without genes
        except KeyError:
            min_coefficient_dict[rxn.id] = min_coefficients[2]
            max_coefficient_dict[rxn.id] = max_coefficients[2]
    
    return min_coefficient_dict, max_coefficient_dict, gene_rxn_dict


# Read in user defined reactions to keep or exclude
def incorporate_user_defined_reactions(rm_rxns, reaction_file):
    
    print('Integrating user definitions...')
    sep = ',' if '.csv' in str(reaction_file) else '\t'
    
    # Check if file actually exists    
    try:
        with open(reaction_file, 'r') as reactions:
            include_rxns = set(reactions.readline().split(sep))
            exclude_rxns = set(reactions.readline().split(sep))
    except FileNotFoundError:
        raise FileNotFoundError('ERROR: Defined reactions file not found! Please correct.')
        
    rm_rxns = rm_rxns.difference(include_rxns)
    rm_rxns |= exclude_rxns

    return rm_rxns


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
# Based on:
# Lewis NE, et al.(2010). Omic data from evolved E. coli are consistent with computed optimal growth from
#       genome-scale models. Molecular Systems Biology. 6, 390.
# Holzhütter, HG. (2004). The principle of flux minimization and its application to estimate 
#       stationary fluxes in metabolic networks. Eur. J. Biochem. 271; 2905–2922.
def constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth):
    
    with model as constrained_model:

        # Apply weigths to new expression
        pfba_expr = Zero
        for rxn in constrained_model.reactions:
            pfba_expr += coefficient_dict[rxn.id] * rxn.forward_variable
            pfba_expr += coefficient_dict[rxn.id] * rxn.reverse_variable
        
        # Calculate fractions of optimum used in each step
        fraction_1 = 1.0 - fraction
        half_fraction = fraction / 2.0
        fraction_2 = 1.0 - half_fraction
        fraction_3 = 1.0 + half_fraction
        
        # Calculate sum of fluxes constraint
        if sampling_depth == 'step_one':
            prev_obj_val = constrained_model.slim_optimize()
            # Set previous objective as a constraint, allow deviation
            prev_obj_constraint = constrained_model.problem.Constraint(constrained_model.objective.expression, lb=prev_obj_val*fraction_1, ub=prev_obj_val)
            constrained_model.add_cons_vars([prev_obj_constraint])
            constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='min', sloppy=True)
            constrained_model.solver.update()
            solution = constrained_model.optimize()
            
            # Determine reactions that do not carry any flux in the constrained model
            inactive_rxns = set([rxn.id for rxn in constrained_model.reactions if abs(solution.fluxes[rxn.id]) < 1e-6])
            return inactive_rxns
        
        else:
            # Explore solution space of constrained model with flux sampling, allow deviation
            constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='max', sloppy=True)
            solution = constrained_model.optimize()
            flux_sum_obj_val = solution.objective_value
            flux_sum_constraint = constrained_model.problem.Constraint(pfba_expr, lb=flux_sum_obj_val*fraction_2, ub=flux_sum_obj_val*fraction_3)
            constrained_model.add_cons_vars([flux_sum_constraint])
            constrained_model.solver.update()
            
            # Perform flux sampling (or FVA)
            flux_object = explore_flux_ranges(constrained_model, sampling_depth)
            return flux_object
    

# Prune model based on blocked reactions from minimization as well as user-defined reactions
def prune_model(new_model, rm_rxns, gene_rxn_dict, defined_rxns):
      
    # Integrate user definitions
    if defined_rxns != False: 
        rm_rxns = incorporate_user_defined_reactions(rm_rxns, defined_rxns)
        
    # Parse elements highlighted for pruning
    rm_genes = []
    nogene_rm_rxns = []
    for rxn in rm_rxns:
        try:
            rm_genes.append(gene_rxn_dict[rxn])
        except KeyError:
            nogene_rm_rxns.append(rxn)
            
    # Screen for duplicates
    rm_genes = list(set(rm_genes))
    nogene_rm_rxns = list(set(nogene_rm_rxns))
    
    # Prune all reactions associated with a deactivated gene
    remove_genes(new_model, rm_genes)
    # Prune inactive reactions without genes
    for rxn in nogene_rm_rxns:
        new_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
    
    # Prune possible residual orphans
    removed = 1
    while removed == 1:
        removed = 0
        for cpd in new_model.metabolites:
            if len(cpd.reactions) == 0:
                cpd.remove_from_model(); removed = 1
        for rxn in new_model.reactions:
            if len(rxn.metabolites) == 0: 
                rxn.remove_from_model(); removed = 1
    
    return new_model


# Analyze the possible ranges of flux in the constrained model
def explore_flux_ranges(model, samples):
    
    try:
        sampling_object = ACHRSampler(model)
        flux_object = sampling_object.sample(samples)        
        analysis = 'flux_sampling'
    except:
        # Handle errors for models that are now too small
        print('Constrained solution space too narrow for sampling, performing FVA instead')        
        flux_object = flux_variability_analysis(model, fraction_of_optimum=0.9)
        analysis = 'fva'
        
    return flux_object, analysis
    

# Constrain bounds for remaining reactions in model based on RIPTiDe results
def apply_bounds(constrained_model, flux_object):
    flux_ranges = {}
    
    # Handle FVA dataframe if necessary
    if len(flux_object.columns) == 2:
        for rxn in constrained_model.reactions:
            min_max = list(flux_object.loc[rxn.id])
            new_lb = min(min_max)
            new_ub = max(min_max)
            constrained_model.reactions.get_by_id(rxn.id).bounds = (new_lb, new_ub)
            flux_ranges[rxn.id] = [new_lb, new_ub]
    
    # Handle flux sampling results
    else:
        for rxn in constrained_model.reactions:
            distribution = list(flux_object[rxn.id])
            new_lb = min(distribution)
            new_ub = max(distribution)
            constrained_model.reactions.get_by_id(rxn.id).bounds = (new_lb, new_ub)
            flux_ranges[rxn.id] = [new_lb, new_ub]

    ellipsoid_vol = calculate_polytope_volume(flux_ranges)
        
    return constrained_model, ellipsoid_vol
    
    
# Calculate approximate volume of solution space, tretyed as ellipsoid
def calculate_polytope_volume(bounds):
    
    # Compile a list of radii from flux ranges
    radii = []
    for rxn in bounds.iterkeys():
        if bounds[rxn] == [0.0, 0.0]:
            continue
        else:
            diameter = abs(bounds[rxn][0]) + abs(bounds[rxn][1])
            radii.append(numpy.median(diameter) / 2.0)
    
    # Calculate volume 
    volume = (4.0/3.0) * numpy.pi * max(radii) * numpy.median(radii) * min(radii)
    volume = round(volume, 3)
    
    return volume


# Reports how long RIPTiDe took to run
def operation_report(start_time, model, riptide, old_vol, new_vol):
    
    # Pruning
    perc_removal = 100.0 - ((float(len(riptide.reactions)) / float(len(model.reactions))) * 100.0)
    perc_removal = round(perc_removal, 1)
    print('\nReactions pruned to ' + str(len(riptide.reactions)) + ' from ' + str(len(model.reactions)) + ' (' + str(perc_removal) + '% reduction)')
    perc_removal = 100.0 - ((float(len(riptide.metabolites)) / float(len(model.metabolites))) * 100.0)
    perc_removal = round(perc_removal, 1)
    print('Metabolites pruned to ' + str(len(riptide.metabolites)) + ' from ' + str(len(model.metabolites)) + ' (' + str(perc_removal) + '% reduction)')
    
    # Flux through objective
    new_ov = round(riptide.slim_optimize(), 3)
    old_ov = round(model.slim_optimize(), 3)
    per_shift = 100.0 - ((new_ov / old_ov) * 100.0)
    if per_shift == 0.0:
        print('\nNo change in flux through the objective')
    elif per_shift > 0.0:
        per_shift = round(abs(per_shift), 2)
        print('\nFlux through the objective REDUCED to ' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(per_shift) + '% shift)')
    elif per_shift < 0.0:
        per_shift = round(abs(per_shift), 2)
        print('\nFlux through the objective INCREASED to ' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(per_shift) + '% shift)')
    
    # Solution space volume
    if new_vol != 'none':
        vol_shift = 100.0 - ((new_vol / old_vol) * 100.0)
        if new_vol > 100000 or old_vol > 100000:
            pass
        elif vol_shift < 0.0:
            vol_shift = round(abs(vol_shift), 2)
            print('Solution space ellipsoid volume INCREASED to ~' + str(new_vol) + ' from ~' + str(old_vol) + ' (' + str(vol_shift) + '% shift)')
        elif vol_shift > 0.0:
            vol_shift = round(vol_shift, 2)
            print('Solution space ellipsoid volume DECREASED to ~' + str(new_vol) + ' from ~' + str(old_vol) + ' (' + str(vol_shift) + '% shift)')
        else:
            print('No change in Solution space volume')
    
    # Check that prune model can still achieve flux through the objective (just in case)
    if riptide.slim_optimize() < 1e-6 or str(riptide.slim_optimize()) == 'nan':
        print('\nWARNING: Contextualized model objective can no longer carry flux')
    
    # Run time
    duration = time.time() - start_time
    if duration < 60.0:
        duration = round(duration)
        print '\nRIPTiDe completed in ' + str(duration) + ' seconds'
    elif duration < 3600.0:
        duration = round((duration / 60.0), 1)
        print '\nRIPTiDe completed in ' + str(duration) + ' minutes'
    else:
        duration = round((duration / 3600.0), 1)
        print '\nRIPTiDe completed in ' + str(duration) + ' hours'


# Create context-specific model based on transcript distribution
def riptide(model, transcription, defined = False, sampling = 10000, percentiles = [50.0, 62.5, 75.0, 87.5], coefficients = [1.0, 0.5, 0.1, 0.01, 0.001], fraction = 0.8):
    '''Reaction Inclusion by Parsimony and Transcriptomic Distribution or RIPTiDe
    
    Creates a contextualized metabolic model based on parsimonious usage of reactions defined
    by their associated transcriptomic abundances. Returns a pruned, context-specific cobra.Model 
    and a pandas.DataFrame of associated flux sampling distributions

    Parameters
    ----------
    model : cobra.Model
        The model to be contextualized
    transcription : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
    defined : False or File
        Text file containing reactions IDs for forced inclusion listed on the first line and exclusion 
        listed on the second line (both .csv and .tsv formats supported)
    sampling : int or False
        Number of flux samples to collect, default is 10000, If False, sampling skipped
    percentiles : list of floats
        Percentile cutoffs of transcript abundance for linear coefficient assignments to associated reactions
        Defaults are [50.0, 62.5, 75.0, 87.5]
    coefficients : list of floats
        Linear coefficients to weight reactions based on distribution placement
        Defaults are [1.0, 0.5, 0.1, 0.01, 0.001]
    fraction : float
        Minimum percent of optimal objective value during FBA steps
        Default is 0.8
    '''
    start_time = time.time()
    
    # Correct some possible user error
    if sampling == False:
        pass
    elif sampling <= 0: 
        sampling = 10000
    else: 
        samples = int(sampling)
    if len(set(transcription.values())) == 1:
        raise ValueError('ERROR: All transcriptomic abundances are identical! Please correct')
    if len(coefficients) != len(percentiles) + 1:
        raise ValueError('ERROR: Invalid ratio of percentile cutoffs to linear coefficients! Please correct')
    fraction = float(fraction)
    if fraction <= 0.0:
        fraction = 0.8
    fraction = 1.0 - fraction
    
    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    print('Initializing model and parsing transcriptome...')
    riptide_model, orig_volume = initialize_model(model)
    min_coefficient_dict, max_coefficient_dict, gene_rxn_dict = assign_coefficients(transcription, riptide_model, percentiles, coefficients)
    
    # Prune now inactive network sections based on coefficients
    print('Pruning zero flux subnetworks...')
    rm_rxns = constrain_and_analyze_model(riptide_model, min_coefficient_dict, fraction, 'step_one')
    riptide_model = prune_model(riptide_model, rm_rxns, gene_rxn_dict, defined)
    
    # Find optimal solution space based on transcription and final constraints
    if sampling != False:
        print('Sampling context-specific solution space (longest step)...')
        flux_object, analysis_type = constrain_and_analyze_model(riptide_model, max_coefficient_dict, fraction, samples)
        
        # Constrain new model
        riptide_model, new_volume = apply_bounds(riptide_model, flux_object)
        operation_report(start_time, model, riptide_model, orig_volume, new_volume)
        return riptide_model, flux_object
    
    else:
        operation_report(start_time, model, riptide_model, orig_volume, 'none')
        return riptide_model
