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


# Create class to house riptide 
class riptideClass:
    model = 'NULL'
    fluxes = 'NULL'
    flux_type = 'NULL'
    quantile_range = 'NULL'
    linear_coefficient_range = 'NULL'
    fraction_of_optimum = 'NULL'


# Read in transcriptomic read abundances, default is tsv with no header 
def read_transcription_file(read_abundances_file, header=False, replicates=False, sep='\t'):
    '''Generates dictionary of transcriptomic abundances from a file.
    
    Parameters
    ----------
    read_abundances_file : string
        User-provided file name which contains gene IDs and associated transcription values
    header : boolean
        Defines if read abundance file has a header that needs to be ignored
        default is no header
    replicates : boolean
        Defines if read abundances contains replicates and medians require calculation
        default is no replicates
    sep: string
        Defines what character separates entries on each line
        defaults to tab (.tsv)
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
    flux_span = flux_variability_analysis(riptide_model, fraction_of_optimum=0.01)
    flux_ranges = {}
    blocked_rxns = []
    for rxn_id, min_max in flux_span.iterrows():
        if max(abs(min_max)) < 1e-6:
            blocked_rxns.append(rxn_id)
        else:
            flux_ranges[rxn_id] = [min(min_max), max(min_max)]
    for rxn in blocked_rxns: 
        riptide_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)

    return riptide_model


# Converts a dictionary of transcript distribution percentiles
def assign_coefficients(raw_transcription_dict, model, percentiles, min_coefficients):
    
    # Screen transcriptomic abundances for genes that are included in model
    transcription_dict = {}
    for gene in model.genes:
        try:
            transcription_dict[gene.id] = raw_transcription_dict[gene.id]
        except KeyError:
            continue
    
    # Calculate transcript abundance cutoffs
    distribution = transcription_dict.values()
    abund_cutoffs = [numpy.percentile(distribution, x) for x in percentiles]
    
    # Screen transcript distribution by newly defined abundance intervals
    coefficient_dict = {}
    for gene in transcription_dict.iterkeys():
        transcription = transcription_dict[gene]
        if transcription in abund_cutoffs:
            index = abund_cutoffs.index(transcription)
            min_coefficient = min_coefficients[index]
        else:
            index = bisect.bisect_right(abund_cutoffs, transcription) - 1
            min_coefficient = min_coefficients[index]
                    
        # Assign corresponding coefficients to reactions associated with each gene
        for rxn in list(model.genes.get_by_any(gene)[0].reactions):            
            if rxn.id in coefficient_dict.keys():
                coefficient_dict[rxn.id].append(min_coefficient)
            else:
                coefficient_dict[rxn.id] = [min_coefficient]
    
    # Assign final coefficients
    nogene_coefficient = numpy.median(min_coefficients)
    for rxn in model.reactions:
        try:
            # Take smallest value for reactions assigned multiple coefficients
            coefficient_dict[rxn.id] = min(coefficient_dict[rxn.id])
        except KeyError:
            coefficient_dict[rxn.id] = nogene_coefficient
            continue
    
    return coefficient_dict


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
def constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth, bound_model):
    
    with model as constrained_model:

        # Apply weigths to new expression
        pfba_expr = Zero
        if sampling_depth == 'minimization':
            for rxn in constrained_model.reactions:
                pfba_expr += coefficient_dict[rxn.id] * rxn.forward_variable
                pfba_expr += coefficient_dict[rxn.id] * rxn.reverse_variable
        else:
            coeff_range = float(max(coefficient_dict.values())) + float(min(coefficient_dict.values()))
            for rxn in constrained_model.reactions:
                max_coeff = coeff_range - float(coefficient_dict[rxn.id])
                pfba_expr += max_coeff * rxn.forward_variable
                pfba_expr += max_coeff * rxn.reverse_variable
                
        # Calculate sum of fluxes constraint
        if sampling_depth == 'minimization':
            prev_obj_val = constrained_model.slim_optimize()
            # Set previous objective as a constraint, allow deviation
            prev_obj_constraint = constrained_model.problem.Constraint(constrained_model.objective.expression, lb=prev_obj_val*fraction, ub=prev_obj_val)
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
            flux_sum_constraint = constrained_model.problem.Constraint(pfba_expr, lb=flux_sum_obj_val*fraction, ub=flux_sum_obj_val)
            constrained_model.add_cons_vars([flux_sum_constraint])
            constrained_model.solver.update()
            
            # Perform flux sampling (or FVA)
            flux_object, flux_type = explore_flux_ranges(constrained_model, sampling_depth, fraction)

            if bound_model == False:
                return flux_object, flux_type
            else:
                constrained_model = apply_bounds(constrained_model, fraction)
                return flux_object, constrained_model, flux_type
    

# Prune model based on blocked reactions from minimization as well as user-defined reactions
def prune_model(new_model, rm_rxns, defined_rxns, conserve):
      
    # Integrate user definitions
    if defined_rxns != False: 
        rm_rxns = incorporate_user_defined_reactions(rm_rxns, defined_rxns)
        
    # Parse elements highlighted for pruning based on GPRs
    if conserve == True:
        final_rm_rxns = []
        for rxn in rm_rxns:
            test = 'pass'
            current_genes = list(new_model.reactions.get_by_id(rxn).genes)
            for gene in current_genes:
                for rxn_sub in gene.reactions:
                    if rxn_sub.id not in rm_rxns:
                        test = 'fail'
                    else:
                        pass
            
        	if test == 'pass': final_rm_rxns.append(rxn)
    else:
    	final_rm_rxns = rm_rxns
                        
    # Screen for duplicates
    final_rm_rxns = list(set(final_rm_rxns))
    
    # Prune inactive reactions
    for rxn in final_rm_rxns:
        new_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
    
    # Prune possible residual orphans, kind of sloppy but it's the only way 
    # I've found for it to actually thoroughly remove orphans
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
def explore_flux_ranges(model, samples, fraction):
    
    try:
        sampling_object = ACHRSampler(model)
        flux_object = sampling_object.sample(samples)        
        analysis = 'flux_sampling'
    except:
        # Handle errors for models that are now too small
        print('Constrained solution space too narrow for sampling, performing FVA instead')        
        flux_object = flux_variability_analysis(model, fraction_of_optimum=fraction)
        analysis = 'fva'
        
    return flux_object, analysis
  

# Apply new bounds to each reaction based on transcript-based constraints
def apply_bounds(constrained_model, fraction):

    fva_result = flux_variability_analysis(model, fraction_of_optimum=fraction)
    for rxn, min_max in fva_result.iterrows(): 
        constrained_model.reactions.get_by_id(rxn).bounds = (min(min_max), max(min_max))

    return constrained_model

  
# Calculate approximate volume of solution space, treated as ellipsoid
def calculate_polytope_volume(model, fraction):
    
    flux_span = flux_variability_analysis(model, fraction_of_optimum=fraction)
    bounds = {}
    for rxn_id, min_max in flux_span.iterrows(): bounds[rxn_id] = [min(min_max), max(min_max)]

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
    perc_removal = round(perc_removal, 2)
    print('\nReactions pruned to ' + str(len(riptide.reactions)) + ' from ' + str(len(model.reactions)) + ' (' + str(perc_removal) + '% change)')
    perc_removal = 100.0 - ((float(len(riptide.metabolites)) / float(len(model.metabolites))) * 100.0)
    perc_removal = round(perc_removal, 2)
    print('Metabolites pruned to ' + str(len(riptide.metabolites)) + ' from ' + str(len(model.metabolites)) + ' (' + str(perc_removal) + '% change)')
    
    # Flux through objective
    new_ov = round(riptide.slim_optimize(), 2)
    old_ov = round(model.slim_optimize(), 2)
    per_shift = 100.0 - ((new_ov / old_ov) * 100.0)
    if per_shift == 0.0:
        pass
    elif per_shift > 0.0:
        per_shift = round(abs(per_shift), 2)
        print('Flux through the objective DECREASED to ~' + str(new_ov) + ' from ~' + str(old_ov) + ' (' + str(per_shift) + '% change)')
    elif per_shift < 0.0:
        per_shift = round(abs(per_shift), 2)
        print('Flux through the objective INCREASED to ~' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(per_shift) + '% change)')
    
    # Solution space volume
    vol_shift = 100.0 - ((new_vol / old_vol) * 100.0)
    new_vol = round(new_vol, 2)
    old_vol = round(old_vol, 2)
    if new_vol > 100000 or old_vol > 100000:
        pass
    elif vol_shift < 0.0:
        vol_shift = round(abs(vol_shift), 2)
        print('Solution space volume INCREASED to ~' + str(new_vol) + ' from ~' + str(old_vol) + ' (' + str(vol_shift) + '% change)')
    elif vol_shift > 0.0:
        vol_shift = round(vol_shift, 2)
        print('Solution space volume DECREASED to ~' + str(new_vol) + ' from ~' + str(old_vol) + ' (' + str(vol_shift) + '% change)')
    
    # Check that prune model can still achieve flux through the objective (just in case)
    if riptide.slim_optimize() < 1e-6 or str(riptide.slim_optimize()) == 'nan':
        print('\nWARNING: Contextualized model objective can no longer carry flux')
    
    # Run time
    seconds = round(time.time() - start_time)
    if seconds < 60:
        print '\nRIPTiDe completed in ' + str(int(seconds)) + ' seconds\n'
    elif seconds < 3600:
        minutes, seconds = divmod(seconds, 60)
        print '\nRIPTiDe completed in ' + str(int(minutes)) + ' minutes and ' + str(int(seconds)) + ' seconds\n'
    else:
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        print '\nRIPTiDe completed in ' + str(int(hours)) + ' hours, ' + str(int(minutes)) + ' minutes, and ' + str(int(seconds)) + ' seconds\n'


# Create context-specific model based on transcript distribution
def riptide(model, transcription, defined = False, samples = 10000, percentiles = [50.0, 62.5, 75.0, 87.5], 
            coefficients = [1.0, 0.5, 0.1, 0.01, 0.001], fraction = 0.8, conservative = False, bound = False):
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
    defined : File
        Text file containing reactions IDs for forced inclusion listed on the first line and exclusion 
        listed on the second line (both .csv and .tsv formats supported)
    samples : int 
        Number of flux samples to collect, default is 10000, If 0, sampling skipped
    percentiles : list
        Percentile cutoffs of transcript abundance for linear coefficient assignments to associated reactions
        Default is [50.0, 62.5, 75.0, 87.5]
    coefficients : list
        Linear coefficients to weight reactions based on distribution placement
        Default is [1.0, 0.5, 0.1, 0.01, 0.001]
    fraction : float
        Minimum percent of optimal objective value during FBA steps
        Default is 0.8
    conservative : bool
        Conservatively remove inactive reactions based on genes
        Default is False
    bound : bool
        Bounds each reaction based on transcriptomic constraints
        Default is False
    '''

    start_time = time.time()
    riptide_object = riptideClass()
    
    # Correct some possible user error
    samples = int(samples)
    if samples < 0: samples = 0
    if len(set(transcription.values())) == 1:
        raise ValueError('ERROR: All transcriptomic abundances are identical! Please correct')
    if len(coefficients) != len(percentiles) + 1:
        raise ValueError('ERROR: Invalid ratio of percentile cutoffs to linear coefficients! Please correct')
    fraction = float(fraction)
    if fraction <= 0.0:
        fraction = 0.8
    percentiles.sort() # sort ascending
    coefficients.sort(reverse=True) # sort descending

    riptide_object.quantile_range = percentiles
    riptide_object.linear_coefficient_range = coefficients
    riptide_object.fraction_of_optimum = fraction

    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    print('\nInitializing model and parsing transcriptome...')
    riptide_model = initialize_model(model)
    orig_volume = calculate_polytope_volume(riptide_model, fraction)
    coefficient_dict = assign_coefficients(transcription, riptide_model, percentiles, coefficients)
    
    # Prune now inactive network sections based on coefficients
    print('Pruning zero flux subnetworks...')
    rm_rxns = constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, 'minimization', bound)
    riptide_model = prune_model(riptide_model, rm_rxns, defined, conservative)
    new_volume = calculate_polytope_volume(riptide_model, fraction)

    # Find optimal solution space based on transcription and final constraints
    if samples != 0:
        print('Sampling context-specific flux distributions (longest step)...')
        if bound == False:
            flux_object, flux_type = constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, samples, bound)
        else:
            flux_object, riptide_model, flux_type = constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, samples, bound)

        riptide_object.model = riptide_model
        riptide_object.fluxes = flux_object
        riptide_object.flux_type = flux_type
    else:
        riptide_object.model = riptide_model
        riptide_object.flux_type = flux_type

    # Analyze changes introduced by RIPTiDe and return results
    operation_report(start_time, model, riptide_model, orig_volume, new_volume)    
    return riptide_object

