#!/usr/bin/python

import sys
import copy
import time
import numpy
import cobra
import pandas
import bisect
import warnings
import symengine
from cobra.util import solver
from optlang.symbolics import Zero
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis import flux_variability_analysis, find_blocked_reactions


# Create a class to house riptide output data
class riptideClass:
    def __init__(self):
        self.model = 'NULL'
        self.transcriptome = 'NULL'
        self.coefficients = 'NULL'
        self.flux_samples = 'NULL'
        self.flux_variability = 'NULL'
        self.quantile_range = 'NULL'
        self.linear_coefficient_range = 'NULL'
        self.fraction_of_optimum = 'NULL'


# Create context-specific model based on transcript distribution
def contextualize(model, transcriptome, defined = False, samples = 500, 
    percentiles = [50.0, 62.5, 75.0, 87.5], coefficients = [1.0, 0.5, 0.1, 0.01, 0.001], 
    fraction = 0.75, conservative = False, objective = True, set_bounds = True):

    '''Reaction Inclusion by Parsimony and Transcriptomic Distribution or RIPTiDe
    
    Creates a contextualized metabolic model based on parsimonious usage of reactions defined
    by their associated transcriptomic abundances. Returns a pruned, context-specific cobra.Model 
    and a pandas.DataFrame of associated flux analysis

    Parameters
    ----------
    model : cobra.Model
        The model to be contextualized
        REQUIRED
    transcriptome : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
        REQUIRED
    defined : File
        Text file containing reactions IDs for forced inclusion listed on the first line and exclusion 
        listed on the second line (both .csv and .tsv formats supported)
    samples : int 
        Number of flux samples to collect, default is 500. If 0, sampling skipped FVA used instead
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
    objective : bool
        Sets previous objective function as a constraint with minimum flux equal to user input fraction
        Default is True
    set_bounds : bool
        Uses flax variability analysis results from constrained model to set new bounds for all equations
        Default is False
    '''

    start_time = time.time()
    riptide_object = riptideClass()
    
    # Correct some possible user error
    samples = int(samples)
    if samples <= 0: samples = 1
    if len(set(transcription.values())) == 1:
        raise ValueError('ERROR: All transcriptomic abundances are identical! Please correct')
    if len(coefficients) != len(percentiles) + 1:
        raise ValueError('ERROR: Invalid ratio of percentile cutoffs to linear coefficients! Please correct')
    fraction = float(fraction)
    if fraction <= 0.0:
        fraction = 0.8
    percentiles.sort() # sort ascending
    coefficients.sort(reverse=True) # sort descending
    solution = model.slim_optimize()
    if model.slim_optimize() < 1e-6 or str(model.slim_optimize()) == 'nan':
        raise ValueError('ERROR: Provided model objective cannot carry flux! Please correct')

    # Save parameters as part of the output object
    riptide_object.quantile_range = percentiles
    riptide_object.linear_coefficient_range = coefficients
    riptide_object.fraction_of_optimum = fraction
    riptide_object.transcriptome = transcriptome

    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    print('\nInitializing model and integrating transcriptomic data...')
    riptide_model = copy.deepcopy(model)
    riptide_model.id = str(riptide_model.id) + '_riptide'

    # Remove totally blocked reactions to speed up subsequent sections
    blocked_rxns = find_blocked_reactions(riptide_model)
    riptide_model = _prune_model(riptide_model, blocked_rxns, defined, conservative)
    coefficient_dict = _assign_coefficients(transcriptome, riptide_model, percentiles, coefficients)
    riptide_object.coefficients = coefficient_dict

    # Prune now inactive network sections based on coefficients
    print('Pruning zero flux subnetworks...')
    #iters = int(round(len(riptide_model.reactions) * 0.05)) + 1 # Adaptive to model size
    #rm_rxns = set([rxn.id for rxn in riptide_model.reactions])
    #for x in range(1, iters):
    #    curr_rxns = _constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, 0, objective)
    #    rm_rxns = rm_rxns.intersection(curr_rxns)
    rm_rxns = _constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, 0, objective)
    riptide_model = _prune_model(riptide_model, rm_rxns, defined, conservative)

    # Find optimal solution space based on transcription and final constraints
    print('Analyzing context-specific flux distributions...')
    flux_samples, fva_result = _constrain_and_analyze_model(riptide_model, coefficient_dict, fraction, samples, objective)
    riptide_object.flux_samples = flux_samples
    riptide_object.flux_variability = fva_result

    # Assign new reaction bounds
    if set_bounds == True:
        riptide_model = _set_new_bounds(riptide_model, fva_result, fraction)
    riptide_object.model = riptide_model

    # Analyze changes introduced by RIPTiDe and return results
    _operation_report(start_time, model, riptide_model)
    return riptide_object


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


# Converts a dictionary of transcript distribution percentiles
def _assign_coefficients(raw_transcription_dict, model, percentiles, min_coefficients):
    
    # Screen transcriptomic abundances for genes that are included in model
    transcription_dict = {}
    for gene in model.genes:
        try:
            transcription_dict[gene.id] = raw_transcription_dict[gene.id]
        except KeyError:
            continue
    
    # Calculate transcript abundance cutoffs
    distribution = list(transcription_dict.values())
    abund_cutoffs = [numpy.percentile(distribution, x) for x in percentiles]
    
    # Screen transcript distribution by newly defined abundance intervals
    coefficient_dict = {}
    for gene in transcription_dict.keys():
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
def _incorporate_user_defined_reactions(rm_rxns, reaction_file):
    
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
    rm_rxns = rm_rxns.union(exclude_rxns)

    return rm_rxns


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth, objective):
    
    with model as constrained_model:

        # Apply weigths to new expression
        pfba_expr = Zero
        if sampling_depth == 0:
            for rxn in constrained_model.reactions:
                pfba_expr += coefficient_dict[rxn.id] * rxn.forward_variable
                pfba_expr += coefficient_dict[rxn.id] * rxn.reverse_variable
        else:
            coeff_range = float(max(list(coefficient_dict.values()))) + float(min(list(coefficient_dict.values())))
            for rxn in constrained_model.reactions:
                max_coeff = coeff_range - float(coefficient_dict[rxn.id])
                pfba_expr += max_coeff * rxn.forward_variable
                pfba_expr += max_coeff * rxn.reverse_variable
        
        # Set previous objective as a constraint, allow deviation
        if objective == True:
            prev_obj_val = constrained_model.slim_optimize()
            prev_obj_constraint = constrained_model.problem.Constraint(constrained_model.objective.expression, lb=prev_obj_val*fraction, ub=prev_obj_val)
            constrained_model.add_cons_vars([prev_obj_constraint])
            constrained_model.solver.update()

        if sampling_depth == 0:
            # Determine reactions that do not carry any flux in highly constrained solution space
            constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='min', sloppy=True)
            constrained_model.solver.update()
            solution = constrained_model.optimize()
            inactive_rxns = set([rxn.id for rxn in constrained_model.reactions if abs(solution.fluxes[rxn.id]) <= 1e-6])
            
            return inactive_rxns
        
        else:
            # Explore solution space of constrained model with flux sampling, allow deviation
            constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='max', sloppy=True)
            constrained_model.solver.update()
            flux_sum_obj_val = constrained_model.slim_optimize()
            flux_sum_constraint = constrained_model.problem.Constraint(pfba_expr, lb=flux_sum_obj_val*fraction, ub=flux_sum_obj_val)
            constrained_model.add_cons_vars([flux_sum_constraint])
            constrained_model.solver.update()
            
            # Analyze flux ranges
            flux_samples = _gapsplit(constrained_model, n=sampling_depth)
            fva = flux_variability_analysis(constrained_model, fraction_of_optimum=fraction)

            return flux_samples, fva


# Prune model based on blocked reactions from minimization as well as user-defined reactions
def _prune_model(new_model, rm_rxns, defined_rxns, conserve):
      
    # Integrate user definitions
    if defined_rxns != False: 
        rm_rxns = _incorporate_user_defined_reactions(rm_rxns, defined_rxns)
        
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


# Use flux variability analysis on the constrained model to set new reaction bounds
def _set_new_bounds(model, fva, fraction):

    # Set new bounds for all reactions
    for rxn in model.reactions:
        fva_result = list(fva.loc[rxn.id])
        rxn.bounds = (min(fva_result), max(fva_result))

    return model


# Reports how long RIPTiDe took to run
def _operation_report(start_time, model, riptide):
    
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
    per_shift = 100.0 - ((float(new_ov) / float(old_ov)) * 100.0)
    if per_shift == 0.0:
        pass
    elif per_shift > 0.0:
        per_shift = round(abs(per_shift), 2)
        print('Flux through the objective DECREASED to ~' + str(new_ov) + ' from ~' + str(old_ov) + ' (' + str(per_shift) + '% change)')
    elif per_shift < 0.0:
        per_shift = round(abs(per_shift), 2)
        print('Flux through the objective INCREASED to ~' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(per_shift) + '% change)')
    
    # Check that prune model can still achieve flux through the objective (just in case)
    try:
        if riptide.slim_optimize() < 1e-6 or str(riptide.slim_optimize()) == 'nan':
            print('\nWARNING: Contextualized model objective can no longer carry flux')
    except:
        pass
    
    # Run time
    seconds = round(time.time() - start_time)
    if seconds < 60:
        print('\nRIPTiDe completed in ' + str(int(seconds)) + ' seconds\n')
    elif seconds < 3600:
        minutes, seconds = divmod(seconds, 60)
        print('\nRIPTiDe completed in ' + str(int(minutes)) + ' minutes and ' + str(int(seconds)) + ' seconds\n')
    else:
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        print('\nRIPTiDe completed in ' + str(int(hours)) + ' hours, ' + str(int(minutes)) + ' minutes, and ' + str(int(seconds)) + ' seconds\n')


#-----------------------------------------------------------------#

'''
gapsplit flux sampler
Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
bioRxiv 652917; doi: https://doi.org/10.1101/652917 
'''

def _gapsplit(
        model, 
        n=500, 
        max_tries=1000,
        primary_tol=0.001,
        secondary_frac=0.05,
        min_range=1e-5,
        enforce_range=True):
    """Randomly sample a COBRA model.

    Parameters
    ----------
    model: cobra.Model
        The model to sample. The model will not be modified during sampling.
    n: integer, default=500
        Number of samples to generate
    max_tries: integer, optional, default=1000
        Sampling attempts that return infeasible or unbounded solutions are
        discarded. Thus the total number of optimizations may exceed `n` for
        difficult models. `max_tries` limits the total number of attempts. If
        None (default), gapsplit will continue until `n` feasible samples are
        found.
    primary_tol: float, optional, default=0.001
        The primary target is split by setting the upper and lower bounds to
        the midway point of the max gap. The bounds are set to within +/-
        `primary_tol` times the width of the gap to avoid infeasible solutions
        due to numerical issues.
    secondary_frac: float, optional, default=0.05
        Fraction of model variables randomly chosen as secondary targets during
        each iteration. Default is 0.05 (5% of reactions). If 0, no secondary
        targeting is used; this may decrease coverage but improves runtime for
        numerically difficult models.
    min_range: float, optional, default=1e-5
        Variables are targeted only if their feasible range is larger than
        this value.
    enforce_range: boolean, optional, default=True
        If true (default), round solutions to fall within the feasible range.
        This prevents small deviations outside the feasible range from causing
        small decreases in coverage.

    Returns
    -------
    pandas.DataFrame
        A data frame with rows = samples and columns = reactions. This is the
        same format as the other cobrapy samplers.
    """

    warnings.filterwarnings('ignore') # Handle uninformative infeasible warning

    reactions = model.reactions
    fva = flux_variability_analysis(model, reactions, fraction_of_optimum=0.001)

    if secondary_frac >= 1.0:
        n_secondary = secondary_frac
    else:
        n_secondary = numpy.floor(secondary_frac * len(model.reactions)).astype(int)

    # only split reactions with feasible range >= min_range
    idxs = (fva.maximum - fva.minimum >= min_range).to_numpy().nonzero()[0]
    weights = (1.0 / (fva.maximum - fva.minimum) ** 2).to_numpy()

    samples = numpy.zeros((n, len(model.reactions)))
    k = 0
    infeasible_count = 0

    for try_ in range(max_tries):
        relative, target, width = _maxgap(samples[0:k,idxs], fva.iloc[idxs,:])
        
        primary_var = numpy.argmax(relative)
        primary_target = target[primary_var]
        primary_lb = primary_target - primary_tol*width[primary_var]
        primary_ub = primary_target + primary_tol*width[primary_var]

        secondary_vars = numpy.random.choice(len(idxs), n_secondary, replace=False)
        secondary_targets = target[secondary_vars]
        secondary_weights = weights[idxs[secondary_vars]]

        new_sample = _generate_sample(
            model, idxs[primary_var], primary_lb, primary_ub,
            idxs[secondary_vars], secondary_targets, secondary_weights)
        if new_sample is not None:
            if enforce_range:
                new_sample[new_sample > fva.maximum] = fva.maximum[new_sample > fva.maximum]
                new_sample[new_sample < fva.minimum] = fva.minimum[new_sample < fva.minimum]

            samples[k,:] = new_sample
            k += 1
        else:
            infeasible_count += 1

        if k >= n: break

    if k < n:
        # max_tries reached; return fewer than n samples
        samples = samples[:k,:]

    warnings.filterwarnings('default') # Return warnings to previous settings
    return pandas.DataFrame(data=samples,columns=fva.maximum.index)


def _generate_sample(
        model, primary_var, primary_lb, primary_ub,
        secondary_vars=None, secondary_targets=None, secondary_weights=None):
    """Formulate a [MI]QP to find a single solution."""
    with model:
        model.reactions[primary_var].lower_bound = primary_lb
        model.reactions[primary_var].upper_bound = primary_ub
        #model.objective = model.problem.Objective(0)
        model.objective = Zero
        solution = model.optimize()
        if solution.status != 'optimal':
            return None
        else:
            return solution.fluxes


def _maxgap(points, fva):
    # points has rows = samples, columns = variables

    # make a copy because we're going to sort the columns
    points = points.copy()
    points = numpy.vstack((fva.minimum, points, fva.maximum))
    points.sort(0)

    gaps = points[1:,:] - points[0:-1,:]
    width = gaps.max(0)
    loc = gaps.argmax(0)
    left = numpy.zeros(width.size)
    for i in range(width.size):
        left[i] = points[loc[i],i]

    relative = width / (points[-1,:] - points[0,:])
    target = left + width / 2.0

    return relative, target, width
