#!/usr/bin/python

import copy
import time
import numpy
import cobra
import pandas
import warnings
import symengine
import itertools
from cobra.util import solver
from scipy.stats import spearmanr
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis import flux_variability_analysis, find_blocked_reactions


# Create a class to house riptide output data
class riptideClass:
    def __init__(self):
        self.model = 'NULL'
        self.transcriptome = 'NULL'
        self.minimization_coefficients = 'NULL'
        self.maximization_coefficients = 'NULL'
        self.flux_samples = 'NULL'
        self.flux_variability = 'NULL'
        self.fraction_of_optimum = 'NULL'
        self.user_defined = 'NULL'
        self.concordance = 'NULL'
        self.GPR_integration = 'NULL'
        self.percent_of_mapping = 'NULL'


# Read in transcriptomic read abundances, default is tsv with no header 
def read_transcription_file(file, header = False, replicates = False, sep = '\t', 
	binning = False, quant_max = 0.9, quant_min = 0.5, step = 0.125):
    '''Generates dictionary of transcriptomic abundances from a file.
    
    Parameters
    ----------
    file : string
        User-provided file name which contains gene IDs and associated transcription values
    header : boolean
        Defines if read abundance file has a header that needs to be ignored
        default is no header
    replicates : boolean
        Defines if read abundances contains replicates and medians require calculation
        default is no replicates
    sep : string
        Defines what character separates entries on each line
        defaults to tab (.tsv)
    binning : boolean
		Perform discrete binning of transcript abundances into quantiles
		OPTIONAL, not advised
		default is False
	quant_max : float
		Largest quantile to consider
		default is 0.9
	quant_min : float
		Largest quantile to consider
		default is 0.5
	step : float
		Step size for parsing quantiles
		default is 0.125
    '''
    abund_dict = {}
    with open(file, 'r') as transcription:
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

    if binning != False:
    	print('Performing transcript abundance binning by quantile...')
    	abund_dict = _assign_quantiles(abund_dict, quant_max, quant_min, step)

    return abund_dict


# Creates transcription abundances catagories based on quantiles - optional
def _assign_quantiles(transcription, quant_max, quant_min, step):

	if quant_max >= 1.0 or quant_min <= 0.0:
		raise ValueError('ERROR: Quantile range values must be between 1.0 and 0.0! Please correct')
	elif step >= 1.0 or step <= 0.0:
		raise ValueError('ERROR: Quantile step must be between 1.0 and 0.0! Please correct')

	abundances = list(transcription.values())
	abundances.sort()

	thresholds = [max(abundances)]
	while quant_max >= quant_min:
		if quant_max <= 0.0: quant_max = 0.01
		thresholds.append(numpy.quantile(abundances, quant_max))
		quant_max -= step
	thresholds.sort()

	# Identifies which quantile each transcript abundance is in and assigns the largest value in that range instead
	for gene in transcription.keys():
		abund = transcription[gene]

		if abund in thresholds:
			index = thresholds.index(abund)
			transcription[gene] = thresholds[index]
		else:
			index = bisect.bisect_right(thresholds, abund)
			transcription[gene] = thresholds[index]

	return transcription


# Create context-specific model based on transcript distribution
def contextualize(model, transcriptome, samples = 500, norm = True,
    fraction = 0.8, minimum = None, conservative = False, objective = True, 
    set_bounds = True, tasks = [], exclude = [], gpr = False, threshold = 1e-6):

    '''Reaction Inclusion by Parsimony and Transcriptomic Distribution or RIPTiDe
    
    Creates a contextualized metabolic model based on parsimonious usage of reactions defined
    by their associated transcriptomic abundances. Returns a pruned, context-specific cobra.Model 
    and a pandas.DataFrame of associated flux analysis along with parameters used in a riptide class object

    Parameters
    ----------
    model : cobra.Model
        The model to be contextualized
        REQUIRED
    transcriptome : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
        REQUIRED
    samples : int 
        Number of flux samples to collect, default is 500
    norm : bool
        Normalize transcript abundances using RPM calculation
        Performed by default
    fraction : float
        Minimum percent of optimal objective value during FBA steps
        Default is 0.8
    minimum : float
        Minimum linear coefficient allowed during weight calculation for pFBA
        Default is None
    conservative : bool
        Conservatively remove inactive reactions based on genes
        Default is False
    objective : bool
        Sets previous objective function as a constraint with minimum flux equal to user input fraction
        Default is True
    set_bounds : bool
        Uses flax variability analysis results from constrained model to set new bounds for all equations
        Default is False
    tasks : list
        List of reaction ID strings for forced inclusion in final model (metabolic tasks)
    exclude : list
        List of reaction ID strings for forced exclusion from final model
    gpr : bool
        Determines if GPR rules will be considered during coefficient assignment
        Default is False
    threshold : float
        Minimum flux a reaction must acheive in order to avoid pruning during flux sum minimization step
        Default is 1e-6
    '''

    start_time = time.time()
    riptide_object = riptideClass()
    
    # Correct some possible user error
    samples = int(samples)
    if samples <= 0: samples = 1
    fraction = float(fraction)
    if fraction <= 0.0: 
        fraction = 0.01
    elif fraction >= 1.0: 
        fraction = 0.99
    if minimum != None:
        if minimum <= 0.0: 
            minimum = 0.0001
        elif minimum > 1.0: 
            minimum = 0.0001
    solution = model.slim_optimize()
    if model.slim_optimize() < 1e-6 or str(model.slim_optimize()) == 'nan':
        raise ValueError('ERROR: Provided model objective cannot carry flux! Please correct')
    minimum_threshold = threshold

    # Save parameters as part of the output object
    riptide_object.fraction_of_optimum = fraction
    riptide_object.transcriptome = transcriptome
    riptide_object.GPR_integration = gpr

    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    print('\nInitializing model and integrating transcriptomic data...')
    riptide_model = copy.deepcopy(model)
    riptide_model.id = str(riptide_model.id) + '_riptide'
    riptide_object.user_defined = {'tasks':tasks, 'excluded':exclude}

    # Remove totally blocked reactions to speed up subsequent sections
    blocked_rxns = set(find_blocked_reactions(riptide_model))
    blocked_rxns = blocked_rxns.difference(set(tasks))
    blocked_rxns = list(blocked_rxns.union(set(exclude)))
    riptide_model = _prune_model(riptide_model, blocked_rxns, conservative)
    min_coefficient_dict, max_coefficient_dict, gene_hits = _assign_coefficients(transcriptome, riptide_model, minimum, norm, gpr)
    riptide_object.minimization_coefficients = min_coefficient_dict
    riptide_object.maximization_coefficients = max_coefficient_dict
    riptide_object.percent_of_mapping = gene_hits

    # Prune now inactive network sections based on coefficients
    print('Pruning zero flux subnetworks...')
    rm_rxns = _constrain_and_analyze_model(riptide_model, min_coefficient_dict, fraction, 0, objective, tasks, minimum_threshold)
    riptide_model = _prune_model(riptide_model, rm_rxns, conservative)

    # Find optimal solution space based on transcription and final constraints
    print('Analyzing context-specific flux distributions...')
    flux_samples, fva_result, concordance = _constrain_and_analyze_model(riptide_model, max_coefficient_dict, fraction, samples, objective, tasks)
    riptide_object.flux_samples = flux_samples
    riptide_object.flux_variability = fva_result
    riptide_object.concordance = concordance

    # Assign new reaction bounds
    if set_bounds == True:
        riptide_model = _set_new_bounds(riptide_model, fva_result)
    riptide_object.model = riptide_model

    # Analyze changes introduced by RIPTiDe and return results
    _operation_report(start_time, model, riptide_model, concordance)
    return riptide_object


# Converts a dictionary of transcript abundances to reaction linear coefficients
def _assign_coefficients(raw_transcription_dict, model, minimum, norm, gpr):
    
    # Screen transcriptomic abundances for genes that are included in model
    transcription_dict = {}
    total = 0.0
    fail = 0.0
    for gene in model.genes:
        total += 1.0
        try:
            transcription_dict[gene.id] = float(raw_transcription_dict[gene.id])            
        except KeyError:
            fail += 1.0
            continue
    # Check if any genes were found
    if total == fail:
    	raise LookupError('ERROR: No gene IDs in transcriptome dictionary found in model.')
    gene_hits = (float(total - fail) / total) * 100.0
    gene_hits = str(round(gene_hits, 2)) + '%'

    # Perform RPM normalization if specified
    if norm == True:
        total_transcript = float(sum(transcription_dict.values()))
        for gene in transcription_dict.keys():
            new_abund = (transcription_dict[gene] / total_transcript) * 1000000.0
            new_abund = round(new_abund, 3)
            transcription_dict[gene] = new_abund

    # Calculate transcript abundance based coefficients, handle divide-by-zero errors
    direct = False
    abund_distribution = list(set(transcription_dict.values()))
    abund_distribution.sort()
    max_transcript = float(max(abund_distribution)) + 1.0
    if direct == True:
        coefficients = [float(x+1.0) / max_transcript for x in abund_distribution]
        coefficient_range = max(coefficients) + min(coefficients)
        coefficients_rev = [(coefficient_range-x) for x in coefficients]
    else:
        coefficients = [float(x+1.0) / max_transcript for x in abund_distribution]
        coefficients_rev = coefficients.copy()
        coefficients_rev.reverse()

    # Assign coefficients to abundances, adhere to minimum if provided by user
    abund_min_coefficient_dict = {}
    abund_max_coefficient_dict = {}
    for index in range(0, len(abund_distribution)):
        curr_min_coefficient = coefficients_rev[index]
        if minimum != None: 
            if coefficients_rev[index] < minimum: 
                curr_min_coefficient = minimum
        curr_max_coefficient = coefficients[index]
        abund_min_coefficient_dict[abund_distribution[index]] = curr_min_coefficient
        abund_max_coefficient_dict[abund_distribution[index]] = curr_max_coefficient

    # Assign coefficients to reactions
    rxn_min_coefficient_dict = {}
    rxn_max_coefficient_dict = {}
    for gene in transcription_dict.keys():
        transcription = transcription_dict[gene]
        min_coefficient = abund_min_coefficient_dict[transcription]
        max_coefficient = abund_max_coefficient_dict[transcription]

        # Assign corresponding coefficients to reactions associated with each gene
        for rxn in list(model.genes.get_by_any(gene)[0].reactions):            
            if rxn.id in rxn_min_coefficient_dict.keys():
                rxn_min_coefficient_dict[rxn.id].append(min_coefficient)
                rxn_max_coefficient_dict[rxn.id].append(max_coefficient)
            else:
                rxn_min_coefficient_dict[rxn.id] = [min_coefficient]
                rxn_max_coefficient_dict[rxn.id] = [max_coefficient]
    
    # Select final coefficients
    for rxn in model.reactions:
        try:
            # Parse GPRs if defined by user
            if gpr == True:
                curr_gpr = str(rxn.gene_reaction_rule).upper()
                if ' AND ' in curr_gpr:
                    rxn_min_coefficient_dict[rxn.id] = max(rxn_min_coefficient_dict[rxn.id])
                    rxn_max_coefficient_dict[rxn.id] = min(rxn_max_coefficient_dict[rxn.id])
                elif ' OR ' in curr_gpr:
                    rxn_min_coefficient_dict[rxn.id] = min(rxn_min_coefficient_dict[rxn.id])
                    rxn_max_coefficient_dict[rxn.id] = max(rxn_max_coefficient_dict[rxn.id])
                else:
                    rxn_min_coefficient_dict[rxn.id] = numpy.median(rxn_min_coefficient_dict[rxn.id])
                    rxn_max_coefficient_dict[rxn.id] = numpy.median(rxn_max_coefficient_dict[rxn.id])
            else:
                rxn_min_coefficient_dict[rxn.id] = min(rxn_min_coefficient_dict[rxn.id])
                rxn_max_coefficient_dict[rxn.id] = max(rxn_max_coefficient_dict[rxn.id])
        # Coefficient if no gene is associated
        except KeyError:
            rxn_min_coefficient_dict[rxn.id] = numpy.median(coefficients)
            rxn_max_coefficient_dict[rxn.id] = numpy.median(coefficients)
            continue
    
    return rxn_min_coefficient_dict, rxn_max_coefficient_dict, gene_hits


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth, objective, tasks, minimum_threshold=1e-6):
    
    with model as constrained_model:

        # Set previous objective as a constraint, allow deviation
        if objective == True:
            prev_obj_val = constrained_model.slim_optimize()
            prev_obj_constraint = constrained_model.problem.Constraint(constrained_model.objective.expression, lb=prev_obj_val*fraction, ub=prev_obj_val)
            constrained_model.add_cons_vars([prev_obj_constraint])
            constrained_model.solver.update()

        # Include metabolic task constraints
        if len(tasks) >= 1:
            for rxn in tasks:
                with constrained_model as m:
                    task_obj_val = m.slim_optimize()
                    task_constraint = m.problem.Constraint(m.objective.expression, lb=task_obj_val*0.01, ub=task_obj_val)
                    task_constraints.append(task_constraint)

            constrained_model.add_cons_vars(task_constraints)
            constrained_model.solver.update()

        # Apply weigths to new expression
        pfba_expr = symengine.RealDouble(0)
        for rxn in constrained_model.reactions:
            pfba_expr += coefficient_dict[rxn.id] * rxn.forward_variable
            pfba_expr += coefficient_dict[rxn.id] * rxn.reverse_variable

        if sampling_depth == 0:
            # Determine reactions that do not carry any flux in highly constrained solution space
            constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='min', sloppy=True)
            constrained_model.solver.update()
            solution = constrained_model.optimize()
            inactive_rxns = set([rxn.id for rxn in constrained_model.reactions if abs(solution.fluxes[rxn.id]) <= minimum_threshold])
            
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

            # Calculate concordance
            concordance = _calc_concordance(flux_samples, coefficient_dict)

            return flux_samples, fva, concordance


# Find level of concordance between contextualized flux and assigned coefficients
def _calc_concordance(flux_samples, coefficient_dict):
    warnings.filterwarnings('ignore') # Handle RuntimeWarning

    flux_medians = []
    coefficients = []

    for rxn in coefficient_dict.keys():
        try:
            flux_medians.append(abs(numpy.median(list(flux_samples[rxn]))))
            coefficients.append(coefficient_dict[rxn])
        except:
            continue
        
    r_val, p_val = spearmanr(coefficients, flux_medians)
    concordance_dict = {'rho':r_val, 'p':p_val}
    
    warnings.filterwarnings('default')
    return concordance_dict


# Prune model based on blocked reactions from minimization as well as user-defined reactions
def _prune_model(new_model, rm_rxns, conserve):
    
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
    warnings.filterwarnings('ignore') # Handle UserWarning
    for rxn in final_rm_rxns: new_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
    warnings.filterwarnings('default')

    # Prune orphaned nodes
    new_model = _complete_orphan_prune(new_model)
    
    return new_model


# Thoroughly remove orphan reactions and metabolites
def _complete_orphan_prune(model):

    removed = 1
    while removed == 1:
        removed = 0

        for cpd in new_model.metabolites:
            if len(cpd.reactions) == 0:
                cpd.remove_from_model(); removed = 1

        for rxn in new_model.reactions:
            if len(rxn.metabolites) == 0: 
                rxn.remove_from_model(); removed = 1

    return model


# Use flux variability analysis on the constrained model to set new reaction bounds
def _set_new_bounds(model, fva):

    # Set new bounds for all reactions
    for rxn in model.reactions:
        fva_result = list(fva.loc[rxn.id])
        rxn.bounds = (min(fva_result), max(fva_result))

    return model


# Reports how long RIPTiDe took to run
def _operation_report(start_time, model, riptide, concordance):
    
    # Pruning
    perc_removal = 100.0 - ((float(len(riptide.reactions)) / float(len(model.reactions))) * 100.0)
    perc_removal = round(perc_removal, 2)
    print('\nReactions pruned to ' + str(len(riptide.reactions)) + ' from ' + str(len(model.reactions)) + ' (' + str(perc_removal) + '% change)')
    perc_removal = 100.0 - ((float(len(riptide.metabolites)) / float(len(model.metabolites))) * 100.0)
    perc_removal = round(perc_removal, 2)
    print('Metabolites pruned to ' + str(len(riptide.metabolites)) + ' from ' + str(len(model.metabolites)) + ' (' + str(perc_removal) + '% change)')
    
    # Flux through objective
    model_check = 'works'
    # Check that prune model can still achieve flux through the objective (just in case)
    try:
        if riptide.slim_optimize() < 1e-6 or str(riptide.slim_optimize()) == 'nan':
            print('WARNING: Contextualized model objective can no longer carry flux')
            model_check = 'broken'
    except:
        pass

    if model_check == 'works':
        new_ov = round(riptide.slim_optimize(), 2)
        old_ov = round(model.slim_optimize(), 2)
        perc_shift = 100.0 - ((float(new_ov) / float(old_ov)) * 100.0)
        if perc_shift == 0.0:
            pass
        elif perc_shift > 0.0:
            perc_shift = round(abs(perc_shift), 2)
            print('Flux through the objective DECREASED to ~' + str(new_ov) + ' from ~' + str(old_ov) + ' (' + str(perc_shift) + '% change)')
        elif perc_shift < 0.0:
            perc_shift = round(abs(perc_shift), 2)
            print('Flux through the objective INCREASED to ~' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(perc_shift) + '% change)')
    
    # Report concordance
    if concordance['rho'] > 0.0 and concordance['p'] <= 0.05:
        if concordance['p'] < 0.0001:
            p_val = 'p<<0.001 ***'
        elif concordance['p'] <=0.001:
            p_val = 'p=' + str(round(concordance['p'], 3)) + ' ***'
        elif concordance['p'] <= 0.01:
            p_val = 'p=' + str(round(concordance['p'], 3)) + ' **'
        elif concordance['p'] <= 0.05:
            p_val = 'p=' + str(round(concordance['p'], 3)) + ' *'
        print('Context-specific metabolism correlates with transcriptome (' + p_val + ')')

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

# gapsplit flux sampler
# Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
# bioRxiv 652917; doi: https://doi.org/10.1101/652917 

def _gapsplit(model, n=500):

    # Define a few more variables
    max_tries=1000
    primary_tol=0.001
    secondary_frac=0.05
    min_range=1e-5
    enforce_range=True

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
        model.objective = symengine.RealDouble(0)
        solution = model.optimize()
        if solution.status != 'optimal':
            return None
        else:
            return solution.fluxes


def _maxgap(points, fva):
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
