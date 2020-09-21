#!/usr/bin/python

import os
import copy
import time
import numpy
import cobra
import bisect
import pandas
import warnings
import symengine
import itertools
from cobra.util import solver
from scipy.stats import spearmanr
from numpy.random import permutation
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis import flux_variability_analysis


# Create a class to house riptide output data
class riptideClass:
    def __init__(self):
        self.model = 'NULL'
        self.transcriptome = 'NULL'
        self.minimization_coefficients = 'NULL'
        self.maximization_coefficients = 'NULL'
        self.pruned = 'NULL'
        self.flux_samples = 'NULL'
        self.flux_variability = 'NULL'
        self.fraction_of_optimum = 'NULL'
        self.metabolic_tasks = 'NULL'
        self.concordance = 'NULL'
        self.gpr_integration = 'NULL'
        self.percent_of_mapping = 'NULL'
        self.defined_coefficients = 'NULL'
        self.included_important = 'NULL'
        self.additional_parameters = 'NULL'
        self.fraction_bounds = 'NULL'
        self.fraction_step = 'NULL'


# Save the output of RIPTiDe in a newly created directory
def save_output(riptide_obj='NULL', path='NULL', file_type='SBML'):
    '''Writes RIPTiDe results to files in a new directory
    
    Parameters
    ----------

    REQUIRED
    riptide_obj : RIPTiDe object
        Class object creared by riptide.contextualize()

    OPTIONAL
    path : str
        New directory to write output files
    file_type : str
        Type of output file for RIPTiDe model
        Accepts either sbml or json
        Default is SBML
    '''

    if riptide_obj == 'NULL':
        raise ValueError('ERROR: Did not provide a RIPTiDe object')
    
    if path == 'NULL':
        print('WARNING: Did not provide an output directory. Using default riptide_files in working directory')
        path = 'riptide_files'
    try:
        os.mkdir(path)
    except:
        print('WARNING: Output path already exists, overwriting previous files')
        pass

    # Write model to file
    if file_type.upper() == 'JSON':
        outFile = path + '/model.json'
        cobra.io.save_json_model(riptide_obj.model, outFile)
    else:
        outFile = path + '/model.sbml'
        cobra.io.write_sbml_model(riptide_obj.model, outFile)

    # Write flux samples and FVA to a tsv
    if isinstance(riptide_obj.flux_samples, str) == False:
        outFile = path + '/flux_samples.tsv'
        riptide_obj.flux_samples.to_csv(outFile, sep='\t')
    outFile = path + '/flux_variability.tsv'
    riptide_obj.flux_variability.to_csv(outFile, sep='\t')

    # Write coefficient dictionaries to tsvs
    outFile = path + '/flux_minimization_coefficients.tsv'
    with open(outFile, 'w') as min_coefficients:
        min_coefficients.write('reaction\tlinear_coefficient\n')
        for rxn in riptide_obj.minimization_coefficients.keys():
            min_coefficients.write(rxn + '\t' + str(riptide_obj.minimization_coefficients[rxn]) + '\n')
    outFile = path + '/flux_maximization_coefficients.tsv'
    with open(outFile, 'w') as max_coefficients:
        max_coefficients.write('reaction\tlinear_coefficient\n')
        for rxn in riptide_obj.maximization_coefficients.keys():
            max_coefficients.write(rxn + '\t' + str(riptide_obj.maximization_coefficients[rxn]) + '\n')

    # Write pruned model component IDs to a tsv
    outFile = path + '/pruned_components.tsv'
    with open(outFile, 'w') as pruned:
    	line = 'Genes:\t' + '\t'.join(riptide_obj.pruned['genes']) + '\n'
    	line += 'Reactions:\t' + '\t'.join(riptide_obj.pruned['reactions']) + '\n'
    	line += 'Metabolites:\t' + '\t'.join(riptide_obj.pruned['metabolites']) + '\n'
    	pruned.write(line)

    # Assemble parameters and output metrics text file
    outFile = path + '/parameters.txt'
    with open(outFile, 'w') as parameters:
        # Parameters
        parameters.write('Fraction of optimum objective value: ' + str(riptide_obj.fraction_of_optimum) + '\n')
        parameters.write('Percent of genes in mapping found in model: ' + str(riptide_obj.percent_of_mapping) + '\n')
        parameters.write('Minimum flux to avoid pruning: ' + str(riptide_obj.additional_parameters['threshold']) + '\n')

        if riptide_obj.fraction_bounds == 'NULL':
            parameters.write('Run in iterative mode: No\n')
        else:
            parameters.write('Run in iterative mode: Yes\n')
        if riptide_obj.gpr_integration == True:
            parameters.write('Differential weighting by GPR: Yes\n')
        else:
            parameters.write('Differential weighting by GPR: No\n')
        if riptide_obj.additional_parameters['exch_weight'] == True:
            parameters.write('Exchanges weighted as adjacent transporters: Yes\n')
        else:
            parameters.write('Exchanges weighted as adjacent transporters: No\n')
        if riptide_obj.additional_parameters['set_bounds'] == True:
            parameters.write('Set reaction bounds based on FVA results: Yes\n')
        else:
            parameters.write('Set reaction bounds based on FVA results: No\n')
        if riptide_obj.additional_parameters['objective'] == True:
            parameters.write('Defined cellular objective: Yes\n')
        else:
            parameters.write('Defined cellular objective: No\n')
        if riptide_obj.additional_parameters['conservative'] == True:
            parameters.write('Conservative pruning based on GPR: Yes\n')
        else:
            parameters.write('Conservative pruning based on GPR: No\n')
        if riptide_obj.additional_parameters['additive'] == True:
            parameters.write('Pooled transcript based on GPR: Yes\n')
        else:
            parameters.write('Pooled transcript based on GPR: No\n')
        if riptide_obj.additional_parameters['open_exchanges'] == True:
            parameters.write('Exchange reactions switched open: Yes\n')
        else:
            parameters.write('Exchange reactions switched open: No\n')
        if riptide_obj.additional_parameters['silent'] == True:
            parameters.write('Run in silent mode: Yes\n')
        else:
            parameters.write('Run in silent mode: No\n\n')
        if riptide_obj.fraction_bounds != 'NULL':
             parameters.write('Iterative optimal fraction lower bound: ' + str(riptide_obj.fraction_bounds[0]) + '\n')
             parameters.write('Iterative optimal fraction upper bound: ' + str(riptide_obj.fraction_bounds[1]) + '\n')
             parameters.write('Iterative optimal fraction step: ' + str(riptide_obj.fraction_step) + '\n')

        # Results
        parameters.write('Percent pruned reactions: ' + str(riptide_obj.additional_parameters['operation']['reactions']) + '%\n')
        parameters.write('Percent pruned metabolites: ' + str(riptide_obj.additional_parameters['operation']['metabolites']) + '%\n')
        if riptide_obj.additional_parameters['operation']['obj_change'] != 'fails':
            parameters.write('Percent change to objective value: ' + str(riptide_obj.additional_parameters['operation']['obj_change']) + '%\n')
        if riptide_obj.concordance != 'Not performed':
            if str(riptide_obj.concordance['r']) != 'nan':
                r_val = round(riptide_obj.concordance['r'], 4)
                p_val = round(riptide_obj.concordance['p'], 4)
                parameters.write('Correlation between activity and transcriptome: R=' + str(r_val) + ', p-value=' + str(p_val) + '\n')
        parameters.write('RIPTiDe run time: ' + str(riptide_obj.additional_parameters['operation']['run_time']) + ' seconds\n')
        

# Read in transcriptomic read abundances, default is tsv with no header 
def read_transcription_file(file, header = False, replicates = False, sep = '\t', rarefy = False, level = 1e5,
    binning = False, quant_max = 0.9, quant_min = 0.5, step = 0.125, norm = True, factor = 1e6):
    '''Generates dictionary of transcriptomic abundances from a file
    
    Parameters
    ----------

    REQUIRED
    file : string
        User-provided file name which contains gene IDs and associated transcription values

    OPTIONAL
    header : boolean
        Defines if read abundance file has a header that needs to be ignored
        Default is no header
    replicates : boolean
        Defines if read abundances contains replicates and medians require calculation
        Default is no replicates
    sep : string
        Defines what character separates entries on each line
        Defaults to tab (.tsv)
    rarefy : bool
        Rarefies rounded transcript abundances to 90% of the smallest replicate
        Default is False
    level : int
        Level by which to rarefy samples
        Default is 100000
    binning : boolean
        Perform discrete binning of transcript abundances into quantiles
        OPTIONAL, not advised
        Default is False
    quant_max : float
        Largest quantile to consider
        Default is 0.9
    quant_min : float
        Largest quantile to consider
        Default is 0.5
    step : float
        Step size for parsing quantiles
        Default is 0.125
    norm : bool
        Normalize transcript abundances using RPM calculation
        Performed by default
    factor : numeric
        Denominator for read normalization calculation
        Default is 1e6 (RPM)
    '''

    # Correct some possible user error
    if quant_max >= 1.0 or quant_max <= 0.0: quant_max = 0.99
    if quant_min <= 0.0 or quant_min >= 1.0: quant_min = 0.01
    if step <= 0.0 or step >= 1.0: step = 0.125
    if factor < 1: factor = 1e3
    level = int(level)

    # Read in file
    abund_dict = {}
    with open(file, 'r') as transcription:
        if header == True: header_line = transcription.readline()

        for line in transcription:
            line = line.split(sep)
            abundances = [float(x) for x in line[1:]]
            abund_dict[str(line[0])] = abundances
            reps = len(abundances)

    # Rarefy abundances
    if rarefy == True:
        genes = list(abund_dict.keys())
        for x in range(0, reps):
            curr_abund = []
            for y in genes:
                curr_abund.append(abund_dict[y][x])
            curr_abund = _rarefy(curr_abund, level)

            for z in range(0, len(genes)):
                abund_dict[genes[z]][x] = curr_abund[z]

    # Calculate median of replicates
    if replicates == True or reps > 1:
        for gene in abund_dict.keys():
            abund_dict[gene] = float(numpy.median(abund_dict[gene]))

    # Perform normalization if specified
    if norm == True:
        total_transcript = float(sum(abund_dict.values()))
        for gene in abund_dict.keys():
            new_abund = (abund_dict[gene] / total_transcript) * float(factor) # RPM by default
            new_abund = round(new_abund, 3)
            abund_dict[gene] = new_abund

    # If user-defined, perform abundance binning by quantile
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


# Rarefies a list of numbers
def _rarefy(counts, n):

    counts = numpy.array([round(x) for x in list(counts)])
    if counts.sum() <= n: return list(result)

    nz = counts.nonzero()[0]
    unpacked = numpy.concatenate([numpy.repeat(numpy.array(i,), counts[i]) for i in nz])
    permuted = permutation(unpacked)[:n]
    result = numpy.zeros(len(counts))
    for p in permuted:
        result[p] += 1

    return list(result)


# Create context-specific model based on transcript distribution
def contextualize(model, transcriptome = 'none', samples = 500, silent = False, exch_weight = False, processes=None,
    fraction = 0.8, minimum = None, conservative = False, objective = True, additive = False, important = [],
    set_bounds = True, tasks = [], exclude = [], gpr = False, threshold = 1e-6, defined = False, open_exchanges = False):

    '''Reaction Inclusion by Parsimony and Transcriptomic Distribution or RIPTiDe
    
    Creates a contextualized metabolic model based on parsimonious usage of reactions defined
    by their associated transcriptomic abundances. Returns a pruned, context-specific cobra.Model 
    and associated flux analyses along with parameters used in a riptide class object

    Parameters
    ----------

    REQUIRED
    model : cobra.Model
        The model to be contextualized

    OPTIONAL
    transcriptome : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
        With default, an artifical transcriptome is generated where all abundances equal 1.0
    samples : int 
        Number of flux samples to collect
        Default is 500
    silent  : bool
        Silences std out 
        Default is False
    exch_weight : bool
        Weight exchange reactions the same as adjacent transporters
        Default is True
    processes : int
        The number of parallel processes to run for FVA. Optional and if not passed,
        will be set to the number of CPUs found. Necessary to change if
        your trying to run paralell instance of RIPTiDe on the same machine
        Default is none
    fraction : float
        Minimum percent of optimal objective value during FBA steps
        Default is 0.8
    minimum : float
        Minimum linear coefficient allowed during weight calculation for pFBA
        Default is None
    conservative : bool
        Conservatively remove inactive reactions based on GPR rules (all member reactions must be inactive to prune)
        Default is False
    objective : bool
        Sets previous objective function as a constraint with minimum flux equal to user input fraction
        Default is True
    additive : bool
        Pool transcription abundances for reactions with multiple contributing gene products
        Default is False
    important : list
        List of gene or reaction ID strings for which the highest weights are assigned regardless of transcription
        Default is False
    set_bounds : bool
        Uses flux variability analysis results from constrained model to set new bounds for all reactions
        Default is True
    tasks : list
        List of gene or reaction ID strings for forced inclusion in final model (metabolic tasks or essential genes)
    exclude : list
        List of reaction ID strings for forced exclusion from final model
    gpr : bool
        Determines if GPR rules will be considered during coefficient assignment
        Default is False
    threshold : float
        Minimum flux a reaction must acheive in order to avoid pruning during flux sum minimization step
        Default is 1e-6
    defined : False or list
        User defined range of linear coeffients, needs to be defined in a list like [1, 0.5, 0.1, 0.01, 0.001]
        Works best paired with binned abundance catagories from riptide.read_transcription_file()
        Default is False
    open_exchanges : bool
        Sets all exchange reactions bounds to (-1000., 1000)
        Default is False
    '''

    start_time = time.time()
    riptide_object = riptideClass()
    riptide_object.additional_parameters = {}
    riptide_object.additional_parameters['threshold'] = threshold
    riptide_object.additional_parameters['silent'] = silent
    riptide_object.additional_parameters['exch_weight'] = exch_weight
    riptide_object.additional_parameters['conservative'] = conservative
    riptide_object.additional_parameters['objective'] = objective
    riptide_object.additional_parameters['additive'] = additive
    riptide_object.additional_parameters['set_bounds'] = set_bounds
    riptide_object.additional_parameters['open_exchanges'] = open_exchanges

    # Correct some possible user error
    samples = int(samples)
    if samples < 1: samples = 500
    if len(important) > 0: samples = 1
    riptide_object.additional_parameters['samples'] = samples
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
    riptide_object.additional_parameters['minimum'] = minimum
    solution = model.slim_optimize()
    if model.slim_optimize() < 1e-6 or str(model.slim_optimize()) == 'nan':
        raise ValueError('ERROR: Provided model objective cannot carry flux! Please correct')
    minimum_threshold = threshold
    if isinstance(tasks, list) == False: tasks = [tasks]
    if defined != False:
        if isinstance(defined, list) == False:
            defined = [defined]

    # Creates artificial transcriptome to identify most parsimonious patterns of metabolism
    if transcriptome == 'none':
        if silent == False:
            print('WARNING: No transcriptome provided. Analyzing most parsimonious state')
        transcriptome = {}
        for gene in model.genes:
            transcriptome[gene.id] = 1.0

    # Save parameters as part of the output object
    riptide_object.fraction_of_optimum = fraction
    riptide_object.transcriptome = transcriptome
    riptide_object.gpr_integration = gpr
    riptide_object.defined_coefficients = defined

    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    if silent == False: print('\nInitializing model and integrating transcriptomic data...')
    riptide_model = copy.deepcopy(model)
    riptide_model.id = str(riptide_model.id) + '_riptide'
    riptide_object.metabolic_tasks = tasks

    # Open exchange reactions
    if open_exchanges == True:
        for rxn in riptide_model.exchanges: rxn.bounds = (-1000., 1000.)

    # Remove totally blocked reactions to speed up subsequent sections
    rm_rxns = list(set(exclude).difference(set(tasks)))
    if len(rm_rxns) > 0:
        riptide_model = _prune_model(riptide_model, rm_rxns, conservative)

    # Define linear coefficients for both steps
    min_coefficient_dict, max_coefficient_dict, gene_hits, important_type = _assign_coefficients(transcriptome, riptide_model, minimum, gpr, defined, additive, exch_weight, important)
    riptide_object.minimization_coefficients = min_coefficient_dict
    riptide_object.maximization_coefficients = max_coefficient_dict
    riptide_object.percent_of_mapping = gene_hits
    
    # Prune now inactive network sections based on coefficients
    if silent == False: print('Pruning zero flux subnetworks...')
    rm_rxns = _constrain_and_analyze_model(riptide_model, min_coefficient_dict, fraction, 0, objective, tasks, minimum_threshold, cpus=processes)
    riptide_model = _prune_model(riptide_model, rm_rxns, conservative)
    riptide_object.pruned = _record_pruned_elements(model, riptide_model)

    # Find optimal solution space based on transcription and final constraints
    if silent == False: print('Analyzing context-specific flux distributions...')
    flux_samples, fva_result, concordance = _constrain_and_analyze_model(riptide_model, max_coefficient_dict, fraction, samples, objective, tasks, cpus=processes)
    riptide_object.flux_samples = flux_samples
    riptide_object.flux_variability = fva_result
    riptide_object.concordance = concordance

    # Assign new reaction bounds
    if set_bounds == True: 
        for rxn in riptide_model.reactions:
            current_fva = list(fva_result.loc[rxn.id])
            rxn.bounds = (min(current_fva), max(current_fva))
    riptide_object.model = riptide_model

    # Report on included subset of user-defined important genes
    if len(important) > 0:
        if important_type == 'genes':
            curr_list = set([x.id for x in riptide_model.genes])
        elif important_type == 'reactions':
            curr_list = set([y.id for y in riptide_model.reactions])
        included_important = set(important).intersection(curr_list)
        riptide_object.included_important = included_important
    else:
        riptide_object.included_important = important

    # Analyze changes introduced by RIPTiDe and return results
    if silent == False: 
        report_dict = _operation_report(start_time, model, riptide_model, concordance)
        riptide_object.additional_parameters['operation'] = report_dict
    else:
        if model.slim_optimize() < 1e-6 or str(model.slim_optimize()) == 'nan':
            print('WARNING: Contextualized model objective can no longer carry flux')
            
    return riptide_object


# Converts a dictionary of transcript abundances to reaction linear coefficients
def _assign_coefficients(raw_transcription_dict, model, minimum, gpr, defined_coefficients, additive, exch_weight, important):
    
    # Screen transcriptomic abundances for genes that are included in model
    rxn_transcript_dict = {}
    total = 0.0
    success = 0.0
    fail = 0.0
    for gene in model.genes:
        total += 1.0
        try:
            current_abund = float(raw_transcription_dict[gene.id]) + 1.0
            current_rxns = list(model.genes.get_by_id(gene.id).reactions)
            success += 1.0
            for rxn in current_rxns:
                try:
                    rxn_transcript_dict[rxn.id].append(current_abund)
                except KeyError:
                    rxn_transcript_dict[rxn.id] = [current_abund]
        except KeyError:
            fail += 1.0
            continue
    # Check if any or very few genes were found
    if total == fail: 
        raise LookupError('ERROR: No gene IDs in transcriptome dictionary found in model.')
    elif success < (len(model.genes) * 0.5):
        print('WARNING: Fewer than half of model genes were found in transcriptome mapping file.')
    gene_hits = (float(total - fail) / total) * 100.0
    gene_hits = str(round(gene_hits, 2)) + '%'
    nogene_abund = numpy.median(list(set(raw_transcription_dict.values())))

    # Reduce transcription to single values per reaction based on max/sum or GPR rules
    exchanges = {}
    all_abundances = set([nogene_abund])
    for rxn in model.reactions:
        # Identify extracellular exchanges
        if set([x.compartment for x in rxn.reactants]) != set([x.compartment for x in rxn.products]):
            for cpd in rxn.metabolites: # Transport reactions
                for sub_rxn in cpd.reactions:
                    if len(sub_rxn.reactants) == 0 or len(sub_rxn.products) == 0:
                        exchanges[rxn.id] = sub_rxn.id
        try:
            transcript_abund = rxn_transcript_dict[rxn.id]
            # Parse GPRs if defined by user
            if gpr == True:
                curr_gpr = str(rxn.gene_reaction_rule).upper()
                if ' AND ' in curr_gpr:
                    current_abund = float(min(transcript_abund))
                    rxn_transcript_dict[rxn.id] = current_abund
                    all_abundances |= set([current_abund])
                elif ' OR ' in curr_gpr:
                    current_abund = float(sum(transcript_abund))
                    rxn_transcript_dict[rxn.id] = current_abund
                    all_abundances |= set([current_abund])
                else:
                    current_abund = float(max(transcript_abund))
                    rxn_transcript_dict[rxn.id] = current_abund
                    all_abundances |= set([current_abund])
            else:
                if additive == True:
                    current_abund = float(sum(transcript_abund))
                    rxn_transcript_dict[rxn.id] = current_abund
                    all_abundances |= set([current_abund])
                else:
                    current_abund = float(max(transcript_abund))
                    rxn_transcript_dict[rxn.id] = current_abund
                    all_abundances |= set([current_abund])
        # Coefficient if no gene is associated
        except KeyError:
            rxn_transcript_dict[rxn.id] = nogene_abund

    # Calculate coefficients
    all_abundances = list(all_abundances)
    all_abundances.sort()
    max_coefficients = [x / max(all_abundances) for x in all_abundances]
    min_coefficients = max_coefficients[::-1]
    coefficient_dict = {}
    for x in range(0, len(all_abundances)):
        coefficient_dict[all_abundances[x]] = [min_coefficients[x], max_coefficients[x]]

    # Assign coefficients to reactions
    rxn_min_coefficient_dict = {}
    rxn_max_coefficient_dict = {}
    for rxn in rxn_transcript_dict.keys():
        rxn_min_coefficient_dict[rxn] = coefficient_dict[rxn_transcript_dict[rxn]][0]
        rxn_max_coefficient_dict[rxn] = coefficient_dict[rxn_transcript_dict[rxn]][1]

    # Set adjacent exchange reactions to the same coefficients
    if exch_weight == True:
        for rxn in exchanges.keys():
            rxn_max_coefficient_dict[exchanges[rxn]] = rxn_max_coefficient_dict[rxn]
            rxn_min_coefficient_dict[exchanges[rxn]] = rxn_min_coefficient_dict[rxn]

    # If user has defined important genes/reactions, integrate new weights here
    important_type = 'NULL'
    if len(important) > 0: 
        rxn_min_coefficient_dict, important_type = _integrate_important(model, important, rxn_min_coefficient_dict)

    return rxn_min_coefficient_dict, rxn_max_coefficient_dict, gene_hits, important_type


# Assemble a corrected list of metabolic tasks based on user
def _integrate_tasks(model, tasks):
    # Check that each task is in the model
    screened_tasks = set()
    for x in tasks:
        # Genes
        try:
            rxns = list(model.genes.get_by_id(x).reactions)
            rxns = set([y.id for y in rxns])
            screened_tasks |= rxns
        except:
            pass
        # Reactions
        try:
            rxn = model.reactions.get_by_id(x)
            screened_tasks |= set([x])
        except:
            continue

    # Iteratively set each as the objective and find new bounds
    task_constraints = []
    for rxn in screened_tasks:
        model.objective = rxn
        task_obj_val = model.slim_optimize()
        # Check sign of objective value, just in case
        if task_obj_val > 0.0:
            task_constraint = model.problem.Constraint(model.objective.expression, lb=task_obj_val*0.01, ub=task_obj_val)
        elif task_obj_val < 0.0:
            task_constraint = model.problem.Constraint(model.objective.expression, lb=task_obj_val, ub=task_obj_val*0.01)
        task_constraints.append(task_constraint)
    
    # Check if any reactions were found in the model that correspond with supplied IDs
    if len(task_constraints) == 0:
        print('WARNING: No reactions found associated with provided task IDs')
    else:
        model.add_cons_vars(task_constraints)
        model.solver.update()

    return model


# Assign heaviest weight during pruning to user-defined important genes
def _integrate_important(model, important, coefficient_dict):

    #weight = numpy.quantile(list(coefficient_dict.values()), 0.25)
    weight = min(list(coefficient_dict.values()))

    # Check if reactions or genes, get reaction IDs
    rxn_ids = set()
    for x in important:
        # Genes
        try:
            rxns = list(model.genes.get_by_id(x).reactions)
            rxn_ids |= set([y.id for y in rxns])
            important_type = 'genes'
        except:
            pass
        # Reactions
        try:
            rxn = model.reactions.get_by_id(x)
            rxn_ids |= set([x])
            important_type = 'reactions'
        except:
            continue

    # Correct weight in coefficient dictionary
    for rxn in rxn_ids: coefficient_dict[rxn] = weight

    return coefficient_dict, important_type


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth, objective, tasks, minimum_threshold=1e-6, cpus):
    
    constrained_model = copy.deepcopy(model)

    # Set previous objective as a constraint, allow deviation
    if objective == True:
        prev_obj_val = constrained_model.slim_optimize()
        prev_obj_constraint = constrained_model.problem.Constraint(constrained_model.objective.expression, lb=prev_obj_val*fraction, ub=prev_obj_val)
        constrained_model.add_cons_vars([prev_obj_constraint])
        constrained_model.solver.update()

    # Apply weigths to new expression
    pfba_expr = symengine.RealDouble(0)
    for rxn in constrained_model.reactions:
        try:
            pfba_expr += coefficient_dict[rxn.id] * rxn.forward_variable
            pfba_expr += coefficient_dict[rxn.id] * rxn.reverse_variable
        except KeyError:
            continue

    if sampling_depth == 0:
        # Include metabolic task constraints
        if len(tasks) >= 1:
            constrained_model = _integrate_tasks(constrained_model, tasks)

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
            
        # Analyze flux ranges and calculate concordance
        if sampling_depth != 1:
            warnings.filterwarnings('ignore') # Handle possible infeasible warnings
            flux_samples = _gapsplit(constrained_model, depth=sampling_depth, cpus=cpus)
            warnings.filterwarnings('default')
            concordance = _calc_concordance(flux_samples, coefficient_dict)
        else:
            flux_samples = 'Not performed'
            concordance = 'Not performed'
            
        fva = flux_variability_analysis(constrained_model, fraction_of_optimum=fraction, processes=cpus)

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
    concordance_dict = {'r':r_val, 'p':p_val}
    
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


# Record IDs of all pruned elements
def _record_pruned_elements(old_model, new_model):
    prunedDict = {}

    # Genes
    old_genes = set([x.id for x in old_model.genes])
    new_genes = set([y.id for y in new_model.genes])
    prunedDict['genes'] = old_genes.difference(new_genes)

    # Reactions
    old_rxns = set([x.id for x in old_model.reactions])
    new_rxns = set([y.id for y in new_model.reactions])
    prunedDict['reactions'] = old_rxns.difference(new_rxns)

    # Metabolites
    old_cpds = set([x.id for x in old_model.metabolites])
    new_cpds = set([y.id for y in new_model.metabolites])
    prunedDict['metabolites'] = old_cpds.difference(new_cpds)

    return prunedDict


# Thoroughly remove orphan reactions and metabolites
def _complete_orphan_prune(model):

    removed = 1
    while removed == 1:
        removed = 0

        # Metabolites
        for cpd in model.metabolites:
            if len(cpd.reactions) == 0:
                cpd.remove_from_model(); removed = 1

        # Reactions
        for rxn in model.reactions:
            if len(rxn.metabolites) == 0: 
                rxn.remove_from_model(); removed = 1

    return model


# Reports how long RIPTiDe took to run
def _operation_report(start_time, model, riptide, concordance):
    report_dict = {}

    # Pruning
    perc_removal = 100.0 - ((float(len(riptide.reactions)) / float(len(model.reactions))) * 100.0)
    perc_removal = round(perc_removal, 2)
    print('\nReactions pruned to ' + str(len(riptide.reactions)) + ' from ' + str(len(model.reactions)) + ' (' + str(perc_removal) + '% change)')
    report_dict['reactions'] = perc_removal
    perc_removal = 100.0 - ((float(len(riptide.metabolites)) / float(len(model.metabolites))) * 100.0)
    perc_removal = round(perc_removal, 2)
    print('Metabolites pruned to ' + str(len(riptide.metabolites)) + ' from ' + str(len(model.metabolites)) + ' (' + str(perc_removal) + '% change)')
    report_dict['metabolites'] = perc_removal

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
        report_dict['obj_change'] = round(perc_shift, 2)
        if perc_shift == 0.0:
            pass
        elif perc_shift > 0.0:
            perc_shift = round(abs(perc_shift), 2)
            print('Flux through the objective DECREASED to ~' + str(new_ov) + ' from ~' + str(old_ov) + ' (' + str(perc_shift) + '% change)')
        elif perc_shift < 0.0:
            perc_shift = round(abs(perc_shift), 2)
            print('Flux through the objective INCREASED to ~' + str(new_ov) + ' from ' + str(old_ov) + ' (' + str(perc_shift) + '% change)')
    else:
        report_dict['obj_change'] = 'fails'

    # Report concordance
    if concordance != 'Not performed':
        rho = 'r=' + str(round(concordance['r'], 3))
        if concordance['p'] <= 0.05:
            p_val = round(concordance['p'], 3)
            if p_val == 0.0:
                p_val = 'p<0.001 *'
            else:
                p_val = 'p=' + str(p_val) + ' *'
            print('Context-specific metabolism correlates with transcriptome (' + rho + ', ' + p_val + ')')
        else:
            print('Context-specific metabolism does not correlate with transcriptome (' + rho + ', n.s.)')

    # Run time
    raw_seconds = int(round(time.time() - start_time))
    report_dict['run_time'] = raw_seconds
    if raw_seconds < 60:
        print('\nRIPTiDe completed in ' + str(raw_seconds) + ' seconds\n')
    else:
        minutes, seconds = divmod(raw_seconds, 60)
        print('\nRIPTiDe completed in ' + str(minutes) + ' minutes and ' + str(int(seconds)) + ' seconds\n')

    return report_dict

#-----------------------------------------------------------------#

# gapsplit flux sampler
# Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
# bioRxiv 652917; doi: https://doi.org/10.1101/652917 

def _gapsplit(model, depth, cpus):
    fva = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.001, processes=cpus)

    # only split reactions with feasible range >= min_range
    idxs = (fva.maximum - fva.minimum >= 1e-5).to_numpy().nonzero()[0]
    weights = (1.0 / (fva.maximum - fva.minimum) ** 2).to_numpy()
    samples = numpy.zeros((depth, len(model.reactions)))
    k = 0
    for try_ in range(1000):
        relative, target, width = _maxgap(samples[0:k,idxs], fva.iloc[idxs,:])
        primary_var = numpy.argmax(relative)
        primary_target = target[primary_var]
        primary_lb = primary_target - 0.001*width[primary_var]
        primary_ub = primary_target + 0.001*width[primary_var]

        new_sample = _generate_sample(model, idxs[primary_var], primary_lb, primary_ub)
        if new_sample is not None:
            new_sample[new_sample > fva.maximum] = fva.maximum[new_sample > fva.maximum]
            new_sample[new_sample < fva.minimum] = fva.minimum[new_sample < fva.minimum]
            samples[k,:] = new_sample
            k += 1
        if k >= depth: break

    if k < depth:
        # max_tries reached; return fewer samples
        samples = samples[:k,:]

    return pandas.DataFrame(data=samples, columns=fva.maximum.index)

def _generate_sample(model, primary_var, primary_lb, primary_ub):
    
    # Formulate a [MI]QP to find a single solution
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


#-----------------------------------------------------------------#


# Iteratively run RIPTiDe over a range of objective minimum fractions
def iterative(model, transcriptome = 'none', frac_min = 0.65, frac_max = 0.85, frac_step = 0.02, samples = 500, exch_weight = False, 
    processes = None, minimum = None, conservative = False, objective = True, additive = False, important = [], set_bounds = True, 
    tasks = [], exclude = [], gpr = False, threshold = 1e-6, defined = False, open_exchanges = False):

    '''Iterative RIPTiDe for a range of minimum objective fluxes
    
    Parameters
    ----------

    REQUIRED
    model : cobra.Model
        The model to be contextualized
    transcriptome : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
    samples : int 
        Number of flux samples to collect
        Default is 500
    frac_min : float
        Lower bound for range of minimal fractions to test
        Default is 0.65
    frac_max : float
        Upper bound for range of minimal fractions to test
        Default is 0.85
    frac_step : float
        Increment to parse input minimal fraction range
        Default is 0.02

    OPTIONAL
    All other optional parameters for riptide.contextualize()
    '''

    iter_start = time.time()

    if samples <= 100:
        raise ValueError('ERROR: Iterative RIPTiDe requires >100 flux sampling depth')
    if transcriptome == 'none':
        raise ValueError('ERROR: Iterative RIPTiDe requires an input transcriptome')

    if frac_min < 0. or frac_min > 1.:
        print('WARNING: Improper minimum fraction provided, setting to default value')
        frac_min = 0.65
    elif frac_min > frac_max:
        print('WARNING: Improper minimum fraction provided, setting to default value')
        frac_min = 0.65

    if frac_max < 0. or frac_max > 1.:
        print('WARNING: Improper maximum fraction provided, setting to default value')
        frac_max = 0.85
    elif frac_max < frac_min:
        print('WARNING: Improper maximum fraction provided, setting to default value')
        frac_max = 0.85

    if frac_step < 0. or frac_min > 1.:
        print('WARNING: Improper fraction step provided, setting to default value')
        frac_max = 0.02
    frac_range = [x for x in range(frac_min, frac_max, frac_step)]
    if len(frac_range) == 1:
        print('WARNING: Only a single fraction is possible in the input bounds and fraction')

    top_rho = 0.
    iters = 0
    for frac in frac_range:
        iters += 1

        iter_riptide = contextualize(model, transcriptome, fraction=frac, silent=True, samples=samples, exch_weight=exch_weight, 
            processes=processes, minimum=minimum, conservative=conservative, 
            objective=objective, additive=additive, important=important, set_bounds=set_bounds, tasks=tasks, exclude=exclude, 
            gpr=gpr, threshold=threshold, defined=defined, open_exchanges=open_exchanges)

        curr_rho = iter_riptide.concordance['r']
        print('Iteration', iters, 'rho:', round(curr_rho, 3))
        if curr_rho > top_rho:
            top_fit = copy.deepcopy(iter_riptide)

    top_fit.fraction_bounds = [frac_min, frac_max]
    top_fit.fraction_step = frac_step

    print('\nBest fit with', top_fit.fraction_of_optimum, 'of optimal objective flux')
    raw_seconds = int(round(time.time() - iter_start))
    report_dict['run_time'] = raw_seconds
    if raw_seconds < 60:
        print('\nIterative RIPTiDe completed in ' + str(raw_seconds) + ' seconds\n')
    else:
        minutes, seconds = divmod(raw_seconds, 60)
        print('\nIterative RIPTiDe completed in ' + str(minutes) + ' minutes and ' + str(int(seconds)) + ' seconds\n')

    return top_fit

