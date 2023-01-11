#!/usr/bin/python

import os
import sys
import copy
import time
import numpy
import cobra
import bisect
import pandas
import warnings
import symengine
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from random import seed
from cobra.util import solver
from scipy.stats import spearmanr
from numpy.random import permutation
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
        self.included_important = 'NULL'
        self.additional_parameters = 'NULL'
        self.fraction_bounds = 'NULL'
        self.maxfit = False


# Save the output of RIPTiDe in a newly created directory
def save_output(riptide_obj='NULL', path='NULL', file_type='JSON', silent=False):
    '''Writes RIPTiDe results to files in a new directory
    
    Parameters
    ----------

    REQUIRED
    riptide_obj : RIPTiDe object
        Class object created by riptide.contextualize() or riptide.maxfit()

    OPTIONAL
    path : str
        New directory to write output files
    file_type : str
        Type of output file for RIPTiDe model
        Accepts either sbml or json
        Default is JSON
    silent : bool
        Silences std out 
        Default is False
    '''

    if riptide_obj == 'NULL':
        raise ValueError('ERROR: Did not provide a RIPTiDe object')
    
    if path == 'NULL':
        if silent == False:
            print('WARNING: Did not provide an output directory. Using default riptide_files in working directory')
        path = 'riptide_files'
    try:
        os.mkdir(path)
    except:
        if silent == False:
            print('WARNING: Output path already exists, overwriting previous files')
        pass

    # Write model to file
    if file_type.upper() == 'JSON':
        outFile = path + '/model.json'
        cobra.io.save_json_model(riptide_obj.model, outFile)
    else:
        outFile = path + '/model.sbml'
        cobra.io.write_sbml_model(riptide_obj.model, outFile)

    # Save transcriptome abundances
    outFile = path + '/transcriptome.tsv'
    with open(outFile, 'w') as transcription:
        for gene in riptide_obj.transcriptome.keys():
            if gene == 'replicates':
                continue
            else:
                abund = '\t'.join([str(x) for x in riptide_obj.transcriptome[gene]])
                transcription.write(gene + '\t' + abund + '\n')

    # Write flux samples and FVA to a tsv
    if isinstance(riptide_obj.flux_samples, str) == False:
        outFile = path + '/flux_samples.tsv'
        riptide_obj.flux_samples.to_csv(outFile, sep='\t')
    if isinstance(riptide_obj.flux_variability, str) == False:
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

    # Save some information from all results from maxfit iterations
    if riptide_obj.maxfit == True:
        outFile = path + '/maxfit_iters.tsv'
        with open(outFile, 'w') as maxfit:
            maxfit.write('opt_fraction\tR_value\tp_value\n')
            for index in riptide_obj.maxfit_report.keys():
                line = str(index) + '\t' + str(riptide_obj.maxfit_report[index]['r']) + '\t' + str(riptide_obj.maxfit_report[index]['p']) + '\n'
                maxfit.write(line)

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
             parameters.write('Max fit optimal fraction lower bound: ' + str(riptide_obj.fraction_bounds[0]) + '\n')
             parameters.write('Max fit optimal fraction upper bound: ' + str(riptide_obj.fraction_bounds[1]) + '\n')

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
def read_transcription_file(file, header = False, sep = '\t', rarefy = False, level = 'default',
    binning = False, quant_max = 0.9, quant_min = 0.5, step = 0.125, norm = True, factor = 1e6, silent=False):
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
    sep : string
        Defines what character separates entries on each line
        Defaults to tab (.tsv)
    rarefy : bool
        Rarefies rounded transcript abundances to desired level
        Default is False
    level : int
        Level by which to rarefy samples
        Default is minimum possible total transcript accross replicates
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
    silent : bool
        Silences std out 
        Default is False
    '''

    # Correct some possible user error
    if quant_max >= 1.0 or quant_max <= 0.0: quant_max = 0.99
    if quant_min <= 0.0 or quant_min >= 1.0: quant_min = 0.01
    if step <= 0.0 or step >= 1.0: step = 0.125
    if factor < 1: factor = 1e3

    # Read in file
    abund_dict = {}
    reps = 1
    total_transcript = 0
    min_transcript = 0
    with open(file, 'r') as transcription:
        if header == True: header_line = transcription.readline()

        for line in transcription:
            line = line.strip().split(sep)
            abundances = line[1:]
            if len(abundances) > 0:
                abundances = [float(x) for x in abundances]
                abund_dict[str(line[0])] = abundances
                total_transcript += float(numpy.median(abundances))
                min_transcript += float(min(abundances))
                if reps < len(abundances): reps = len(abundances)

        abund_dict['replicates'] = reps

    # Rarefy abundances
    if level == 'default': 
        level = int(min_transcript)
    else:
        level = int(level)
    if rarefy == True:
        genes = list(abund_dict.keys())
        for x in range(0, reps):
            curr_abund = []
            for y in genes:
                curr_abund.append(abund_dict[y][x])
            curr_abund = _rarefy(curr_abund, level)

            for z in range(0, len(genes)):
                abund_dict[genes[z]][x] = curr_abund[z]

    # Perform reads per million normalization if specified, RPM by default
    if norm == True:
        for gene in abund_dict.keys():
            for x in range(0, reps):
                if gene == 'replicates':
                    continue
                else:
                    abund_dict[gene][x] = (abund_dict[gene][x] / total_transcript) * float(factor)

    # If user-defined, perform abundance binning by quantile
    if binning == True:
        if silent == False:
            print('Performing transcript abundance binning by quantile...')
        abund_dict = _assign_quantiles(abund_dict, quant_max, quant_min, step)

    return abund_dict


# Creates transcription abundances catagories based on quantiles - optional
def _assign_quantiles(transcription, quant_max, quant_min, step):

    if quant_max >= 1.0 or quant_min <= 0.0:
        raise ValueError('ERROR: Quantile range values must be between 1.0 and 0.0! Please correct')
    elif step >= 1.0 or step <= 0.0:
        raise ValueError('ERROR: Quantile step must be between 1.0 and 0.0! Please correct')

    for gene in transcription.keys(): transcription[gene] = float(numpy.median(transcription[gene]))
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
def _rarefy(abunds, n):

    counts = numpy.array([round(x) for x in list(abunds)])
    if counts.sum() <= n: 
        return list(abunds)

    nz = counts.nonzero()[0]
    unpacked = numpy.concatenate([numpy.repeat(numpy.array(i,), counts[i]) for i in nz])
    permuted = permutation(unpacked)[:n]
    result = numpy.zeros(len(counts))
    for p in permuted:
        result[p] += 1

    return list(result)


# Version of riptide.contextualize compatible with multiprocessing
def _iter_riptide(frac, argDict):
    
    iter = contextualize(model=argDict['model'], transcriptome=argDict['transcriptome'], fraction=frac, 
                         silent=argDict['silent'], samples=argDict['samples'], 
                         minimum=argDict['minimum'], conservative=argDict['conservative'], objective=argDict['objective'], 
                         additive=argDict['additive'], important=argDict['important'], set_bounds=argDict['set_bounds'], 
                         tasks=argDict['tasks'], exclude=argDict['exclude'], gpr=argDict['gpr'], threshold=argDict['threshold'], 
                         open_exchanges=argDict['open_exchanges'], phase=argDict['phase'])
    return iter
    

# Finds the best Rho value from a list of riptide objects
def _find_best_fit(frac_range, argDict):
    increment = 100.0 / float(len(frac_range))
    progress = 0.0
    if argDict['silent'] == False:
        sys.stdout.write('\rProgress: 0%     ')
    
    maxfit_report = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=argDict['cpus']) as executor:

        maxfit_jobs = [executor.submit(_iter_riptide, frac, argDict) for frac in frac_range]
        first_iter = 1
        warnings.filterwarnings('ignore', category=UserWarning)
        for job in concurrent.futures.as_completed(maxfit_jobs):
            progress += increment
            try:
                result = job.result()
            except:
                continue

            progress = float("%.3f" % progress)
            if argDict['silent'] == False:
                sys.stdout.write('\rProgress: ' + str(progress) + '%       ')
                sys.stdout.flush()

            maxfit_report[result.fraction_of_optimum] = result.concordance

            if first_iter == 1:
                best_fit = result
                first_iter = 0
            elif result.concordance['r'] > best_fit.concordance['r']:
                best_fit = result

    warnings.filterwarnings('default', category=UserWarning)

    if argDict['silent'] == False:
        sys.stdout.write('\rProgress: 100%         \n\n')
        sys.stdout.flush()
    best_fit.maxfit_report = maxfit_report

    return best_fit


# Iteratively run RIPTiDe over a range of objective minimum fractions
def maxfit(model, transcriptome = 'none', frac_min = 0.25, frac_max = 0.85, frac_step = 0.1, prune = True,
    samples = 500, cpus = 'all', minimum = False, conservative = False, objective = True, additive = False, 
    important = [], set_bounds = True, silent = False, tasks = [], exclude = [], gpr = False, threshold = 1e-5, open_exchanges = False):

    '''
    Iterative RIPTiDe for a range of minimum objective fluxes, returns model with best correlation 
    with flux sampling results to transcriptome abundances.
    
    Parameters
    ----------

    REQUIRED
    model : cobra.Model
        The model to be contextualized
    transcriptome : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
    
    OPTIONAL
    cpus : int
        CPUs number for parallelization
        Default is all available 
    frac_min : float
        Lower bound for range of minimal fractions to test
        Default is 0.25
    frac_max : float
        Upper bound for range of minimal fractions to test
        Default is 0.85
    frac_step : float
        Starting interval size within fraction range
        Default is 0.1
    prune : bool
        Perform pruning step
        Default is True
    samples : int 
        Number of flux samples to collect
        Default is 500
    silent : bool
        Silences std out 
        Default is False
    minimum : float
        Minimum linear coefficient allowed during weight calculation for pFBA
        Default is False
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
    direct : bool
        Assigns both minimization and maximization step coefficents directly, instead of relying on abundance distribution
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
        Default is 1e-5
    open_exchanges : bool
        Sets all exchange reactions bounds to (-1000., 1000)
        Default is False
    skip_fva : bool
        Skip final flux variability analysis step
        Default is False
    '''

    iter_start = time.time()
    seed(937162211)

    if samples <= 100:
        raise ValueError('ERROR: Iterative RIPTiDe requires >100 flux sampling depth')
    if transcriptome == 'none':
        raise ValueError('ERROR: Iterative RIPTiDe requires an input transcriptome')
    if cpus == 'all':
        cpus = os.cpu_count()
    elif isinstance(cpus, int) == False:
        raise ValueError('ERROR: Invalid number of threads provided')

    if frac_min < 0. or frac_min > 1.:
        if silent == False:
            print('WARNING: Improper minimum fraction provided, setting to default value')
        frac_min = 0.25
    elif frac_min > frac_max:
        if silent == False:
            print('WARNING: Improper minimum fraction provided, setting to default value')
        frac_min = 0.25

    if frac_max < 0. or frac_max > 1.:
        print('WARNING: Improper maximum fraction provided, setting to default value')
        frac_max = 0.85
    elif frac_max < frac_min:
        if silent == False:
            print('WARNING: Improper maximum fraction provided, setting to default value')
        frac_max = 0.85

    frac_range = [round(x, 3) for x in list(numpy.arange(frac_min, frac_max+frac_step, frac_step))]
    if len(frac_range) == 1:
        if silent == False:
            print('WARNING: Only a single fraction is possible in the input bounds and fraction')

    if silent == False:
        print('Running max fit RIPTiDe for objective fraction range:', frac_min, 'to', frac_max, '...')

    argDict = {'model':model, 'transcriptome':transcriptome, 'silent':silent, 'cpus':cpus,
               'samples':samples, 'minimum':minimum, 
               'conservative':conservative, 'objective':objective, 'additive':additive, 
               'important':important, 'set_bounds':set_bounds, 'tasks':tasks, 'exclude':exclude, 
               'gpr':gpr, 'threshold':threshold, 'open_exchanges':open_exchanges, 'phase':2}

    if silent == False:
        print('\nParallelizing analysis of context-specific subnetwork flux...')
    
    top_fit = _find_best_fit(frac_range, argDict)
    all_maxfit = top_fit.maxfit_report
    
    if silent == False:
        print('Testing local objective fractions to ' + str(top_fit.fraction_of_optimum) + '...')
    small_frac_range = []
    while frac_step >= 0.025:
        frac_step = round(frac_step / 2.0, 3)
        small_frac_range += [round(top_fit.fraction_of_optimum - frac_step, 3), round(top_fit.fraction_of_optimum + frac_step, 3)]
    
    new_fit = _find_best_fit(list(set(small_frac_range)), argDict)
    if new_fit.concordance['r'] > top_fit.concordance['r']:
        top_fit = new_fit
        for frac in new_fit.maxfit_report.keys():
             all_maxfit[frac] = new_fit.maxfit_report[frac]

    top_fit.maxfit_report = all_maxfit
    top_fit.fraction_bounds = [frac_min, frac_max]
    top_fit.transcriptome = transcriptome
    if silent == False:
        print('Context-specific metabolism fit with', top_fit.fraction_of_optimum, 'of optimal objective flux')

    # Analyze changes introduced by RIPTiDe and return results
    report_dict = _operation_report(iter_start, model, top_fit.model, top_fit.concordance, silent, phase=1)
    top_fit.additional_parameters['operation'] = report_dict
    top_fit.maxfit = True

    return top_fit


# Create context-specific model based on transcript distribution
def contextualize(model, transcriptome = 'none', samples = 500, silent = False, prune = True,
    fraction = 0.8, minimum = False, conservative = False, objective = True, additive = False, important = [], direct = False,
    set_bounds = True, tasks = [], exclude = [], gpr = False, threshold = 1e-5, open_exchanges = False, skip_fva = False, phase=1):

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
    fraction : float
        Minimum objective fraction used during single run setting
        Default is 0.8

    * Other arguments from iterative implementation are carried over
    '''

    start_time = time.time()
    seed(937162211)

    riptide_object = riptideClass()
    riptide_object.additional_parameters = {}
    riptide_object.additional_parameters['threshold'] = threshold
    riptide_object.additional_parameters['silent'] = silent
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
    if minimum != False:
        if minimum <= 0.0: 
            minimum = 0.0001
        elif minimum > 1.0: 
            minimum = 0.0001
    elif minimum == False:
        minimum = None
    riptide_object.additional_parameters['minimum'] = minimum
    solution = model.slim_optimize(error_value=0.)
    if solution < 1e-6:
        raise ValueError('ERROR: Provided model objective cannot carry flux! Please correct')
    minimum_threshold = threshold
    if isinstance(tasks, list) == False: tasks = [tasks]

    # Creates artificial transcriptome to identify most parsimonious patterns of metabolism
    if transcriptome == 'none':
        if silent == False:
            print('WARNING: No transcriptome provided. Analyzing most parsimonious state')
        transcriptome = {}
        transcriptome['replicates'] = 1
        for gene in model.genes: transcriptome[gene.id] = [1.0]
    else:
        # Make sure transcriptome entries are lists
        for x in transcriptome.keys():
            if isinstance(transcriptome[x], list):
                continue
            else:
                transcriptome[x] = [transcriptome[x]]

        transcriptome['replicates'] = int(len(transcriptome[list(transcriptome.keys())[0]]))
        for gene in transcriptome.keys():
            if gene == 'replicates':
                continue
            else:
                transcriptome[gene] = list(transcriptome[gene])

    # Save parameters as part of the output object
    riptide_object.fraction_of_optimum = fraction
    riptide_object.transcriptome = transcriptome
    riptide_object.gpr_integration = gpr

    # Check original model functionality
    # Partition reactions based on transcription percentile intervals, assign corresponding reaction coefficients
    if silent == False: 
        if phase == 1:
            print('\nInitializing model and integrating transcriptomic data...')
    riptide_model = copy.deepcopy(model)
    riptide_model.id = str(riptide_model.id) + '_riptide'
    riptide_object.metabolic_tasks = tasks

    # Open exchange reactions
    if open_exchanges == True:
        for rxn in riptide_model.exchanges: rxn.bounds = (-1000., 1000.)

    # Remove totally blocked reactions to speed up subsequent sections
    rm_rxns = list(set(exclude).difference(set(tasks)))
    if len(rm_rxns) > 0 and prune == True:
        riptide_model = _prune_model(riptide_model, rm_rxns, conservative)

    # Define linear coefficients for both steps
    rxn_transcriptome, gene_hits = _transcript_to_reactions(transcriptome, riptide_model, gpr, additive)
    keep_rxns = set()
    all_min_coefficient_dict = {}
    if silent == False:
        if phase == 1:
            print('Pruning zero flux subnetworks...')
    for x in range(0, transcriptome['replicates']):
        current_rxn_transcriptome = {}
        for index in rxn_transcriptome.keys():
            if index == 'replicates':
                continue
            else:
                current_rxn_transcriptome[index] = rxn_transcriptome[index][x]
        min_coefficient_dict, max_coefficient_dict, important_type = _assign_coefficients(current_rxn_transcriptome, riptide_model, important, direct)
        for x in min_coefficient_dict.keys():
            try:
                all_min_coefficient_dict[x].append(min_coefficient_dict[x])
            except KeyError:
                all_min_coefficient_dict[x] = [min_coefficient_dict[x]]
        
        # Determine active network sections based on coefficients
        active_rxns = _constrain_and_analyze_model(model=riptide_model, coefficient_dict=min_coefficient_dict, fraction=fraction, sampling_depth=0, 
            objective=objective, tasks=tasks, minimum_threshold=minimum_threshold, skip_fva=skip_fva)
        keep_rxns = keep_rxns.union(active_rxns)

    # Determine inactive reactions and prune model
    rm_rxns = set([x.id for x in riptide_model.reactions]).difference(keep_rxns)
    if len(rm_rxns) > 0 and prune == True:
        riptide_model = _prune_model(riptide_model, rm_rxns, conservative)
        riptide_object.pruned = _record_pruned_elements(model, riptide_model)

    riptide_object.minimization_coefficients = all_min_coefficient_dict
    riptide_object.percent_of_mapping = gene_hits

    # Find optimal solution space based on transcription and final constraints
    if silent == False:
        if phase == 1:
            print('Analyzing context-specific flux distributions...')
    median_rxn_transcriptome = {}
    for x in rxn_transcriptome.keys(): 
        if x == 'replicates':
            continue
        else:
            median_rxn_transcriptome[x] = numpy.median(rxn_transcriptome[x])
    min_coefficient_dict, max_coefficient_dict, important_type = _assign_coefficients(median_rxn_transcriptome, riptide_model, important, direct)
    flux_samples, fva_result, concordance = _constrain_and_analyze_model(model=riptide_model, coefficient_dict=max_coefficient_dict, fraction=fraction, 
        sampling_depth=samples, objective=objective, tasks=tasks, minimum_threshold=minimum_threshold, skip_fva=skip_fva)
    riptide_object.maximization_coefficients = max_coefficient_dict
    riptide_object.flux_samples = flux_samples
    riptide_object.flux_variability = fva_result
    riptide_object.concordance = concordance

    # Assign new reaction bounds
    if set_bounds == True and skip_fva == False: 
        for rxn in riptide_model.reactions:
            current_fva = list(fva_result.loc[rxn.id])
            current_lb = min(current_fva)
            current_ub = max(current_fva)
            if current_lb == current_ub:
                current_lb = current_lb * 0.9
                current_ub = current_ub * 1.1
            rxn.bounds = (current_lb, current_ub)
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
    report_dict = _operation_report(start_time, model, riptide_model, concordance, silent, phase)
    riptide_object.additional_parameters['operation'] = report_dict

    return riptide_object
 

# Converts a dictionary of transcript abundances to reaction linear coefficients
def _transcript_to_reactions(reps_transcription_dict, model, gpr, additive):
    
    # Screen transcriptomic abundances for genes that are included in model
    rxn_transcript_dict = {}
    for rxn in model.reactions: rxn_transcript_dict[rxn.id] = []
    total = 0.0
    success = 0.0
    fail = 0.0

    reps = reps_transcription_dict['replicates']
    if type(reps) is list: reps = reps[0]
    for x in range(0, reps):
        current_transcription_dict = {}
        for index in reps_transcription_dict.keys():
            if index == 'replicates':
                continue
            else:
                current_transcription_dict[index] = reps_transcription_dict[index][x]

        current_reaction_dict = {}
        for gene in model.genes:
            total += 1.0
            try:
                current_abund = float(current_transcription_dict[gene.id]) + 1.0 
                current_rxns = list(gene.reactions)
                success += 1.0
                for rxn in current_rxns:
                    try:
                        current_reaction_dict[rxn.id].append(current_abund)
                    except KeyError:
                        current_reaction_dict[rxn.id] = [current_abund]
            except KeyError:
                fail += 1.0
                continue
        # Check if any or very few genes were found
        if total == fail: 
            raise LookupError('ERROR: No gene IDs in transcriptome dictionary found in model.')
        gene_hits = (float(total - fail) / total) * 100.0
        gene_hits = str(round(gene_hits, 2)) + '%'
        nogene_abund = numpy.median(list(set(current_transcription_dict.values())))

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
                transcript_abund = current_reaction_dict[rxn.id]
                # Parse GPRs if defined by user
                if gpr == True:
                    curr_gpr = str(rxn.gene_reaction_rule).upper()
                    if ' AND ' in curr_gpr:
                        current_abund = float(min(transcript_abund))
                        rxn_transcript_dict[rxn.id].append(current_abund)
                        all_abundances |= set([current_abund])
                    elif ' OR ' in curr_gpr:
                        current_abund = float(sum(transcript_abund))
                        rxn_transcript_dict[rxn.id].append(current_abund)
                        all_abundances |= set([current_abund])
                    else:
                        current_abund = float(max(transcript_abund))
                        rxn_transcript_dict[rxn.id].append(current_abund)
                        all_abundances |= set([current_abund])
                else:
                    if additive == True:
                        current_abund = float(sum(transcript_abund))
                        rxn_transcript_dict[rxn.id].append(current_abund)
                        all_abundances |= set([current_abund])
                    else:
                        current_abund = float(max(transcript_abund))
                        rxn_transcript_dict[rxn.id].append(current_abund)
                        all_abundances |= set([current_abund])
            # Coefficient if no gene is associated
            except KeyError:
                rxn_transcript_dict[rxn.id].append(nogene_abund)
    
    return rxn_transcript_dict, gene_hits



# Converts a dictionary of transcript abundances to reaction linear coefficients
def _assign_coefficients(rxn_transcript_dict, model, important, direct):
    
    # Calculate coefficients
    all_abundances = list(rxn_transcript_dict.values())
    denom = max(all_abundances)
    all_abundances.sort()
    max_coefficients = [x / denom for x in all_abundances]
    if direct == True:
        fctr = 1.0 + min(max_coefficients)
        min_coefficients = [fctr - (x / denom) for x in all_abundances]
    else:
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

    # If user has defined important genes/reactions, integrate new weights here
    important_type = 'NULL'
    if len(important) > 0: 
        rxn_min_coefficient_dict, important_type = _integrate_important(model, important, rxn_min_coefficient_dict)

    return rxn_min_coefficient_dict, rxn_max_coefficient_dict, important_type


# Assemble a corrected list of metabolic tasks based on user
def _integrate_tasks(model, tasks):
    # Check that each task is in the model
    screened_tasks = []
    for x in set(tasks):
        # Genes
        try:
            rxns = list(model.genes.get_by_id(x).reactions)
            rxns = [y.id for y in rxns]
            screened_tasks += rxns
        except:
            pass
        # Reactions
        try:
            rxn = model.reactions.get_by_id(x)
            screened_tasks.append(rxn.id)
        except:
            continue

    # Check if any reactions were found in the model that correspond with supplied IDs
    if len(screened_tasks) == 0:
        print('WARNING: No reactions found associated with provided task IDs')

    # Iteratively set each as the objective and find new bounds
    for rxn in screened_tasks:
        model.objective = rxn
        task_obj_val = model.slim_optimize(error_value=0.)
        # Check sign of objective value, just in case
        if task_obj_val > 0.0:
            task_constraint = model.problem.Constraint(model.objective.expression, lb=task_obj_val*0.01, ub=task_obj_val)
            model.add_cons_vars(task_constraint)
            model.solver.update()
        elif task_obj_val < 0.0:
            task_constraint = model.problem.Constraint(model.objective.expression, lb=task_obj_val, ub=task_obj_val*0.01)
            model.add_cons_vars(task_constraint)
            model.solver.update()
    
    return model


# Assign heaviest weight during pruning to user-defined important genes
def _integrate_important(model, important, coefficient_dict):

    #weight = numpy.quantile(list(coefficient_dict.values()), 0.25)
    weight = min(list(coefficient_dict.values()))

    # Check if reactions or genes, get reaction IDs
    rxn_ids = []
    for x in important:
        # Genes
        try:
            rxns = list(model.genes.get_by_id(x).reactions)
            rxn_ids += [y.id for y in rxns]
            important_type = 'genes'
        except:
            pass
        # Reactions
        try:
            rxn = model.reactions.get_by_id(x)
            rxn_ids.append(rxn.id)
            important_type = 'reactions'
        except:
            continue

    # Correct weight in coefficient dictionary
    for rxn in set(rxn_ids): coefficient_dict[rxn] = weight

    return coefficient_dict, important_type


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_and_analyze_model(model, coefficient_dict, fraction, sampling_depth, objective, tasks, minimum_threshold, skip_fva):
    
    constrained_model = copy.deepcopy(model)

    # Set previous objective as a constraint, allow deviation
    if objective == True:
        prev_obj_val = constrained_model.slim_optimize(error_value=0.)
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
        active_rxns = set([rxn.id for rxn in constrained_model.reactions if abs(solution.fluxes[rxn.id]) > minimum_threshold])
            
        return active_rxns
        
    else:
        # Explore solution space of constrained model with flux sampling, allow deviation
        constrained_model.objective = constrained_model.problem.Objective(pfba_expr, direction='max', sloppy=True)
        constrained_model.solver.update()
        flux_sum_obj_val = constrained_model.slim_optimize(error_value=0.)
        flux_sum_constraint = constrained_model.problem.Constraint(pfba_expr, lb=flux_sum_obj_val*fraction, ub=flux_sum_obj_val)
        constrained_model.add_cons_vars([flux_sum_constraint])
        constrained_model.solver.update()
            
        # Analyze flux ranges and calculate concordance
        if sampling_depth != 1:
            warnings.filterwarnings('ignore') # Handle possible infeasible warnings
            flux_samples = _gapsplit(constrained_model, depth=sampling_depth)
            warnings.filterwarnings('default')
            concordance = _calc_concordance(flux_samples, coefficient_dict)
        else:
            flux_samples = 'Not performed'
            concordance = 'Not performed'
        
        if skip_fva == False:
            fva = flux_variability_analysis(constrained_model, fraction_of_optimum=fraction)
        else:
            fva = 'Not performed'

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
def _operation_report(start_time, model, riptide, concordance, silent, phase):
    report_dict = {}

    # Pruning
    perc_removal = 100.0 - ((float(len(riptide.reactions)) / float(len(model.reactions))) * 100.0)
    perc_removal = round(perc_removal, 2)
    if silent == False:
        if phase == 1:
            print('\nReactions pruned to ' + str(len(riptide.reactions)) + ' from ' + str(len(model.reactions)) + ' (' + str(perc_removal) + '% change)')
    report_dict['reactions'] = perc_removal
    perc_removal = 100.0 - ((float(len(riptide.metabolites)) / float(len(model.metabolites))) * 100.0)
    perc_removal = round(perc_removal, 2)
    if silent == False:
        if phase == 1:
            print('Metabolites pruned to ' + str(len(riptide.metabolites)) + ' from ' + str(len(model.metabolites)) + ' (' + str(perc_removal) + '% change)')
    report_dict['metabolites'] = perc_removal

    # Flux through objective
    model_check = 'works'
    # Check that prune model can still achieve flux through the objective (just in case)
    try:
        if riptide.slim_optimize(error_value=0.) < 1e-6:
            model_check = 'broken'
    except:
        pass

    if model_check == 'works':
        new_ov = round(riptide.slim_optimize(error_value=0.), 2)
        old_ov = round(model.slim_optimize(error_value=0.), 2)
        perc_shift = 100.0 - ((float(new_ov) / float(old_ov)) * 100.0)
        report_dict['obj_change'] = round(perc_shift, 2)
        if perc_shift == 0.0:
            pass
        elif perc_shift > 0.0:
            perc_shift = round(abs(perc_shift), 2)
            if silent == False:
                if phase == 1:
                    print('Flux through the objective DECREASED to ~' + str(new_ov) + ' from ~' + str(old_ov) + ' (' + str(perc_shift) + '% change)')
        elif perc_shift < 0.0:
            perc_shift = round(abs(perc_shift), 2)
            if silent == False:
                if phase == 1:
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
            if silent == False:
                if phase == 1:
                    print('Context-specific metabolism correlates with transcriptome (' + rho + ', ' + p_val + ')')
        else:
            if silent == False:
                if phase == 1:
                    print('Context-specific metabolism does not correlate with transcriptome (' + rho + ', n.s.)')

    # Run time
    raw_seconds = int(round(time.time() - start_time))
    report_dict['run_time'] = raw_seconds
    if silent == False:
        if raw_seconds < 60:
            if phase == 1:
                print('\nRIPTiDe completed in ' + str(raw_seconds) + ' seconds\n')
        else:
            minutes, seconds = divmod(raw_seconds, 60)
            mins = 'minute'
            if minutes > 1:
                mins = 'minutes'
            secs = 'second'
            if seconds > 1:
                secs = 'seconds'
            
            if phase == 1:
                print('\nRIPTiDe completed in,', str(minutes), mins, 'and', str(int(seconds)), secs, '\n')
            
    return report_dict

#-----------------------------------------------------------------#

# gapsplit flux sampler
# Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
# bioRxiv 652917; doi: https://doi.org/10.1101/652917 

def _gapsplit(model, depth):
    fva = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.001)

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

