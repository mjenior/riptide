#!/usr/bin/env python3

import os
import sys
import time
import numpy
import bisect
import pandas
import warnings
from random import seed
from copy import deepcopy
from datetime import datetime

import cobra
import symengine
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
        self.fraction_of_optimum = 'NULL'
        self.metabolic_tasks = 'NULL'
        self.concordance = 'NULL'
        self.gpr_integration = 'NULL'
        self.percent_of_mapping = 'NULL'
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

    warnings.filterwarnings('ignore', category=DeprecationWarning) # cobra/io/dict.py np.float 

    if riptide_obj == 'NULL':
        raise ValueError('ERROR: Did not provide a RIPTiDe object')
    
    if path == 'NULL':
        if silent == False:
            print('WARNING: Did not provide an output directory. Using default riptide_files in working directory')
        path = riptide_obj.model.id + '_' + riptide_obj.run_start
    try:
        path = path + '_' + riptide_obj.run_start
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


# Simplified to run in parallel
def _multi_contextualize(args, frac):
    current_fit = contextualize(model=args['model'], transcriptome=args['transcriptome'], fraction=frac, samples=args['samples'], 
                         minimum=args['minimum'], conservative=args['conservative'], objective=args['objective'], 
                         additive=args['additive'], set_bounds=args['set_bounds'], task_frac=args['task_frac'],
                         tasks=args['tasks'], exclude=args['exclude'], gpr=args['gpr'], threshold=args['threshold'], 
                         phase=2, silent=True)

    return current_fit


# Finds the best Rho value from a list of riptide objects
def _find_best_fit(frac_range, argDict, prev_best=None):
    increment = 100.0 / float(len(frac_range)+1)
    progress = 0.0

    if prev_best == None:
        best_fit_concordance = 0.
        if argDict['silent'] == False:
            print('Analyzing context-specific subnetwork flux ranges...')
    else:
        best_fit_concordance = prev_best.concordance
        if argDict['silent'] == False:
            print('Testing local objective fractions to ' + str(best_fit.fraction_of_optimum) + '...')

    if argDict['silent'] == False: sys.stdout.write('\rProgress: 0%')

    maxfit_report = {}
    best_fit = _multi_contextualize(argDict, frac_range[0])
    maxfit_report[best_fit.fraction_of_optimum] = best_fit.concordance
    
    progress += increment
    improved = 0
    if isinstance(best_fit.concordance, str) != True:
        best_fit_concordance = 0.
        if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%          ')
    elif best_fit.concordance['r'] > best_fit_concordance:
        best_fit_concordance = best_fit.concordance
        improved += 1
        if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%  -  improved fit identified (' + str(improved) + ')          ')
    
    # Identify best fit of flux sample to transcriptome
    for frac in frac_range[1:]:
        try:
            fit = _multi_contextualize(argDict, frac)
            improvement = False
        except:
            progress += increment
            if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%          ')
            continue
            
        maxfit_report[fit.fraction_of_optimum] = fit.concordance
        if isinstance(fit.concordance, str) != True:
            if fit.concordance['r'] > best_fit_concordance: 
                improvement = True
                best_fit = fit
                best_fit_concordance = fit.concordance['r']
        
        progress += increment
        if argDict['silent'] == False:
            if improvement == False:
                if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%          ')
            elif improved == 0 and prev_best == None:
                improved += 1
                if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%  -  fit identified          ')
            else:
                improved += 1
                if argDict['silent'] == False: sys.stdout.write('\rProgress: ' + str(float("%.2f" % progress)) + '%  -  improved fit identified (' + str(improved) + ')          ')

    if argDict['silent'] == False:
        if improvement == False and improved == 0:
            if argDict['silent'] == False: 
                sys.stdout.write('\rProgress: 100%          \n\n')
                sys.stdout.flush()
        elif improved == 0 and prev_best == None:
            if argDict['silent'] == False:
                sys.stdout.write('\rProgress: 100%  -  fit identified          \n\n')
                sys.stdout.flush()
        else:
            improved += 1
            if argDict['silent'] == False: 
                sys.stdout.write('\rProgress: 100%  -  improved fit identified (' + str(improved) + ')          \n\n')
                sys.stdout.flush()

    best_fit.maxfit_report = maxfit_report

    return best_fit


# Iteratively run RIPTiDe over a range of objective minimum fractions
def maxfit(model, transcriptome = 'none', frac_min = 0.1, frac_max = 0.9, frac_step = 0.1, prune = True,
    samples = 1000, minimum = False, conservative = False, objective = True, additive = False, 
    set_bounds = True, silent = False, tasks = [], task_frac = 0.01, exclude = [], gpr = False, threshold = 1e-6):

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
        Default is 1000
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
    direct : bool
        Assigns both minimization and maximization step coefficents directly, instead of relying on abundance distribution
        Default is False
    set_bounds : bool
        Uses flux variability analysis results from constrained model to set new bounds for all reactions
        Default is True
    tasks : list
        List of gene or reaction ID strings for forced inclusion in final model (metabolic tasks or essential genes)
    task_frac : float
        Minimum fraction of optimal flux for metabolic task reactions during pruning
        Default is 0.01
    exclude : list
        List of reaction ID strings for forced exclusion from final model
    gpr : bool
        Determines if GPR rules will be considered during coefficient assignment
        Default is False
    threshold : float
        Minimum flux a reaction must acheive in order to avoid pruning during flux sum minimization step
        Default is 1e-5
    '''

    iter_start = time.time()
    curr_run = str(datetime.now()).replace(' ','_').split('.')[0]
    seed(937162211)

    if samples <= 100:
        raise ValueError('ERROR: Iterative RIPTiDe requires >100 flux sampling depth')
    if transcriptome == 'none':
        raise ValueError('ERROR: Iterative RIPTiDe requires an input transcriptome')

    if _test_exchange_space(model, minimum=1e5):
        raise ValueError('ERROR: Solution space is too constrained to analyze by current exchange bounds')

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

    argDict = {'model':model, 'transcriptome':transcriptome, 'silent':silent, 
               'samples':samples, 'minimum':minimum, 'task_frac':task_frac,
               'conservative':conservative, 'objective':objective, 'additive':additive, 
               'set_bounds':set_bounds, 'tasks':tasks, 'exclude':exclude, 
               'gpr':gpr, 'threshold':threshold, 'phase':2}
    
    warnings.filterwarnings('ignore', category=UserWarning)
    top_fit = _find_best_fit(frac_range, argDict)
    all_maxfit = top_fit.maxfit_report
    
    if isinstance(top_fit.concordance, str) != True:
        small_frac_range = []
        while frac_step >= 0.025:
            frac_step = round(frac_step / 2.0, 3)
            small_frac_range += [round(top_fit.fraction_of_optimum - frac_step, 3), round(top_fit.fraction_of_optimum + frac_step, 3)]

        new_fit = _find_best_fit(list(set(small_frac_range)), argDict, top_fit)
        if isinstance(new_fit.concordance, str) != True:
            if new_fit.concordance['r'] > top_fit.concordance['r']:
                top_fit = new_fit
                for frac in new_fit.maxfit_report.keys():
                     all_maxfit[frac] = new_fit.maxfit_report[frac]
    
    warnings.filterwarnings('default', category=UserWarning)
    top_fit.maxfit_report = all_maxfit
    top_fit.fraction_bounds = [frac_min, frac_max]
    top_fit.transcriptome = transcriptome
    
    # Check best fit
    if silent == False:
        if top_fit.concordance == 'Not performed':
            print('No-level of optimal objective flux correlated with transcriptome')
            print(top_fit.fraction_of_optimum, 'of optimum-associated model returned')
        else:
            print('Context-specific metabolism fit with', top_fit.fraction_of_optimum, 'of optimal objective flux')

    # Analyze changes introduced by RIPTiDe and return results
    report_dict = _operation_report(iter_start, model, top_fit.model, top_fit.concordance, silent, phase=1)
    top_fit.additional_parameters['operation'] = report_dict
    top_fit.maxfit = True
    top_fit.run_start = curr_run

    return top_fit


# Assemble a corrected list of metabolic tasks based on user
def _screen_tasks(model, tasks, silent):
    # Check that each task is in the model
    screened_tasks = []
    for x in tasks:
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
            pass

    # Check if any reactions were found in the model that correspond with supplied IDs
    screened_tasks = set(screened_tasks)
    if len(screened_tasks) == 0:
        if silent == False:
            print('WARNING: No reactions found associated with provided task IDs')
    
    return screened_tasks


# Create context-specific model based on transcript distribution
def contextualize(model, transcriptome = 'none', samples = 1000, silent = False, prune = True,
    fraction = 0.8, minimum = False, conservative = False, objective = True, additive = False, direct = False,
    set_bounds = True, tasks = [], task_frac = 0.01, exclude = [], gpr = False, threshold = 1e-6, phase=1):

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
    riptide_object.start_time = str(start_time)
    curr_run = str(datetime.now()).replace(' ','_').split('.')[0]
    riptide_object.run_start = curr_run
    
    seed(937162211)

    if _test_exchange_space(model, minimum=1e5):
        raise ValueError('ERROR: Solution space is too constrained to analyze by current exchange bounds')

    riptide_object = riptideClass()
    riptide_object.additional_parameters = {}
    riptide_object.additional_parameters['threshold'] = threshold
    riptide_object.additional_parameters['silent'] = silent
    riptide_object.additional_parameters['conservative'] = conservative
    riptide_object.additional_parameters['objective'] = objective
    riptide_object.additional_parameters['additive'] = additive
    riptide_object.additional_parameters['set_bounds'] = set_bounds

    # Correct some possible user error
    samples = int(samples)
    if samples < 1: samples = 500
    riptide_object.additional_parameters['samples'] = samples
    fraction = float(fraction)
    if fraction <= 0.0: 
        fraction = 0.01
    elif fraction >= 1.0: 
        fraction = 0.99
    if task_frac <= 0.0: 
        task_frac = 0.01
    elif task_frac >= 1.0: 
        task_frac = 0.01
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
    if phase == 1:
        if silent == False: 
            print('\nInitializing model and integrating transcriptomic data...')
    
    # Define linear coefficients for both steps
    rxn_transcriptome, gene_hits = _transcript_to_reactions(transcriptome, model, gpr, additive)
    rm_rxns = set([x.id for x in model.reactions])
    all_min_coefficient_dict = {}
    all_max_coefficient_dict = {}
    riptide_object.percent_of_mapping = gene_hits
    
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
        min_coefficient_dict, max_coefficient_dict = _assign_coefficients(current_rxn_transcriptome, model, direct)
        for x in min_coefficient_dict.keys():
            try:
                all_min_coefficient_dict[x].append(min_coefficient_dict[x])
                all_max_coefficient_dict[x].append(max_coefficient_dict[x])
            except KeyError:
                all_min_coefficient_dict[x] = [min_coefficient_dict[x]]
                all_max_coefficient_dict[x] = [max_coefficient_dict[x]]
        min_coefficient = min(list(min_coefficient_dict.values()))
        max_coefficient = max(list(max_coefficient_dict.values()))
        
        # Correct task coefficients
        if objective == True: 
            objective = str(model.objective.expression).split()[0].split('*')[-1]
            min_coefficient_dict[objective] = min_coefficient
            max_coefficient_dict[objective] = max_coefficient
            try:
                all_min_coefficient_dict[objective].append(min_coefficient_dict[objective])
                all_max_coefficient_dict[objective].append(max_coefficient_dict[objective])
            except KeyError:
                all_min_coefficient_dict[objective] = [min_coefficient_dict[objective]]
                all_max_coefficient_dict[objective] = [max_coefficient_dict[objective]]
        if len(tasks) >= 1: 
            screened_tasks = _screen_tasks(model, tasks, silent)
            riptide_object.metabolic_tasks = screened_tasks
            for task in screened_tasks:
                min_coefficient_dict[objective] = min_coefficient
                max_coefficient_dict[objective] = max_coefficient
                try:
                    all_min_coefficient_dict[objective].append(min_coefficient_dict[objective])
                    all_max_coefficient_dict[objective].append(max_coefficient_dict[objective])
                except KeyError:
                    all_min_coefficient_dict[objective] = [min_coefficient_dict[objective]]
                    all_max_coefficient_dict[objective] = [max_coefficient_dict[objective]]
        else:
            screened_tasks = []
            
        # Determine active network sections based on coefficients
        active_rxns = _constrain_for_pruning(model=model, min_coefficients=min_coefficient_dict, fraction=fraction,
                                             objective=objective, minimum_flux=minimum_threshold, silent=silent)
        rm_rxns = rm_rxns.difference(active_rxns)
        if len(screened_tasks) >= 1:
            for task in screened_tasks:
                active_rxns = _constrain_for_pruning(model=model, min_coefficients=min_coefficient_dict, fraction=task_frac,
                                                     objective=task, minimum_flux=minimum_threshold, silent=silent)
                rm_rxns = rm_rxns.difference(active_rxns)
        riptide_object.minimization_coefficients = all_min_coefficient_dict
                
    # Prune inactive model components
    rm_rxns = rm_rxns.union(set(exclude))
    if len(rm_rxns) > 0 and prune == True:
        riptide_model = _prune_model(model, rm_rxns, conservative)
        riptide_object.pruned = _record_pruned_elements(model, riptide_model)
    else:
        riptide_model = deepcopy(model)

    # Find optimal solution space based on transcription and final constraints
    if silent == False:
        if phase == 1:
            print('Analyzing context-specific flux distributions...')

    # Find best sampling weights
    for x in all_max_coefficient_dict.keys():
        max_coefficient_dict[x] = max(all_max_coefficient_dict[x])
    riptide_object.maximization_coefficients = max_coefficient_dict

    # Sample weighted flux distributions
    flux_samples, concordance = _constrain_for_sampling(model=riptide_model, max_coefficients=max_coefficient_dict, sampling_depth=samples,
                                                       objective=objective, frac=fraction, minimum_flux=minimum_threshold, silent=silent)
    riptide_object.flux_samples = flux_samples
    riptide_object.concordance = concordance
        
    # Assign new reaction bounds
    if set_bounds == True and isinstance(flux_samples, str) == False:
        for rxn in riptide_model.reactions:
            current_dist = list(flux_samples[rxn.id])

            try:
                old_bounds = rxn.bounds
                rxn.lower_bound = min(current_dist)
                rxn.upper_bound = max(current_dist)
                if riptide_model.slim_optimize(error_value=0.) == 0.:
                    rxn.bounds = old_bounds
            except:
                continue
    riptide_object.model = riptide_model

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
def _assign_coefficients(rxn_transcript_dict, model, direct):
    
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

    return rxn_min_coefficient_dict, rxn_max_coefficient_dict


# Constrain task minimum flux
def _constrain_task(model, task, frac, min_val, silent):
    
    model.objective = task
    model.objective.direction = 'max'
    
    task_name = task + '_constraint'
    task_bound1 = model.slim_optimize(error_value=0.)
    if silent == False and task_bound1 <= min_val:
        print('WARNING:', task, 'minimum flux infeasible')
        task_bounds = [model.reactions.get_by_id(task).lower_bound, model.reactions.get_by_id(task).upper_bound]
    else:
        task_bound2 = task_bound1 * frac
        task_bounds = [task_bound1, task_bound2]
    task_expression = model.reactions.get_by_id(task).flux_expression
    task_constraint = model.problem.Constraint(task_expression, name=task_name, ub=max(task_bounds), lb=min(task_bounds))
    model.add_cons_vars(task_constraint)
    model.solver.update()
    
    return model


# Create weighted expression - default is pFBA
def _weighted_expression(model, coefficients={}):
    
    if len(list(coefficients.keys())) == 0:
        for rxn in model.reactions: 
            coefficients[rxn.id] = 1.0
    
    new_expr = symengine.RealDouble(0)
    for rxn in model.reactions:
        try:
            new_expr += coefficients[rxn.id] * rxn.forward_variable
            new_expr += coefficients[rxn.id] * rxn.reverse_variable
        except KeyError:
            continue
            
    return new_expr


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_for_pruning(model, min_coefficients, fraction, objective, minimum_flux, silent):
    
    constrained_model = deepcopy(model)
    
    # Add objective/task constraints
    if isinstance(objective, str) == True:
        min_coefficient = min(list(min_coefficients.values()))
        min_coefficients[objective] = min_coefficient
        constrained_model = _constrain_task(constrained_model, objective, fraction, minimum_flux, silent)
            
    pruning_expr = _weighted_expression(model=constrained_model, coefficients=min_coefficients)
    constrained_model.objective = model.problem.Objective(pruning_expr, direction='min', sloppy=True)
    constrained_model.solver.update()
    
    constrained_fluxes = constrained_model.optimize().fluxes
    active_rxns = set([rxn.id for rxn in constrained_model.reactions if abs(constrained_fluxes[rxn.id]) >= minimum_flux]) 
    
    return active_rxns


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
def _constrain_for_sampling(model, max_coefficients, sampling_depth, objective, frac, minimum_flux, silent):
    
    constrained_model = deepcopy(model)
    
    # Apply weigths to new expression, allow deviation
    if isinstance(objective, str) == True:
        constrained_model = _constrain_task(constrained_model, objective, frac, minimum_flux, silent)
    sampling_expr = _weighted_expression(model=constrained_model, coefficients=max_coefficients)
    constrained_model.objective = model.problem.Objective(sampling_expr, direction='max', sloppy=True)
    constrained_model.solver.update()
    
    flux_sum = constrained_model.slim_optimize(error_value=0.)
    flux_sum_constraint = constrained_model.problem.Constraint(sampling_expr, lb=flux_sum*0.01, ub=flux_sum)
    model.add_cons_vars(flux_sum_constraint)
    constrained_model.solver.update()
    
    # Analyze flux ranges and calculate concordance
    try:
        flux_samples = _gapsplit(constrained_model, depth=sampling_depth)
        concordance = _calc_concordance(flux_samples, max_coefficients)
    except:
        print('WARNING: Flux sampling not completed with given constraints.')
        flux_samples = 'Not performed'
        concordance = 'Not performed'

    return flux_samples, concordance

    
# Find level of concordance between contextualized flux and assigned coefficients
def _calc_concordance(flux_samples, coefficient_dict):
    warnings.filterwarnings('ignore', category=RuntimeWarning)

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
    
    warnings.filterwarnings('default', category=RuntimeWarning)
    return concordance_dict


# Prune model based on blocked reactions from minimization as well as user-defined reactions
def _prune_model(model, rm_rxns, conserve):
    
    new_model = deepcopy(model)
    if str(new_model.id).strip() == '':
        new_model.id = 'model_riptide'
    else:
        new_model.id = str(new_model.id) + '_riptide'
    
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
    warnings.filterwarnings('ignore', category=UserWarning)
    for rxn in final_rm_rxns: new_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
    warnings.filterwarnings('default', category=UserWarning)

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
        new_ov = round(riptide.slim_optimize(error_value=0.), 4)
        old_ov = round(model.slim_optimize(error_value=0.), 4)
        perc_shift = 100.0 - ((float(new_ov) / float(old_ov)) * 100.0)
        report_dict['obj_change'] = round(perc_shift, 2)
        if perc_shift == 0.0:
            if silent == False:
                print('No change in objective flux of ~' + str(new_ov))
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
        if silent == False:
            print('WARNING: Context-specific model has no objective flux.')

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


# Calculate rough estimate of solution space based on uptake exchange bounds
def _test_exchange_space(model, minimum=1e5):
    
    with model as m:
        curr_medium = model.medium
        lbs = list(curr_medium.values())
        lbs.sort(reverse=True)
        
        space = numpy.prod(lbs[:3])

    if space >= minimum and len(lbs) < 4:
        return True
    else:
        return False


#-----------------------------------------------------------------#

# gapsplit flux sampler
# Adapted from:
# Keaty TC & Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models.
# bioRxiv 652917; doi: https://doi.org/10.1101/652917 

def _gapsplit(model, depth):
    warnings.filterwarnings('ignore') # Handle some infeasible warnings

    fva = flux_variability_analysis(model, list(model.reactions), fraction_of_optimum=0.001, processes=1)

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

    warnings.filterwarnings('default')
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

