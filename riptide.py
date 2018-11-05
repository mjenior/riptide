#!/usr/bin/python
'''
Reaction
Inclusion by
Parsimonious usage,
Transcript
Distribution, and
Exploration of topology

or
RIPTiDE

Please cite when using:
Jenior ML and Papin JA. (2018) RIPTiDE: Metabolic model contextualization through exploration 
of mutliple transcription and flux minimization levels. BioRxiv. DOI

Example usage 1:
my_model = cobra.io.read_sbml_model('my_model.sbml')
transcript_abundances = read_transcription_file('reads_to_genes.tsv', replicates=True)
reaction_bins = create_reaction_partitions(my_model, transcript_abundances)
contextualized_model = contextualize(my_model, reaction_bins, include_rxns=rxn1-rxn2, remove_rxns=rxn3-rxn4)

Example usage 2:
my_model = cobra.io.read_sbml_model('my_model.sbml')
gene_bins = read_binning_file('gene_bins.tsv')
reaction_bins = create_reaction_partitions(my_model, gene_bins)
contextualized_model = contextualize(my_model, reaction_bins)
'''

# Dependencies
import numpy
import copy
import cobra
from cobra.util import solver
from optlang.symbolics import Zero
from itertools import chain
from cobra.manipulation.delete import *

# Read in transcriptomic read abundances, default is tsv with no header 
def read_transcription_file(read_abundances_file, replicates=False, header=False, sep='\t'):
    
    abund_dict = {}
    with open(read_abundances_file, 'r') as transcription:
        
        if header == True:
            header_line = transcription.readline()

        for line in transcription:
            line = line.split(sep)

            gene = str(line[0])
            if replicates == True:
                abundance = float(numpy.median([float(x) for x in line[1:]]))
            else:
                abundance = float(line[1])

            abund_dict[gene] = abundance

    return abund_dict


# Converts a dictionary of transcript distribution percentiles (defaults: 95th, 75th, & 50th)
# Based on:
# Schultz, A, & Qutub, AA (2016). Reconstruction of Tissue-Specific Metabolic Networks Using CORDA. 
#      PLoS Computational Biology. https://doi.org/10.1371/journal.pcbi.1004808
def create_reaction_partitions(model, transcription, save_bins_as=False, high_cutoff=95, mid_cutoff=75, low_cutoff=50):

    # Define transcript abundance cutoffs
    distribution = transcription.values()
    perc_hi = numpy.percentile(distribution, high_cutoff)
    perc_mid = numpy.percentile(distribution, mid_cutoff)
    perc_lo = numpy.percentile(distribution, low_cutoff)
        
    # Write to file if defined by user
    if save_bins_as != False:
        perc_file = open(save_bins_as, 'w')

    # Keep a count of genes included in each partition
    count_4 = 0
    count_3 = 0
    count_2 = 0
    count_1 = 0
    count_0 = 0
    
    percentile_dict = {}
    for gene in model.genes:
        
        # If current gene not in transciptomic profile, default to inclusion
        failed = 0
        try:
            test = transcription[gene.id]
        except KeyError:
            curr_percentile = 0
            failed = 1
            count_0 += 1
            pass

        # Assign percentile grouping scores
        if failed == 0:
            if transcription[gene.id] >= perc_hi:
                curr_percentile = 1
                count_1 += 1
            elif transcription[gene.id] < perc_hi and transcription[gene.id] >= perc_mid:
                curr_percentile = 2
                count_2 += 1
            elif transcription[gene.id] < perc_mid and transcription[gene.id] >= perc_lo:
                curr_percentile = 3
                count_3 += 1
            elif transcription[gene.id] < perc_lo and transcription[gene.id] >= 0.0:
                curr_percentile = 4
                count_4 += 1

        # Write the binning results to a file if requested by the user
        if save_bins_as != False:
            entry = str(gene.id) + '\t' + str(curr_percentile) + '\n'
            perc_file.write(entry)

        # Converts genes into corresponding reactions and adds them to a dictionary
        for rxn in list(gene.reactions): 
            percentile_dict[rxn.id] = curr_percentile

    if save_bins_as != False:
        perc_file.close()
    
    # Report sizes of each bin
    print('Bin 0 (not found): ' + str(count_0))
    print('Bin 1 (highest): ' + str(count_1))
    print('Bin 2: ' + str(count_2))
    print('Bin 3: ' + str(count_3))
    print('Bin 4 (lowest): ' + str(count_4))
    
    return percentile_dict


# Read in user-defined/editted gene priority binning file
def read_binning_file(partition_file):

    bin_dict = {}
    with open(partition_file, 'r') as percentiles:
        for line in percentiles:
            line = line.split()
            bin_dict[line[0]] = float(line[1])

    return bin_dict


# Bin reactions based on their percentile transcription
def parse_reaction_binning(model, percentiles, keep_rxns=[], remove_rxns=[]):

    include = set()
    exclude = set()
    perc_top = set()
    perc_hi = set()
    perc_mid = set()
    perc_lo = set()
    
    # Assign bins associated with each percentile of the read abundance distribution
    for rxn in model.reactions:

        # Screen lists based on user-defined reactions
        if len(keep_rxns) > 0:
            if rxn.id in keep_rxns:
                include |= set([rxn.id])
                continue
        elif len(remove_rxns) > 0:
            if rxn.id in remove_rxns:
                exclude |= set([rxn.id])
                continue
                
        # Check if reaction in partition dictionary
        try: 
            test = percentiles[rxn.id]
        except KeyError:
            include |= set([rxn.id])
            continue
            
        # Define those reactions not considered in minimization steps
        if percentiles[rxn.id] == 0:
            include |= set([rxn.id])
            continue
        elif percentiles[rxn.id] == 5:
            exclude |= set([rxn.id])
            continue

        # Assess remainder of reactions
        if percentiles[rxn.id] == 1:
            perc_top |= set([rxn.id])
            continue
        elif percentiles[rxn.id] == 2:
            perc_hi |= set([rxn.id])
            continue
        elif percentiles[rxn.id] == 3:
            perc_mid |= set([rxn.id])
            continue
        elif percentiles[rxn.id] == 4:
            perc_lo |= set([rxn.id])
            continue
        else:
            include |= set([rxn.id])  # If not found, automatically include in final model

    return include, exclude, perc_top, perc_hi, perc_mid, perc_lo


# Read in user defined reactions to keep or exclude
def read_user_defined_reactions(reaction_file):

    with open(reaction_file, 'r') as reactions:
        include_rxns = reactions.readline().split(',')
        exclude_rxns = reactions.readline().split(',')

    return include_rxns, exclude_rxns


# Identify active reactions in a model
def find_active_reactions(model):

    solution = model.optimize()
    fluxes = solution.fluxes
    active_rxns = set()
    for reaction, flux in fluxes.items():
        if abs(flux) > 1e-6:
            active_rxns |= set([reaction])

    return active_rxns


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
# Based on:
# Lewis NE, et al.(2010). Omic data from evolved E. coli are consistent with computed optimal growth from
#       genome-scale models. Molecular Systems Biology. 6, 390.
def pFBA_by_percent_of_optimum(model, rxn_ids, optimum_fraction, exclude_from_min=True):
    
    # Minimize flux through all reactions such that the fraction of objective optimum is still achieved
    remove_ids = set()
    keep_ids = set()
    with model as m:
        
        # Fix previous objective as constraint with threshold of predefined fraction of uncontexualized flux
        solver.fix_objective_as_constraint(m, fraction=optimum_fraction)
        
        # Formulate pFBA objective
        if exclude_from_min == True:
            rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id not in rxn_ids)
            excl_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in rxn_ids)
            excl_rxn_vars = chain(*excl_rxn_vars)
        else:
            rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions)
        rxn_vars = chain(*rxn_vars)

        # Set new objective for minimum sum of fluxes
        pfba_obj = m.problem.Objective(Zero, direction='min', sloppy=True)
        m.objective = pfba_obj
        
        # Set linear coefficients based on if they are to be excluded from minimization
        m.objective.set_linear_coefficients({x: 1.0 for x in rxn_vars})
        if exclude_from_min == True:
            m.objective.set_linear_coefficients({y: 0.0 for y in excl_rxn_vars})
        
        # Identify those reactions of interest that carry flux in pFBA solution
        active_rxns = find_active_reactions(m)
        keep_ids |= rxn_ids.intersection(active_rxns)
        remove_ids |= rxn_ids.difference(active_rxns)

    return remove_ids, keep_ids, active_rxns


# Parse reactions to be removed and find adjacent active reactions 
# Based on:
# Jensen PA, et al. (2017). Antibiotics Disrupt Coordination between Transcriptional 
#       and Phenotypic Stress Responses in Pathogenic Bacteria. Cell Reports. 20, 1705-1716.
def reduce_pruning_by_topology(model, keep_ids, remove_ids):

    save_proximal_rxns = set()
    for rxn in model.reactions:

        # Skip reactions not labeled in the keep reactions
        if rxn.id not in keep_ids: 
            continue
        
        # Parse all substrates of the current reactions and those reactions they also participate in
        substrates = rxn.metabolites
        for cpd in substrates:

            # Find adjacent reactions that are currently designated for removal
            adjacent_rxns = cpd.reactions
            save_proximal_rxns |= set([x.id for x in adjacent_rxns]).intersection(remove_ids)
            
    # Amend removal set based on topology parsing
    new_remove_ids = copy.deepcopy(remove_ids)
    for rxn in save_proximal_rxns: 
        new_remove_ids.remove(rxn)
            
    return new_remove_ids


# Prune model and text that contextualized model is still able to grow
def prune_and_test(model, remove_rxn_ids):

    # Prune highlighted reactions from model, removing newly orphaned genes and metabolites
    for rxn in remove_rxn_ids:
        try:
            model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
        except:
            pass
    
    # Remove residual orphaned reactions and metabolites (just in case)
    unused_current_cpd = 1
    unused_current_rxn = 1
    while unused_current_cpd != 0 or unused_current_rxn != 0:
        unused_cpd = prune_unused_metabolites(model)
        unused_rxn = prune_unused_reactions(model)        
        unused_current_cpd = len(unused_cpd)
        unused_current_rxn = len(unused_rxn)
    
    # Check that prune model can still achieve flux through the objective (just in case)
    if model.slim_optimize() < 1e-6: 
        print('WARNING: Pruned model objective can no longer carry flux')
        pass

    return model


# Create context-specific model based on transcript distribution
def contextualize(model, binning_dict, min_fraction=0.01, max_fraction=0.95, report=True, defined_rxns=False):

    # Format user defined reactions for keeping or removal
    if defined_rxns != False:
        user_keep_list, user_remove_list = read_user_defined_reactions(defined_rxns)
    else:
        user_keep_list = []
        user_remove_list = []

    # Partition reactions based on transcription percentile intervals
    ignore, exclude, perc_top, perc_hi, perc_mid, perc_lo = parse_reaction_binning(model, binning_dict, user_keep_list, user_remove_list)
    
    # Generate copy of old model for contextualization
    contextualized_model = copy.deepcopy(model)
    
    # Highest percentile reactions, excluding them in pFBA flux minization at 95% of optimal growth,
    # with the parsing local reactions and saving those with active neighbors
    remove_top, keep_top, active_top = pFBA_by_percent_of_optimum(contextualized_model, perc_top, optimum_fraction=max_fraction)
    final_remove_top = reduce_pruning_by_topology(contextualized_model, keep_top, remove_top)
    contextualized_model = prune_and_test(contextualized_model, final_remove_top)
 
    # Screen second highest percentile reactions, now excluding them in pFBA flux minization with only 1% of optimal growth,
    # with the added step of parsing local reactions and saving those with active neighbors
    remove_hi, keep_hi, active_hi = pFBA_by_percent_of_optimum(contextualized_model, perc_hi, optimum_fraction=min_fraction)
    final_remove_hi = reduce_pruning_by_topology(contextualized_model, keep_hi, remove_hi)
    contextualized_model = prune_and_test(contextualized_model, final_remove_hi)

    # Screen next lowest percentile reactions, including them in pFBA flux minization but with 95% of optimal growth
    remove_mid, keep_mid, active_mid = pFBA_by_percent_of_optimum(contextualized_model, perc_mid, optimum_fraction=max_fraction, exclude_from_min=False)
    contextualized_model = prune_and_test(contextualized_model, remove_mid)

    # Screen lowest percentile reactions, including them in pFBA flux minization with only 1% of optimal growth
    remove_lo, keep_lo, active_lo = pFBA_by_percent_of_optimum(contextualized_model, perc_lo, optimum_fraction=min_fraction, exclude_from_min=False)
    contextualized_model = prune_and_test(contextualized_model, remove_lo)

    # Reactions not considered by minimization
    if len(exclude) > 0: contextualized_model = prune_and_test(contextualized_model, exclude)

    # Cross-reference all active reactions and all removed reactions and add back those where flux is sometimes parsimonious
    total_active = set(list(active_top) + list(active_hi) + list(active_mid) + list(active_lo) + list(ignore))
    total_remove = set(list(final_remove_top) + list(final_remove_hi) + list(remove_mid) + list(remove_lo) + list(exclude))
    overlap_rxns = total_active.intersection(total_remove)
    salvaged_rxns = [model.reactions.get_by_id(rxn) for rxn in overlap_rxns]
    contextualized_model.add_reactions(salvaged_rxns)

    # Report statistics on pruning steps
    contextualized_model.id = str(contextualized_model.id) + '_riptide'
    output_string = """
Contextualized model ID:             {model_id}
Pruned reactions:                    {pruned}
Included reactions:                  {kept}

Highest percentile interval
Total reactions:                     {total_top}
Highlighted for initial pruning:     {excluded_top}
Saved by adjacent reaction:          {saved_top}
Included in pruned model:            {included_top}
-------------------------------------------------
Mid-High percentile interval
Total reactions:                     {total_hi}
Highlighted for initial pruning:     {excluded_hi}
Saved by adjacent reaction:          {saved_hi}
Included in pruned model:            {included_hi}
-------------------------------------------------
Mid-Low percentile interval
Total reactions:                     {total_mid}
Pruned from model:                   {excluded_mid}
Included in pruned model:            {included_mid}
-------------------------------------------------
Lowest percentile interval
Total reactions:                     {total_lo}
Pruned from model:                   {excluded_lo}
Included in pruned model:            {included_lo}
-------------------------------------------------
Reactions recovered by final screen  {salvaged}
-------------------------------------------------
Ommitted from pFBA (partially user-defined)
Pruned from final model:             {excluded}
Included in final model:             {included}
""".format(model_id = contextualized_model.id, 
        pruned = (len(exclude)+len(remove_lo)+len(remove_mid)+len(final_remove_hi)+len(final_remove_top)),
        kept = len(contextualized_model.reactions),
        included = len(ignore),
        excluded = len(exclude),
        total_lo = len(perc_lo), 
        included_lo = len(keep_lo),
        excluded_lo = len(remove_lo),
        total_mid = len(perc_mid),
        included_mid = len(keep_mid), 
        excluded_mid = len(remove_mid), 
        total_hi = len(perc_hi),
        included_hi = (len(keep_hi)+(len(remove_hi)-len(final_remove_hi))),
        excluded_hi = len(remove_hi),
        saved_hi = (len(remove_hi)-len(final_remove_hi)),
        total_top = len(perc_top),
        included_top = (len(keep_top)+(len(remove_top)-len(final_remove_top))),
        excluded_top = len(remove_top),
        saved_top = (len(remove_top)-len(final_remove_top)),
        salvaged = len(salvaged_rxns))

    if report == True: print(output_string)
    
    return contextualized_model


# Compute relative doubling time (in minutes) from objective value
def doubling_time(model):
    
    with model as m:
        ov = m.slim_optimize()
        growth = (1.0 / float(ov)) * 3600.0

    return growth


# Run whole contextualization protocol mutliple times and return fastest growing model
def growth_optimize_by_context(model, binning_dict, min_fraction=0.01, max_fraction=0.95, defined_rxns=False, iters=5):

    curr_iter = 1
    doubling_rate = 100000000.0
    while curr_iter <= iters:
        curr_iter += 1
        current_model = contextualize(model, binning_dict, min_fraction=min_fraction, max_fraction=max_fraction, report=False, defined_rxns=defined_rxns)

        # Test growth rate
        current_rate = doubling_time(current_model)
        if current_rate < doubling_rate:
            fast_model = copy.deepcopy(current_model)
            doubling_rate = current_rate
        else:
            continue

    return fast_model
