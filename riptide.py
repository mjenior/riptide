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
import time
import cobra
from cobra.util import solver
from optlang.symbolics import Zero
from itertools import chain
from cobra.manipulation.delete import *
from cobra.flux_analysis.sampling import OptGPSampler

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


# Read in user defined reactions to keep or exclude
def read_user_defined_reactions(reaction_file):

    with open(reaction_file, 'r') as reactions:
        include_rxns = reactions.readline().split(',')
        exclude_rxns = reactions.readline().split(',')

    return include_rxns, exclude_rxns


# Converts a dictionary of transcript distribution percentiles
# Based on:
# Schultz, A, & Qutub, AA (2016). Reconstruction of Tissue-Specific Metabolic Networks Using CORDA. 
#      PLoS Computational Biology. https://doi.org/10.1371/journal.pcbi.1004808
def create_reaction_partitions(model, transcription, save_bins_as=False, high_cutoff=90, mid_cutoff=75, low_cutoff=50, defined_rxns=False):
    
    if defined_rxns != False:
        user_keep_list, user_remove_list = read_user_defined_reactions(defined_rxns)
    else:
        user_keep_list = []
        user_remove_list = []
        
    # Define transcript abundance cutoffs
    distribution = transcription.values()
    perc_hi = numpy.percentile(distribution, high_cutoff)
    perc_mid = numpy.percentile(distribution, mid_cutoff)
    perc_lo = numpy.percentile(distribution, low_cutoff)
    percentile_dict = {'cutoffs': [high_cutoff, mid_cutoff, low_cutoff]}
    
    # Write to file if defined by user
    if save_bins_as != False:
        perc_file = open(save_bins_as, 'w')

    # Keep a count of genes included in each partition
    count_4 = 0
    count_3 = 0
    count_2 = 0
    count_1 = 0
    
    for gene in model.genes:
        
        # If current gene not in transciptomic profile, default to lowest included group
        failed = 0
        try:
            test = transcription[gene.id]
        except KeyError:
            curr_percentile = 4
            failed = 1
            count_4 += 1
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
        user_defined = user_keep_list + user_remove_list
        for rxn in list(gene.reactions): 
            if rxn in user_defined:
                continue
            else:
                percentile_dict[rxn.id] = curr_percentile
            
    for rxn in user_keep_list:
        percentile_dict[rxn] = 0
    for rxn in user_remove_list:
        percentile_dict[rxn] = 5

    if save_bins_as != False:
        perc_file.close()
    
    # Report sizes of each bin
    print('Highest: ' + str(count_1))
    print('High-Mid: ' + str(count_2))
    print('Mid-Low: ' + str(count_3))
    print('Lowest: ' + str(count_4))
    
    return percentile_dict
       
        
# Bin reactions based on their percentile transcription
def parse_reaction_partitions(model, percentiles):

    included = set()
    excluded = set()
    perc_top = set()
    perc_hi = set()
    perc_mid = set()
    perc_lo = set()
    
    # Assign bins associated with each percentile of the read abundance distribution
    for rxn in model.reactions:

        # Check if reaction in partition dictionary
        try: 
            test = percentiles[rxn.id]
        except KeyError:
            included |= set([rxn.id])
            continue
            
        # Define those reactions not considered in minimization steps
        if percentiles[rxn.id] == 0:
            included |= set([rxn.id])
            continue
        elif percentiles[rxn.id] == 5:
            excluded |= set([rxn.id])
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

    return included, excluded, perc_top, perc_hi, perc_mid, perc_lo


# Read in user-defined/editted gene priority binning file
def read_binning_file(partition_file):

    bin_dict = {}
    with open(partition_file, 'r') as percentiles:
        for line in percentiles:
            line = line.split()
            bin_dict[line[0]] = float(line[1])

    return bin_dict


# Identify active reactions in a model
def find_active_reactions(model, sample):

    active_rxns = set()
    
    # Use flux sampling to identify all active reactions
    if sample == True:
        optgp_object = OptGPSampler(iCdJ794_clinda, processes=4)
        flux_samples = optgp_object.sample(1000)

        for distribution in flux_samples.iterrows():
            for rxn in list(flux_samples.columns):
                if abs(distribution[1][rxn]) > 1e-6:
                    active_rxns |= set([rxn])
    
    # Use a single FBA solution to find active reactions
    else:
        solution = model.optimize()
        fluxes = solution.fluxes
        for reaction, flux in fluxes.items():
            if abs(flux) > 1e-6:
                active_rxns |= set([reaction])
            
    return active_rxns


# Determine those reactions that carry flux in a pFBA objective set to a threshold of maximum
# Based on:
# Lewis NE, et al.(2010). Omic data from evolved E. coli are consistent with computed optimal growth from
#       genome-scale models. Molecular Systems Biology. 6, 390.
def minimize_flux_by_percentile(model, partition_dict, adjacent, sample):
    
    # Minimize flux through all reactions such that the fraction of objective optimum is still achieved
    remove_ids = set()
    keep_ids = set()
    with model as m:
        
        # Partition reactions based on transcription percentile intervals
        included, excluded, top_ids, hi_ids, mid_ids, lo_ids = parse_reaction_partitions(model, partition_dict)
        
        # Identify adjacent reactions to hightranscription reactions
        if adjacent == True:
            adj_top_ids = reduce_pruning_by_topology(model, top_ids)
        
        # Fix previous objective as constraint with threshold of predefined fraction of uncontexualized flux
        solver.fix_objective_as_constraint(m, fraction=0.97)
        
        # Formulate pFBA objective
        rxn_ids = set([str(rxn.id) for rxn in m.reactions])
        other_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in rxn_ids)
        other_rxn_vars = chain(*other_rxn_vars)
        top_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in top_ids)
        top_rxn_vars = chain(*top_rxn_vars)
        hi_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in hi_ids)
        hi_rxn_vars = chain(*hi_rxn_vars)
        mid_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in mid_ids)
        mid_rxn_vars = chain(*mid_rxn_vars)
        lo_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in lo_ids)
        lo_rxn_vars = chain(*lo_rxn_vars)
        
        if adjacent == True:
            adj_rxn_vars = ((rxn.forward_variable, rxn.reverse_variable) for rxn in m.reactions if rxn.id in adj_top_ids)
            adj_rxn_vars = chain(*adj_rxn_vars)
        
        # Set new objective for minimum sum of fluxes
        pfba_obj = m.problem.Objective(Zero, direction='min', sloppy=True)
        m.objective = pfba_obj
        
        # Set linear coefficients based on if they are to be excluded from minimization
        m.objective.set_linear_coefficients({x: 1.0 for x in other_rxn_vars})
        m.objective.set_linear_coefficients({x: 0.0 for x in top_rxn_vars})
        m.objective.set_linear_coefficients({x: 0.25 for x in hi_rxn_vars})
        m.objective.set_linear_coefficients({x: 0.75 for x in mid_rxn_vars})
        m.objective.set_linear_coefficients({x: 1.0 for x in lo_rxn_vars})
        
        if adjacent == True:
            m.objective.set_linear_coefficients({x: 0.0 for x in adj_rxn_vars})

        # Get active sections of network  
        active_rxns = find_active_reactions(m, sample)
        remove_ids = rxn_ids.difference(active_rxns)

        # Screen reactions based on user definitions
        remove_ids = remove_ids.difference(included)
        remove_ids |= excluded
        
    return remove_ids


# Prune model and text that contextualized model is still able to grow
def prune_and_test(model, remove_rxn_ids):
    
    new_model = copy.deepcopy(model)
    
    # Prune highlighted reactions from model, removing newly orphaned genes and metabolites
    for rxn in remove_rxn_ids:
        try:
            new_model.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
        except:
            pass
    
    # Remove residual orphaned reactions and metabolites (just in case)
    unused_current_cpd = 1
    unused_current_rxn = 1
    while unused_current_cpd != 0 or unused_current_rxn != 0:
        unused_cpd = prune_unused_metabolites(new_model)
        unused_rxn = prune_unused_reactions(new_model)        
        unused_current_cpd = len(unused_cpd)
        unused_current_rxn = len(unused_rxn)
    
    # Check that prune model can still achieve flux through the objective (just in case)
    if new_model.slim_optimize() < 1e-6: 
        print('WARNING: Pruned model objective can no longer carry flux')

    return new_model


# Parse reactions to be removed and find adjacent active reactions 
# Based on:
# Jensen PA, et al. (2017). Antibiotics Disrupt Coordination between Transcriptional 
#       and Phenotypic Stress Responses in Pathogenic Bacteria. Cell Reports. 20, 1705-1716.
def reduce_pruning_by_topology(model, rxn_ids):

    adjacent_rxns = set()
    for rxn in rxn_ids:
        rxn = model.reactions.get_by_id(rxn)
        
        # Parse all substrates of the current reactions
        substrates = rxn.metabolites
        for cpd in substrates:

            # Find adjacent reactions
            substrate_rxns = cpd.reactions
            adjacent_rxns |= set([x.id for x in substrate_rxns]).difference(rxn_ids)
            
    return adjacent_rxns


# Create context-specific model based on transcript distribution
def contextualize(model, partitions, adjacent=True, sample=False):
    
    # Identify reactions for pruning
    remove_rxns = minimize_flux_by_percentile(model, partitions, adjacent, sample)
    
    # Generate contextualized model
    contextualized_model = prune_and_test(model, remove_rxns)
    contextualized_model.id = str(contextualized_model.id) + '_riptide'
    
    # Report on pruning
    print('Reactions pruned from originial model: ' + str(len(remove_rxns)))
    print('Reactions included in contextualized model: ' + str(len(contextualized_model.reactions)))
    
    return contextualized_model
