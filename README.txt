=======
RIPTiDe
=======
Reaction Inclusion by Parsimonious usage and Transcript Distribution


# Please cite when using:
Jenior ML and Papin JA. (2018) RIPTiDe: Metabolic model contextualization through exploration 
of mutliple transcription and flux minimization levels. BioRxiv. DOI

# Dependencies
cobrapy (>=version 0.14)

symengine

# Example usage:
my_model = cobra.io.read_sbml_model('my_model.sbml')

transcript_abundances = read_transcription_file('reads_to_genes.tsv', replicates=True)

contextualized_model, flux_samples = riptide(my_model, transcript_abundances)
