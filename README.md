RIPTiDe
=======
Reaction Inclusion by Parsimonious usage and Transcript Distribution


### Please cite when using:
Jenior ML, Moutinha TJ, and Papin JA. (2019). Parsimonious transcript data integration improves context-specific predictions of bacterial metabolism in complex environments. BioRxiv. DOI

### Dependencies
cobrapy (>=version 0.13)

symengine

### Installation:
pip install riptide

### Example usage:
from riptide import *

my_model = cobra.io.read_sbml_model('my_model.sbml')

transcript_abundances = read_transcription_file('reads_to_genes.tsv', replicates=True)

contextualized_model, flux_samples = riptide(my_model, transcript_abundances)
