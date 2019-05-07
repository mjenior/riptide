RIPTiDe
=======
Reaction Inclusion by Parsimonious usage and Transcript Distribution


### Please cite when using:
Jenior ML, Moutinho TJ, and Papin JA. (2019). Parsimonious transcript data integration improves context-specific predictions of bacterial metabolism in complex environments. BioRxiv. DOI 

### Dependencies
cobrapy (>=version 0.13)

symengine

### Installation:
pip install riptide

### Example usage:
from riptide import *

my_model = cobra.io.read_sbml_model('examples/genre.sbml')

transcript_abundances_1 = read_transcription_file('examples/transcriptome1.tsv')

transcript_abundances_2 = read_transcription_file('examples/transcriptome2.tsv')

contextualized_model_1, flux_samples_1 = riptide(my_model, transcript_abundances_1)

contextualized_model_2, flux_samples_2 = riptide(my_model, transcript_abundances_2)


    read_transcription_file() parameters
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


    riptide() parameters
    ----------
    model : cobra.Model
        The model to be contextualized
    transcription : dictionary
        Dictionary of transcript abundances, output of read_transcription_file()
    defined : False or File
        Text file containing reactions IDs for forced inclusion listed on the first line and exclusion 
        listed on the second line (both .csv and .tsv formats supported)
    sampling : int or False
        Number of flux samples to collect, default is 10000, If False, sampling skipped
    percentiles : list of floats
        Percentile cutoffs of transcript abundance for linear coefficient assignments to associated reactions
        Defaults are [50.0, 62.5, 75.0, 87.5]
    coefficients : list of floats
        Linear coefficients to weight reactions based on distribution placement
        Defaults are [1.0, 0.5, 0.1, 0.01, 0.001]
    fraction : float
        Minimum percent of optimal objective value during FBA steps
        Default is 0.8
    conservative : str
    	Conservatively remove inactive reactions based on GPR rules
    	Either 'y' or 'n', default in 'n' (no)

