# RIPTiDe

**R**eaction **I**nclusion by **P**arsimony and **T**ranscript **D**istribution

v3.4.2
----

Transcriptomic analyses of bacteria have become instrumental to our understanding of their responses to changes in their environment. While traditional analyses have been informative, leveraging these datasets within genome-scale metabolic network reconstructions (GENREs) can provide greatly improved context for shifts in pathway utilization and downstream/upstream ramifications for changes in metabolic regulation. Many previous techniques for GENRE transcript integration have focused on creating maximum consensus with input datasets, but these approaches have been shown to generate less accurate metabolic predictions than a transcript-agnostic method of flux minimization (pFBA), which identifies the most efficient/economic patterns of metabolism given certain growth constraints. Despite this success, growth conditions are not always easily quantifiable and highlights the need for novel platforms that build from these findings. This method, known as RIPTiDe, combines these concepts and utilizes overall minimization of flux weighted by transcriptomic analysis to identify the most energy efficient pathways to achieve growth that include more highly transcribed enzymes, without previous insight into extracellular conditions. This platform could be important for revealing context-specific bacterial phenotypes in line with governing principles of adaptive evolution, that drive disease manifestation or interactions between microbes.

Please cite when using:
```
Jenior ML, Moutinho Jr TJ, Dougherty BV, & Papin JA. (2020). Transcriptome-guided parsimonious flux analysis improves predictions with metabolic networks in complex environments. PLOS Comp Biol. https://doi.org/10.1371/journal.pcbi.1007099.
```

## Dependencies
```
>=python-3.6.4
>=cobra-0.15.3
>=pandas-0.24.1
>=symengine-0.4.0
>=scipy-1.3.0
```

## Installation

Installation is:
```
$ pip install riptide
```

### Arguments for core RIPTiDe functions:

**riptide.read_transcription_file() - Generates dictionary of transcriptomic abundances from a file**
```
REQUIRED
file : string
    User-provided file name which contains gene IDs and associated transcription values

OPTIONAL
header : boolean
    Defines if read abundance file has a header that needs to be ignored
    Default is no header
sep: string
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
    Smallest quantile to consider
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
```

**riptide.contextualize() - Create context-specific model based on transcript distribution with maximum fit of flux distribution to input transctiptome**
```
REQUIRED
model : cobra.Model
    The model to be contextualized
transcriptome : dictionary
    Dictionary of transcript abundances, output of read_transcription_file()

OPTIONAL
cpus  : int
    CPUs number for parallelization
    Default is all available 
frac_min : float
    Lower bound for range of minimal fractions to test
    Default is 0.25
frac_max : float
    Upper bound for range of minimal fractions to test
    Default is 0.85
samples : int 
    Number of flux samples to collect
    Default is 500
silent  : bool
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
    Default is 1e-6
defined : False or list
    User defined range of linear coeffients, needs to be defined in a list like [1, 0.5, 0.1, 0.01, 0.001]
    Works best paired with binned abundance catagories from riptide.read_transcription_file()
    Default is False
open_exchanges : bool
    Sets all exchange reactions bounds to (-1000., 1000.)
    Default is False
skip_fva : bool
    Skips final flux variability analysis
    Default is False
```

**riptide.single_contextualize() - Create context-specific model based on transcript distribution with user-defined objective flux minimum**
```
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

* Other arguments from iterative implementation are carried over (excluding frac_min and frac_max)
```

**riptide.save_output() - Writes RIPTiDe results to files in a new directory**
```
REQUIRED
riptide_obj : RIPTiDe object
    Class object creared by riptide.contextualize()

OPTIONAL
path : str
    New directory to write output files
file_type : str
    Type of output file for RIPTiDe model
    Accepts either sbml or json
    Default is JSON
```

## Usage

**Comments before starting:** 
- Make sure that genes in the transcriptome file matches those that are in your model.
- Check the example files for proper data formatting
- Binning genes into discrete thresholds for coefficient assignment is available in riptide.read_transcription_file() (not recommended)
- Opening the majority of exchange reactions (bounds = +/- 1000) may yeild better prediction when extracellular conditions are unknown
- The resulting RIPTiDe object has multiple properties including the context-specific model and flux analyses, accessing each is described below

```python
import riptide

my_model = cobra.io.read_sbml_model('tests/genre.sbml')

transcript_abundances_1 = riptide.read_transcription_file('tests/transcriptome1.tsv')
transcript_abundances_2 = riptide.read_transcription_file('tests/transcriptome2.tsv') # has replicates

riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)
riptide_object_1_b = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1, tasks=['rxn1'], exclude=['rxn2','rxn3'])

riptide.save_output(riptide_obj=riptide_object_1_a, path='~/Desktop/example_riptide_output')
``` 

### Example riptide.contextualize() stdout report:
```

Initializing model and integrating transcriptomic data...
Pruning zero flux subnetworks...
Analyzing context-specific flux distributions...

Running max fit RIPTiDe for objective fraction range: 0.65 to 0.85
Progress: 100%

Testing local fractions to 0.3...
Progress: 100%

Context-specific metabolism fit with 0.35 of optimal objective flux

Reactions pruned to 285 from 1129 (74.76% change)
Metabolites pruned to 285 from 1132 (74.82% change)
Flux through the objective DECREASED to ~54.71 from ~65.43 (16.38% change)
Context-specific metabolism correlates with transcriptome (r=0.149, p=0.011 *)

Max fit RIPTiDe completed in, 4 minutes and 33 seconds  

```

In the final step, RIPTiDe assesses the fit of transcriptomic data for the calculated network activity through correlation of transcript abundance and median flux value for each corresponding reaction. The Spearman correlation coefficient and associated p-value are the reported following the fraction of network topology that is pruned during the flux minimization step.

Max fit RIPTiDe tests all minimum objective flux fractions over the provided range and returns only the model with the best Spearman correlation between context-specific flux for reactions and the associated transcriptomic values. Note, terminating search if a subsequent iteration has a lower correlation coefficient than the last could result from a local maxima and must be considered if an exhaustive analysis is preferred.

### Resulting RIPTiDe object (class) properties:
The resulting object is a container for the following data structures.

- **model** - Contextualized genome-scale metabolic network reconstruction
- **transcriptome** - Transcriptomic replicate abundances provided by user
- **percent_of_mapping** - Percent of genes in mapping file found in input GENRE
- **minimization_coefficients** - Linear coefficients used during flux sum minimization (based on transcriptome replicates)
- **maximization_coefficients** - Linear coefficients for each reaction based used during flux sampling
- **pruned** - Dictionary containing the IDs of genes, reactions, and metabolites pruned by RIPTiDe
- **flux_samples** - Flux samples from constrained model
- **flux_variability** - Flux variability analysis from constrained model
- **fraction_of_optimum** - Minimum specified percentage of optimal objective flux during contextualization
- **metabolic_tasks** - User defined reactions whose activity is saved from pruning
- **concordance** - Spearman correlation results between linear coefficients and median fluxes from sampling
- **gpr_integration** - Whether GPR rules were considered during assignment of linear coefficients
- **defined_coefficients** - Range of linear coefficients RIPTiDe is allowed to utilize provided as a list
- **included_important** - Reactions or Genes included in the final model which the user defined as important
- **additional_parameters** - Dictionary of additional parameters RIPTiDe uses
- **fraction_bounds** - Minimum and maximum values for the range of objective flux minimum fractions tested
- **maxfit_iters** - Objective flux and fit to transcriptome for each minimum flux fraction tested

**Examples of accessing components of RIPTiDe output:**
```python
context_specific_GENRE = riptide_object.model
context_specific_FVA = riptide_object.flux_variability
context_specific_flux_samples = riptide_object.flux_samples
```

## Additional Information

Thank you for your interest in RIPTiDe!

If you encounter any problems, please [file an issue](https://github.com/mjenior/riptide/issues) along with a detailed description.

Distributed under the terms of the [MIT](http://opensource.org/licenses/MIT) license, "riptide" is free and open source software
