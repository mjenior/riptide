# RIPTiDe

**R**eaction **I**nclusion by **P**arsimony and **T**ranscript **D**istribution

----

Transcriptomic analyses of bacteria have become instrumental to our understanding of their responses to changes in their environment. While traditional analyses have been informative, leveraging these datasets within genome-scale metabolic network reconstructions (GENREs) can provide greatly improved context for shifts in pathway utilization and downstream/upstream ramifications for changes in metabolic regulation. Many previous techniques for GENRE transcript integration have focused on creating maximum consensus with input datasets, but these approaches have been shown to generate less accurate metabolic predictions than a transcript-agnostic method of flux minimization (pFBA), which identifies the most efficient/economic patterns of metabolism given certain growth constraints. Despite this success, growth conditions are not always easily quantifiable and highlights the need for novel platforms that build from these findings. This method, known as RIPTiDe, combines these concepts and utilizes overall minimization of flux weighted by transcriptomic analysis to identify the most energy efficient pathways to achieve growth that include more highly transcribed enzymes, without previous insight into extracellular conditions. This platform could be important for revealing context-specific bacterial phenotypes in line with governing principles of adaptive evolution, that drive disease manifestation or interactions between microbes.

Please cite when using:
```
Jenior ML, Moutinho TJ, and Papin JA. (2019). Transcriptome-guided parsimonious flux analysis improves predictions with metabolic networks in complex environments. bioRxiv 637124; doi: https://doi.org/10.1101/637124
```

Utilizes python implementation of the gapsplit flux sampler. Please also cite:
```
Keaty TC and Jensen PA (2019). gapsplit: Efficient random sampling for non-convex constraint-based models. bioRxiv 652917; doi: https://doi.org/10.1101/652917
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

**riptide.read_transcription_file()**
```
file : string
    User-provided file name which contains gene IDs and associated transcription values
header : boolean
    Defines if read abundance file has a header that needs to be ignored
    Default is no header
replicates : boolean
    Defines if read abundances contains replicates and medians require calculation
    Default is no replicates (False)
sep: string
    Defines what character separates entries on each line
    Defaults to tab (.tsv)
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

**riptide.contextualize()**
```
model : cobra.Model
    The model to be contextualized (REQUIRED)
transcriptome : dictionary
    Dictionary of transcript abundances, output of read_transcription_file (REQUIRED)
samples : int
    Number of flux samples to collect
    Default is 500
silent  : bool
    Silences std out 
    Default is False
exch_weight : bool
    Weight exchange reactions the same ase adjacent transporters
fraction : float
    Minimum percent of optimal objective value during FBA steps
    Default is 0.8
minimum : float
    Minimum linear coefficient allowed during weight calculation for pFBA
    Default is None
conservative : bool
    Conservatively remove inactive reactions based on genes
    Default is False
bound : bool
    Bounds each reaction based on transcriptomic constraints
    Default is False
objective : bool
    Sets previous objective function as a constraint with minimum flux equal to user input fraction
    Default is True
set_bounds : bool
    Uses flax variability analysis results from constrained model to set new bounds for all reactions
    Default is True
tasks : list
    List of reaction ID strings for metabolic tasks to be included in final model 
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
    OPTIONAL, default is False
```

## Usage

```python
from riptide import *

my_model = cobra.io.read_sbml_model('examples/genre.sbml')

transcript_abundances_1 = riptide.read_transcription_file('examples/transcriptome1.tsv')
transcript_abundances_2 = riptide.read_transcription_file('examples/transcriptome2.tsv', replicates=True)

riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)
riptide_object_1_b = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1, tasks=['rxn1'], exclude=['rxn2','rxn3'])
riptide_object_2 = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_2)
``` 

### Example stdout report:
```

Initializing model and integrating transcriptomic data...
Pruning zero flux subnetworks...
Analyzing context-specific flux distributions...

Reactions pruned to 285 from 1129 (74.76% change)
Metabolites pruned to 285 from 1132 (74.82% change)
Flux through the objective DECREASED to ~54.71 from ~65.43 (16.38% change)
Contextualized GENRE is concordant with the transcriptome (p=0.003 *)

RIPTiDe completed in 29 seconds

```

### Resulting RIPTiDe object (class) properties:

- **model** - Contextualized genome-scale metabolic network reconstruction
- **transcriptome** - Transcriptomic abundances provided by user
- **percent_of_mapping** - Percent of genes in mapping file found in input GENRE
- **minimization_coefficients** - Linear coefficients used during flux sum minimization
- **maximization_coefficients** - Linear coefficients for each reaction based used during flux sampling
- **flux_samples** - Flux samples from constrained model
- **flux_variability** - Flux variability analysis from constrained model
- **fraction_of_optimum** - Minimum specified percentage of optimal objective flux during contextualization
- **metabolic_tasks** - User defined reactions whose activity is saved from pruning
- **concordance** - Spearman correlation results between linear coefficients and median fluxes from sampling
- **gpr_integration** - Whether GPR rules were considered during assignment of linear coefficients
- **defined_coefficients** - Range of linear coefficients RIPTiDe is allowed to utilize provided as a list


## Additional Information

Thank you for your interest in RIPTiDe, for additional questions please email mljenior@virginia.edu.

If you encounter any problems, please [file an issue](https://github.com/mjenior/riptide/issues) along with a detailed description.

Distributed under the terms of the [MIT](http://opensource.org/licenses/MIT) license, "riptide" is free and open source software
