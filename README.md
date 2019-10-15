# RIPTiDe

**R**eaction **I**nclusion by **P**arsimony and **T**ranscript **D**istribution

----

Transcriptomic analyses of bacteria have become instrumental to our understanding of their responses to changes in their environment. While traditional analyses have been informative, leveraging these datasets within genome-scale metabolic network reconstructions can provide greatly improved context for shifts in pathway utilization and downstream/upstream ramifications for changes in metabolic regulation. Previous techniques for transcript integration have focused on creating maximum consensus with the input datasets. However, these approaches have collectively performed poorly for metabolic predictions even compared to transcript-agnostic approaches of flux minimization that identifies the most efficient patterns of metabolism given certain growth constraints. Our new method, RIPTiDe, combines these concepts and utilizes overall minimization of flux in conjunction with transcriptomic analysis to identify the most energy efficient pathways to achieve growth that include more highly transcribed enzymes. RIPTiDe requires a low level of manual intervention which leads to reduced bias in predictions. 

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

### Additional parameters for main RIPTiDe functions:

**riptide.read_transcription_file()**
```
file : string
    User-provided file name which contains gene IDs and associated transcription values
header : boolean
    Defines if read abundance file has a header that needs to be ignored
    default is no header
replicates : boolean
    Defines if read abundances contains replicates and medians require calculation
    default is no replicates (False)
sep: string
    Defines what character separates entries on each line
    defaults to tab (.tsv)
binning : boolean
    Perform discrete binning of transcript abundances into quantiles
    OPTIONAL, not advised
    default is False
quant_max : float
    Largest quantile to consider
    default is 0.9
quant_min : float
    Smallest quantile to consider
    default is 0.5
step : float
    Step size for parsing quantiles
    default is 0.125
```

**riptide.contextualize()**
```
model : cobra.Model
    The model to be contextualized (REQUIRED)
transcriptome : dictionary
    Dictionary of transcript abundances, output of read_transcription_file (REQUIRED)
samples : int
    Number of flux samples to collect, default is 500
norm : bool
    Normalize transcript abundances using RPM calculation
    Performed by default
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
    List of reaction ID strings for forced inclusion in final model (metabolic tasks)
exclude : list
    List of reaction ID strings for forced exclusion from final model
gpr : bool
    Determines if GPR rules will be considered during coefficient assignment
    Default is False
threshold : float
    Minimum flux a reaction must acheive in order to avoid pruning during flux sum minimization step
    Default is 1e-6
```

### Example stdout report:
```

Initializing model and integrating transcriptomic data...
Pruning zero flux subnetworks...
Analyzing context-specific flux distributions...

Reactions pruned to 285 from 1129 (74.76% change)
Metabolites pruned to 285 from 1132 (74.82% change)
Flux through the objective DECREASED to ~54.71 from ~65.43 (16.38% change)
Contextualized GENRE is concordant with the transcriptome (p=0.003)

RIPTiDe completed in 1 minute and 15 seconds

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
- **user_defined** - User defined reactions in a 2 element dictionary that either were included or excluded
- **concordance** - Spearman correlation results between linear coefficients and median fluxes from sampling


## Additional Information

Thank you for your interest in RIPTiDe, for additional questions please email mljenior@virginia.edu.

If you encounter any problems, please [file an issue](https://github.com/mjenior/riptide/issues) along with a detailed description.

Distributed under the terms of the [MIT](http://opensource.org/licenses/MIT) license, "riptide" is free and open source software
