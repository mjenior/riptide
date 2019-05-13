# RIPTiDe

**R**eaction **I**nclusion by **P**arsimony and **T**ranscript **D**istribution

----

Transcriptomic analyses of bacteria have become instrumental to our understanding of their responses to changes in their environment. While traditional analyses have been informative, leveraging these datasets within genome-scale metabolic network reconstructions can provide greatly improved context for shifts in pathway utilization and downstream/upstream ramifications for changes in metabolic regulation. Previous techniques for transcript integration have focused on creating maximum consensus with the input datasets. However, these approaches have collectively performed poorly for metabolic predictions even compared to transcript-agnostic approaches of flux minimization that identifies the most efficient patterns of metabolism given certain growth constraints. Our new method, RIPTiDe, combines these concepts and utilizes overall minimization of flux in conjunction with transcriptomic analysis to identify the most energy efficient pathways to achieve growth that include more highly transcribed enzymes. RIPTiDe requires a low level of manual intervention which leads to reduced bias in predictions. 

Please cite when using:
```
Jenior ML, Moutinho TJ, and Papin JA. (2019). Parsimonious transcript data integration improves context-specific predictions of bacterial metabolism in complex environments. BioRxiv.
```

## Installation

Installation is simply:
```
$ pip install riptide
```

Or from github:
```
$ pip install git+https://github.com/mjenior/riptide
```

## Usage

```python
from riptide import *

my_model = cobra.io.read_sbml_model('examples/genre.sbml')

transcript_abundances_1 = read_transcription_file('examples/transcriptome1.tsv')
transcript_abundances_2 = read_transcription_file('examples/transcriptome2.tsv')

riptide_object_1 = riptide(my_model, transcript_abundances_1)
riptide_object_2 = riptide(my_model, transcript_abundances_2)
``` 

### Additional parameters for main RIPTiDe functions:

**read_transcription_file()**
```
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
```

**riptide()**
```
model : cobra.Model
    The model to be contextualized
transcription : dictionary
    Dictionary of transcript abundances, output of read_transcription_file()
defined : False or File
    Text file containing reactions IDs for forced inclusion listed on the first line and exclusion 
    listed on the second line (both .csv and .tsv formats supported)
samples : int 
    Number of flux samples to collect, default is 10000, If 0, sampling skipped
percentiles : list of floats
    Percentile cutoffs of transcript abundance for linear coefficient assignments to associated reactions
    Default is [50.0, 62.5, 75.0, 87.5]
coefficients : list of floats
    Linear coefficients to weight reactions based on distribution placement
    Default is [1.0, 0.5, 0.1, 0.01, 0.001]
fraction : float
    Minimum percent of optimal objective value during FBA steps
    Default is 0.8
conservative : bool
    Conservatively remove inactive reactions based on genes
    Default is False
bound : bool
    Bounds each reaction based on transcriptomic constraints
    Default is False
```

### Example stdout report:
```

Initializing model and parsing transcriptome...
Pruning zero flux subnetworks...
Sampling context-specific flux distributions (longest step)...

Reactions pruned to 291 from 1129 (74.22% change)
Metabolites pruned to 289 from 1134 (74.51% change)
Flux through the objective DECREASED to ~76.48 from ~89.77 (14.8% change)
Solution space volume DECREASED to ~1785.89 from ~8460.51 (78.89% change)

RIPTiDe completed in 2 minutes and 41 seconds

```

### Resulting RIPTiDe object properties:

- **model** - contextualized genome-scale metabolic network reconstruction
- **fluxes** - Flux sampling or flux variability analysis pandas object
- **flux_type** - Type of flux analysis performed
- **quantile_range** - percentile intervals by which to parse transcript abundance distribution
- **linear_coefficient_range** - linear coeeficients assigned to corresponding quantile
- **fraction_of_optimum** - minimum percentage of optimal allowable flux through the objective during contextualization


## Additional Information

Thank you for your interest in RIPTiDe, for additional questions please email mljenior@virginia.edu.

If you encounter any problems, please [file an issue](https://github.com/mjenior/riptide/issues) along with a detailed description.

Distributed under the terms of the [MIT](http://opensource.org/licenses/MIT) license, "riptide" is free and open source software
