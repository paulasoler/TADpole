# TADpole

TADpole is a computational tool designed to identify and analyze the entire hierarchy of topologically associated domains (TADs) in intra-chromosomal interaction matrices.

## 1) Installation

<!--
### 1.1) Using the _devtools_ package

This is the recommended installation procedure.

- First, install the devtools package from CRAN

```
install.packages("devtools")
```

- Then, Install the TADpole package from GitHub

```
devtools::install_github("paulasoler/TADpole")
```

### 1.2) Manual installation from source
-->

- First, install the required dependencies in R

```
install.packages(c('bigmemory', 'data.table', 'reshape2', 'pryr',
'ggpubr','ggplot2','ggdendro','plyr','zoo','cowplot','gridExtra', 'viridis', 'purrr',
'dendextend', 'doParallel', 'foreach', 'fpc', 'Matrix', 'rioja'))
```

- Then, get the latest version of the source code from Github

by using _wget_:

```
wget https://github.com/paulasoler/TADpole/archive/master.zip
unzip master.zip
mv master TADpole
```

or by cloning the repository:

```
git clone https://github.com/paulasoler/TADpole.git
```

- Finally, install TADpole package


```
R CMD INSTALL TADpole
```

## 2) Getting started

In this repository, we provide a test case from a publicly available Hi-C data set (SRA: [SRR1658572](https://www.ebi.ac.uk/ena/data/view/SRR1658572)) (1).

In the `inst/extdata/` directory, we provided a 3Mb-region (chr18:9,200,000-12,120,000) of a human Hi-C dataset at 20kb resolution. 

```
- inst/extdata/raw_chr18:460-606_20kb.mat
```

To obtain this interaction matrix, we processed the Hi-C data using the [TADbit](https://github.com/3DGenomes/TADbit) (2) Python library, that deals with all the necessary steps to analyze and normalize Hi-C data.

### 2.1) Input data
To run the main function `TADpole`, you need to provide an intrachromosomal interaction matrix, representing an entire chromosome, or a continuous chromosome region. The input is a generic tab-separated file containing the interaction matrix (M) with N rows and N columns, where N is the number of bins in which the chromosome region is divided. Each position of the matrix (Mij) contains the number of interaction values (raw or normalized) between the corresponding pair of genomic bins i and j. We recommend [ONED (https://github.com/qenvio/dryhic) (3) normalization, as it effectively corrects for known experimental biases.


### 2.2) Running the algorithm

Schematic overview of the TADpole algorithm (for further details, refer to Soler-Vila et.al [4] (https://github.com/paulasoler/TADpole) 

![Zoom](https://github.com/paulasoler/TADpole/blob/master/misc/Figure1.png)

The basic usage is the following:

```
library(TADpole)
chr18_460-606_20kb <- system.file("extdata", "raw_chr18:460-606_20kb.mat", package = "TADpole")

tadpole <- TADpole(mat_file = chr18_460-606_20kb, 
chr = "chr18", start = 9200000, end = 12120000, resol = 20000,bad_frac = 0.01, centromere_search = FALSE)
```

#### 2.2.1) Parameters
- **mat_file**: `path` to the input file. Must be in a tab-delimited matrix format.
- **chr**: `string` with the chromosome name.
- **start**: `numeric` initial position of the chromosomal region or the chromosome.
- **end**: `numeric` final position of the chromosomal region or the chromosome.
- **resol**: `numeric` binning-size of the Hi-C experiment (in base-pairs)
- **max_pcs**: `numeric` the maximum number of principal components to retain for the analysis. Default value of 200 is recommended.
- **min_clusters**: `numeric` minimum number of chromatin partitions.
- **bad_frac**: `numeric` fraction of the matrix to flag as bad columns.
- **hist_bad_columns**: `logical` plot the distribution of column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging purposes.
- **centromere_search**: `logical` split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) matrices.


## 3) Output
The function `TADpole` returns a `tadpole` object containing the following descriptors:

- ***n_pcs***: optimal number of principal components (NPCs*).
- ***optimal_n_clusters***: optimal number of chromatin partitions (that is the index of the optimal level (ND*) plus 1).
- ***dendro***: hierarchical tree-like structure cut at the maximum significant number of levels identified by the broken-stick model (max(ND)).
- ***clusters***: a list containing the chromatin partitions per each hierarchical level _(x)_ defined by the broken stick model.
  + ***clusters$`x`***: start and end coordinades of all chromatin partitions.
- ***score***: CH index associated to each dendrogram.
- ***merging_arms***: if `centromere_search` is `TRUE`, contains the start and end coordinates of the TADs of the full chromosome.

```
head(tadpole)

$n_pcs
[1] 20

$optimal_n_clusters
[1] 12

$dendro

Call:
rioja::chclust(d = dist(pcs))

Cluster method   : coniss
Distance         : euclidean
Number of objects: 198


$clusters
$clusters$`2`
  start end
1     1 110
2   111 200

attr(,"scores")
     1        2        3        4        5        6        7        8        9
1   NA 47,90916 42,22857 39,40353 43,61547 41,24569  0,00000  0,00000  0,00000
2   NA 44,47879 43,28183 45,06219 44,02830 45,38542 49,09032  0,00000  0,00000

```
### 3.1) Plotting the results

#### 3.1.1) Raw Hi-C plot and histogram of the interaction values.
Automatically, TADpole generates a heatmap of the intra-chromosomal interaction matrix under study, together with a histogram to know the distribution of the Hi-C interaction values. In the latter, we can see a dashed line that delimits the number of excluded columns (and the corresponding rows) of the analysis by presenting a low number of interactions (called as bad columns).Specifically, the columns (rows) that contain an empty cell at the main diagonal, and those whose cumulative interactions are below the first (by default) percentile, are excluded from the analysis.

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure2.png" width="70%">
</p>

#### 3.1.2) Hierarchical plot
**Left**, complete dendrogram of the Hi-C matrix cut at a maximum significant number of levels (max(ND)) reported by 
the broken-stick model (containing from 2 to 16 partitions) and, from them, the highest scoring level according to the CH index is selected. **Right**, Hi-C contact map showing the complete hierarchy of the significant levels selected by the BS model (black lines) along with the optimal one in 12 specific partitions, as identified by the highest CH index (blue line).

```
hierarchical_plot(mat_file = chr18_460-606_20kb, chr = "chr18", start = 9200000, end = 12120000, resol = 20000,
tadpole = tadpole, centromere_search=FALSE)
```
##### 3.1.2.1) Parameters
- **mat_file**: `path` to the input file. Must be in a tab-delimited matrix format.
- **tadpole**: `tadpole` object
- **chr**: `string` with the chromosome name.
- **start**: `numeric` initial position of the chromosomal region or the chromosome.
- **end**: `numeric` final position of the chromosomal region or the chromosome.
- **resol**: `numeric` binning-size of the Hi-C experiment (in base-pairs).
- **centromere_search**: `logical` split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) matrices.


<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure3.png" width="70%">
</p>

#### 3.1.3) Matrix of Calinski-Harabasz  indexes 

```
CH_map(tadpole)
```

##### 3.1.3.1) Parameters
- **tadpole**: `tadpole` object.

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure4.png" width="50%" align="center">
</p>

# DiffT Score
To compare pairs of topological partitions, P and Q, identified by TADpole at a fixed level of the hierarchy, we defined a difference topology score (DiffT). Specifically, the partitioned matrices were transformed into binary forms p for P, and analogously q for Q, in which each entry pij (qij) is equal to 1 if the bins i and j are in the same TAD and 0 otherwise. Then, DiffT is computed as the normalized (from 0 to 1) difference between the binarized matrices as a function of the bin index b as:

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/DiffT_formula.png" width="30%" align="center">
</p>

where N is the total number of bins.
<br>

### 1) Input data
The DiffT score analysis was used to compare the chromatin partitions obtanied from a fixed hierarchical level determined in two different experiments, control and case.

In the `data/` directory, there are 2 files, control and case in a BED-like `data.frame`.

```
- data/control.bed
- data/case.bed
```

### 2) Computing the DiffT score

```
control <- read.table(system.file("extdata", "control.bed", package = "TADpole"))
case <- read.table(system.file("extdata", "case.bed", package = "TADpole"))

difft_control_case = diffT(as.data.frame(control[,c(1,2,3)]),
                                as.data.frame(case[,c(1,2,3)]))
```

#### 2.1) Parameters
- **bed_x**, **bed_y**: two `data.frame`s with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.

### 3) Output
The function `diffT` returns a `numeric` vector representing the cumulative DiffT score score profiles as a function of the matrix bins.
The highest local differences between the two matrices can be identified by the sharpest changes in the slope of the function.

```
difft_melt = melt(difft_control_case)
difft_melt$bin = seq(nrow(difft_melt))
difft_melt$level = level

difft_melt$level = as.factor(difft_melt$level)
ggline(difft_melt, x = "bin", y = "value",color = "level", plot_type = "l") 
``````

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/DiffT_score.png" width="60%" align="center">
</p>

## Authors

- **Paula Soler Vila** - (https://github.com/paulasoler/)
- **Pol Cuscó Pons** - (https://github.com/nanakiksc/)
- **Marco Di Stefano** - (https://github.com/MarcoDiS)

## References

1. RAO, Suhas SP, et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell, 2014, 159.7: 1665-1680.
2. SERRA, François, et al. Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors. PLoS computational biology, 2017, 13.7: e1005665.
3. VIDAL, Enrique, et al. OneD: increasing reproducibility of Hi-C samples with abnormal karyotypes. Nucleic acids research, 2018, 46.8: e49-e49.
