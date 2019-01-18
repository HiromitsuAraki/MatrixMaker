# MatrixMaker
## GUI tool for converting methyole data to a data matrix compatible with [Methylica](https://github.com/HiromitsuAraki/Methylica)
*MatrixMaker* is a GUI-based tool for converts user's methylome data to a data matrix compatible with *Methylica*. *Methylica* not only provides ICA-based sample clustering but identifies independent components (ICs), or methylomic signatures, specific to various subsets of the samples. Major contributors to a subset-specific IC serve as methylation markers of the subset and imply biological processes underlying the IC. *Methylica* would thus be a powerful tool to analyze samples composed of multiple subtypes, such as those of cancers.
<br>

Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features, respectively. We provide MatrixMaker, which converts user's methylome data to a data matrix compatible with Methylica.


## Install/Launch MatrixMaker
1.  Install [R environment](https://www.r-project.org/)
2.  Install [shiny](https://shiny.rstudio.com).  
`install.packages("shiny")`
3.  Launch *Methylica*  
The following R code will launch *Methylica*.  
`shiny::runGitHub("HiromitsuAraki/MatrixMaker")`
<br>

## Input file format
*Methylica* requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features, respectively. We provide [*MatrixMaker*](https://github.com/HiromitsuAraki/MatrixMaker), which converts user's methylome data to a data matrix compatible with *Methylica*. 
- Methylome data
  - 1st column: Chr
  - 2nd column: Start
  - 3rd column: End
  - 4th column: Gene symbol
  - 5th column ~ : Methylome data of each sample
  <br>
- Sample meta data
  - 1st column: Sample ID
  - 2nd column ~ : Status of the features (e.g. cancer subtype, stage, gender)  
  **NOTE: The status of the features should be discrete, as *Methylica* cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).**  
<br>

## Implementations
### Data uploading
*Methylica* requires methylome data and sample metadata as its inputs. Please refer **Input file format** about the file format of methylome data and sample metadata for *Methylica*. Users need to assign the file location from a browser as below.

<img src="./README_files/Figures/DataUpload.png" width=300x300>
<br>

### Parameter setting
Following data upload, *Methylica* requests its users to select species with its reference genome version, genomic elements to be analyzed (CpG island, gene body, first intron and promoter), and k or the number of ICs (minimum = 2; maximum = the number of samples). *Methylica* provides a default setting of k, defined as the first k components whose cumulative contribution ratio exceeds 80% in principal component analysis. When users select all parameters, users need to press "Run" button to start analysis.

<img src="./README_files/Figures/Parameters.png" width=300x300>
<br>

### Visualizing ICs
When ICA is completed, *Methylica* provides a heatmap of clustering based on the mixing coefficient matrix A, boxplots to compare weightings among the user-defined group categories, and a table of loading values comprising each IC.

#### Heatmap clustering of mixing coefficient matrix
The vertical axis corresponds to mixing coefficients (sample weight in each IC) and the horizontal axis corresponds to samples. Heatmap clustering can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/HeatmapClustering.png" width=500x500>
<br>

#### Boxplots of mixing coefficient matrix
The x axis is sample weights in each IC and the y axis is group categolized based on sample feature. Users can select another sample feature in *Sample property* box as a group defined in sample meta data box. Boxplots can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/Boxplots.png" width=500x500>
<br>

#### Loadings table
Regions with absolute loadings >=2 are appeared in this table. Users can retrieve user's interesting genes in *Search* text box.

<img src="./README_files/Figures/ICtable.png" width=500x500>
<br>


### highlighting highly contributed regions
Once an IC of interest is identified, *Methylica* can display its high loading genomic regions as a methylation heatmap and a clickable list to select a region whose methylation levels in group categories are indicated as boxplots. The list can be exported for various analyses (e.g. GO enrichment and sequence motif search) to facilitate biological interpretation of the IC.

#### Heatmap clustering of highly contributed regions
Users can assign IC and threshold of loadings to be visualized in "IC" and "Loadings" text box respectively. The vertical axis corresponds to genomic regions with high loadings and the horizontal axis corresponds to samples. The genomic regions are sorted by loadings. Heatmap clustering can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/HeatmapClustering_highLF.png" width=500x500>
<br>

#### Loadings table of highly contributed regions
Boxplot of methylome data in each region is appeared by clicking gene symbol (see **Boxplot of highly contributed region**). Whole list of loadings table can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/ICtable_highLF.png" width=500x500>
<br>

#### Boxplot of highly contributed region
The x axis is methylation level and the y axis is each group categolized based on sample feature. Users can select another sample feature as a group defined in sample meta data in *Sample property* box. 

<img src="./README_files/Figures/Boxplots_highLF.png" width=500x500>
<br>
