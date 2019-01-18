# MatrixMaker
## GUI tool for converting methyole data to a data matrix compatible with Methylica
*MatrixMaker* is a GUI-based tool for converts user's methylome data to a data matrix compatible with [*Methylica*](https://github.com/HiromitsuAraki/Methylica).

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

