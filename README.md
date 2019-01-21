# MatrixMaker
## GUI tool for converting methyole data to a data matrix compatible with Methylica
*MatrixMaker* is a GUI-based tool for converts methylome data to a data matrix compatible with [*Methylica*](https://github.com/HiromitsuAraki/Methylica). 

## Install/Launch MatrixMaker
1.  Install [R environment](https://www.r-project.org/)
2.  Install [shiny](https://shiny.rstudio.com).  
`install.packages("shiny")`
3.  Launch *MatrixMaker*  
The following R code will launch *MatrixMaker*.  
`shiny::runGitHub("HiromitsuAraki/MatrixMaker")`
<br>

## Input file format
*MatrixMaker* accepts two types of methylome data, which are sequence-based methylome data, such as whole-genome bisulfite sequencing (WGBS) or Reduced Representation Bisulfite Sequencing (RRBS), and Infinium methylation array data (MethylationEPIC and 450k). Each data format is shown as below.
- Seuqence-basd methylome data
  - 1st column: Chr
  - 2nd column: Start
  - 3rd column: End
  - 4th column: Strand
  - 5th column: Methylated reads
  - 6th column: Total reads
  <br>
- Infinium methylation array data
  - 1st column: Target ID
  - 2nd column ~ : beta value
  
<br>

## Implementations
### Data uploading
*MatrixMaker* requires methylome data and sample metadata as its inputs. Please refer **Input file format** about the file format of methylome data and sample metadata for *Methylica*. Users need to assign the file location from a browser as below.

<img src="./README_files/Figures/MatrixMaker.png" width=500x500>
<br>

### Parameter setting
Following data upload, *Methylica* requests its users to select species with its reference genome version, genomic elements to be analyzed (CpG island, gene body, first intron and promoter), and k or the number of ICs (minimum = 2; maximum = the number of samples). *Methylica* provides a default setting of k, defined as the first k components whose cumulative contribution ratio exceeds 80% in principal component analysis. When users select all parameters, users need to press "Run" button to start analysis.
