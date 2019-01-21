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
*MatrixMaker* accepts two types of methylome data, which are sequence-based methylome data, such as whole-genome bisulfite sequencing (WGBS) or Reduced Representation Bisulfite Sequencing (RRBS), and Infinium methylation array data (MethylationEPIC and 450k). Users for seuqence-basd methylome data need to prepare BED format files of each sample generated by general methylation calling tools, such as [BSMAP](https://www.ncbi.nlm.nih.gov/pubmed/19635165) or [Bismark](https://www.ncbi.nlm.nih.gov/pubmed/21493656). Users for infinium methylation array need to prepare *Methylation Profile table* generated by  [GenomeStudio](http://jp.support.illumina.com/array/array_software/genomestudio.html). Each data format is shown as below.
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
### Data uploading and parameter setting
*MatrixMaker* requires methylome data as its inputs. Please refer **Input file format** about the file format of methylome data. First, users need to assign a directory location in which only BED format files or *Methylation Profile table* are stored. Following the assigment of directory location, *MatrixMaker* requests its users to choose platform users used, species with its reference genome version used in sequence alignment, and user's interesting genomic elements (CpG island, gene body, first intron and promoter). When users choose promoter as genomic element, users can decide promoter regions from -10,000 to -1 and from 0 to 10,000 as promoter start and end position respectively, where the zero position represents transcription start site (TSS). Finaly, data matrix can be downloaded  
need to press "Submit" button to generate methylome data matrix for *Methylica*.

<img src="./README_files/Figures/MatrixMaker.png" width=500x500>
<br>
