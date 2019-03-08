# ReFernment
[![Build Status](https://travis-ci.com/TARobison/ReFernment.svg?token=q2xtkzWBws7Qp6deSQsN&branch=master)](https://travis-ci.com/TARobison/ReFernment)

## Why ReFernment
ReFernment was created in response to the drudgery of manually annotating a plastome sequence that has levels of RNA editing. These editing sites result in a genomic sequence that contains features that look like errors (e.g., internal stop codons) so the annotations will be rejected by NCBI unless they are annotated correctly. ReFernment corrects the annotation to produce files ready for GenBank submission. ReFernment takes as input a GenBank flat file (.gb) and a GFF3 file (no sequence) and generates new annotations for nonsense mutations that can reasonably be explained by RNA editing, and provides conceptual translations for coding sequences with RNA editing. 
It should be made clear that ReFernment is not intended to predict every RNA editing site, as can other available software packages. These packages rely on RNA-Seq data to compare the genomic DNA, and thus determine if a nucleotide has been edited, but these data are not always available to researchers. We made ReFernment as a simple tool to save time and ease the GenBank submission process for those working with plastids that have high levels of RNA editing. 

## How it works and how to cite
Please see a breif description of how ReFernment operates in our short [Software Note](https://doi.org/10.1002/aps3.1216) in Applications in Plant Sciences. Please, if you use ReFernment for published work, cite the above manuscript.  

## How to use it
ReFernment *requires* only a GFF3 file (with sequence embedded) to work, but you can optionally provide a GenBank flat file and seperate fasta file, if your workflow requires this. ReFernment takes as input the following:

* `gffFolderPath`(required): directory containing the GFF files 
* `outputFolderPath`(required): the directory you would like to output the newly generated annotation files
* `genomes`(required): the names of the files being used
* `gbFolderPath`(optional): directory containing the GB files
* `fastaFolderPath`(optional): directory containing the FASTA files

The `gbFolderPath` and `gffFolderPath` can refer to the the same directory, but unless you are *totally comfortable* with ReFernment overwriting the original GB files, then `outputFolderPath` should refer to a directory separate from the other two. When preparing to use ReFernment you should keep the names of the GB files and GFF files identical, except for their extensions. For example `Hemionitis_subcordata.gb` and `Hemionitis_subcordata.gff`. This makes it easy for ReFernment to loop through many files all at once if you have large numbers of plastomes that need to be annotated. Along those lines, `genomes` should **not** contain the file extension. ReFernment will add the file extensions later.

Below we want ReFernment to annotate the plastomes of Asplenium pekinense and Woodwardia unigemmata, which are both found in the `examples` folder. To start, we open an R console (either in your preferred command line interface or in RStudio), load ReFernment into the workspace. Then we declare a vector named `genomes` which contains all of the plastomes which you would like to annotate. 

```r
install.packages("path/to/ReFernment/")
#add the above command if this is your first time using ReFernment
library("ReFernment")
genomes <- c("Asplenium_pek", "Woodwardia_uni")
```
Next, we'll declare a vector containing the path to the input GFF3 files and another vector to the location where we would like the output files to be saved. 

```r
gffFolder <- "/Users/Me/Location/Of/gff/Files/"

outputFolder <- "/Users/Me/Location/Of/Output/Folder/"
```

Next, we simply call the `ReFernment` function and wait for it to finish. This can take several minutes if you have a large number of plastomes. Note that the operation of ReFernment is the same if you are also providing gb and/or fasta files, just add the required variables. 

```r
ReFernment(gffFolder, fastaFolderPath=NULL, gbFolderPath=NULL, outputFolder, genomes)
```
As ReFernment runs, it may produce one or more of the following warning messages:

```
* There are a high number of edited Stops ( [number] ) in [geneName] manually check to make sure frame is correct
* There seems to be a problem with the protein sequence of ( [protein name] ) please check its sequence in the protein fasta manually
```

The first message is intended to warn the user of two common problems:

1. potential assembly errors that would result in an entire gene appearing to be frameshifted
2. annotation errors where the proper frame for a gene was not selected to begin with

The warning will appear if a given coding sequence has more than 5 internal stops. There are, of course, cases where this happens and it is not the result of assembly error, but these are relatively rare. Especially if there are more than 10 internal stops the user should take a close look at this coding sequecne to check for assembly error or annotation errors. 

The second message is produced when there are errors translating a sequence. This is usually due to the coding sequence spanning the orgin of the genome. I suggest that you shift the origin of the genome so that it doesn't have a CDS overlapping it. If the problem persists after doing this, please feel free to contact me. 

## Getting Help
If  you have any issues using ReFernment, please get in touch with me and I will do my best to address them. You can either do this by eithter creating an issue on GitHub, or emailing me directly at 

robison dot tanner at gmail dot com. 

Likewise, if you have any suggestions for improving ReFernment, I would love to hear them. I want to make this tool as useful as possible. 
