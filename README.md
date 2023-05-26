#############################################################

This is the repository for the paper "Origin and evolutionary malleability of T cell receptor α diversity"
by Giorgetti, O.B., O’Meara, C.P., Schorpp, M., and Boehm, T.

#############################################################

The R code will generate the figures in the main figures and extended data figures, 
from repertoire data tables containing CDR3 counts.

We also provide below a few examples on how to explore the data using custom functions.

If you require additional help, contact OBG at orlandogiorgetti@gmail.com.

#############################################################



Raw data
========

In all cases raw data consists of Illumina sequencing of TCRa and TCRb cDNA amplicons with UMI barcoding (see paper for more details).
The raw data can be found in PRJNA612865 (Minifish), PRJNA865512 (all other wild species) and PRJNA865921 (Crispr mutant animals).
For using the tables we generated for the analysis done in the paper see below.


Tables with extracted UMI counts
================================

These tables, along with the genomes, and the genes used for the extraction are located at Zenodo in a zipped file.
Unzipping the file in the folder containing this file should create a folder named "Databases" in it. Databases should contain 8 files and 2 folders.
The files TCR.a.spp and TCR.b.spp can be navigated in R after running the pipeline, and contain the extracted table for the repertoire sequencing of 10 species.


Running the code
================

1) Change the first line of code in the file "Species pipeline.R" so that spp.root.folder points to folder containing this readme, where "Databases" was previously unzipped(see above).
2) Load all the libraries listed, you might need to install some. A list of functioning versions of these packages at the moment of this writing is can be found below.
3) Run the rest of the code in "Species pipeline.R". Note internet connection is needed for the code to run, due to the "Open Tree of Life" package.
4) This code will create additional folders for the Figures and tables.
It takes around 2,5 to 3 minutes on a Mac Mini 2023 to run the code.



Exploring the tables and extracting data from the genomes for repertoire data
=============================================================================

Examples for further exploration of the data: 
  
The list TCR.a.spp.VJ was used in the final version of the paper; it contains only the molecules where it was possible to identify V and J genes.
TCR.a.spp includes all V and J, including ambiguous or not mapped variants. TCR.a.spp.rules is just a dummy variable for generating plots.
    
The dataframe with the wildtype Zebrafish individuals TCR alpha can be accessed by:
  TCR.a.spp.VJ$ZF

The label of each individual specimen is the variable Sample.bio.name, which matches the numbering of the raw data:
  TCR.a.spp.VJ$ZF$Sample.bio.name

To make a plot in the style of Fig1 c use:
  VDJ.plot(TCR.a.spp.VJ$ZF)

The list of germline J segments (genes) with GRanges format, relative to the extracted genomes is:
  spp.TCR.dict$ZF$a$J.GR.dict

The RSS for those elements, where the index is the position of the last letter of the RSS is:
  spp.TCR.dict$ZF$a$J.GR.dict$heptamer.end

Where the extracted genome is:
  spp.genome.data.list$ZF$a

To get the sequences:
  getSeq(spp.genome.data.list$ZF$a,spp.TCR.dict$ZF$a$J.GR.dict)
or just use:
  get.J.Seq("ZF","a")

Note we included allelic variants for wild species, but the sequencing for model species (D. rerio and M. musculus) was done with the IMGT reference genes.
CRISPR constructs were injected on the TLEK line, where we found some allelic variants respect to Tübingen fish. The J gene variants are:
  mut.TCR.dict$ZF_TLEK$a$J.alleles

A quick way of viewing the entropy and mutual information between CDR3 length, V and J genes:
  general.entropy.report(TCR.a.spp.VJ$ZF)
  Note the results in our paper imply there is a high mutual information between J genes and CDR3 length for TCR alpha, especially in teleosts.

To calculate the entropy of the repertoire, and generate a plot for 42 nt long CDR3s:
Note we only used this formula for TCR alpha in this paper.
  h.ZFa = TCR.entropy.VorJ(database = TCR.a.spp.VJ$ZF)
  h.ZFa$statistics
  do.call('entropy.output.VorJ',list(h.ZFa$nt[[42]],plot = 4))
  

Analysis of germline elements in genome assemblies
==================================================

The consensus matrix for the RSS of Zebrafish J alpha can be extracted as follows:
  get.J.RSS.matrix("ZF","a")
  
Where the RSS can be found as follows:  
  RSS.scorer.iterative(get.J.RSS.matrix("ZF","a"),get.J.Seq("ZF","a"),iter.max = 5)
  where RSS.index is the list of positions for the start of the matrix (note this always positive and is NOT the index we use for RSS in the case of J, corresponding to the last nucleotide of   the RSS, but one can be easily translated into the other by substracting the total length + 1 from the position)

The extraction of J genes in the vicinity of C alpha locus was described in the text and is basically a regular expression search. It is not exhaustive but it covers a majority of the J genes in our test cases. The extracted Js can be inspected by looking at:
  selected.Js

The Calpha and Cdelta elements that delimit the J elements used for the analysis can be accessed from the variables:
  Calphas
  Cdeltas

The phylogeny for the species is in:
  spp.classification

And a translation between common name and taxonomy:
  spp.cl.translate
