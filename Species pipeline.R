spp.root.folder = "/Users/giorgetti/R/TCRalpha-master"

{
library(tidyverse)
library(Biostrings)
library(ShortRead)
library(igraph)
library(parallel)
library(BSgenome)
library(plyr)
library(rotl)
library(taxize)
library(ape)
library(ggplot2)
library(ggseqlogo)
library(gridExtra)
library(ggrepel)
library(msa)
}
{
  inicio = Sys.time()
  spp.code.folder = file.path(spp.root.folder,"")
  spp.data.folder = file.path(spp.root.folder,"Databases/")
  spp.maps.folder = file.path(spp.root.folder,"Genome maps/")
  spp.tables.folder = file.path(spp.root.folder,"Tables/")
  if(!file.exists(spp.tables.folder))dir.create(spp.tables.folder)
  spp.fig.folder = file.path(spp.root.folder,"Figures/")
  if(!file.exists(spp.fig.folder))dir.create(spp.fig.folder)
  
  spp.short.names = c("C. punctatum","P. senegalus","A. ruthenus Ca 1","A. ruthenus Ca 2","P. progenetica","D. rerio","O. mykiss","O. latipes","P. annectens","M. musculus","L. africana")

  source(paste0(spp.code.folder,"Data extraction functions.R"))
  source(paste0(spp.code.folder,"Sequence manipulation functions.R"))
  source(paste0(spp.data.folder,"C alpha database.R"))
  # Generation of original parsed files is in VDJ pipeline.R
  # To load aggregated lists:
  TCR.a.spp <- readRDS(paste0(spp.data.folder,"TCR.a.spp.RDS"))
  TCR.b.spp <- readRDS(paste0(spp.data.folder,"TCR.b.spp.RDS"))
  spp.genome.data.list <- readRDS(paste0(spp.data.folder,"spp.genome.data.list.RDS"))
  ZF.Ca.mut.df.list <- readRDS(paste0(spp.data.folder,"ZF.Ca.mut.df.list.RDS"))
  ZF.Ca.wt.df.list <- readRDS(paste0(spp.data.folder,"ZF.Ca.wt.df.list.RDS"))
  ZF.23.25.dla.a <- readRDS(paste0(spp.data.folder,"ZF.23.25.dla.a.RDS"))
  MM.a.mut.a32.dfl.MUT <- readRDS(paste0(spp.data.folder,"MM.a.mut.a32.dfl.MUT.RDS"))
  spp.RSS <- readRDS(paste0(spp.data.folder,"RSS_data/spp.RSS.RDS"))
  
  spp.RSS.xt <- readRDS(paste0(spp.data.folder,"RSS_data/spp.RSS.xt.RDS"))
  
  for(x in list.files(paste0(spp.data.folder,"Dictionaries"),full.names = T)) load(x)
  for(x in list.files(paste0(spp.data.folder,"RSS/"))) load(file = paste0(paste0(spp.data.folder,"RSS/"),x))
  source(paste0(spp.code.folder,"IT functions.R"))
  source(paste0(spp.code.folder,"Repertoire data frames generation.R"))
  source(paste0(spp.code.folder,"Plot functions.R"))
  source(paste0(spp.code.folder,"Phylogeny functions.R"))
  source(paste0(spp.code.folder,"Phylogenetic tree generation.R"))
  
  
  load(paste0(spp.data.folder,"RSS_data/selected.Js"))
  
  
  source(file.path(spp.code.folder,"/Figures code/Fig_1.R"))
  source(file.path(spp.code.folder,"/Figures code/Fig_2.R"))
  source(file.path(spp.code.folder,"/Figures code/Fig_3.R"))
  source(file.path(spp.code.folder,"/Figures code/Fig_4.R"))
  source(file.path(spp.code.folder,"/Figures code/Fig_Sup.R"))
  print("Tiempo total")
  print(difftime(Sys.time(),inicio))
  
  }
