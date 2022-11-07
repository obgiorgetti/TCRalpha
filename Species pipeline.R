{
spp.code.folder = "/Users/giorgetti/Desktop/Species_project/Pipeline/"
spp.data.folder = "/Users/giorgetti/Desktop/Species_project/Databases/"
spp.maps.folder = "/Users/giorgetti/Desktop/Species_project/Genome maps/"
spp.tables.folder = "/Users/giorgetti/Desktop/Species_project/Tables/"
spp.fig.folder = "/Users/giorgetti/Desktop/Species_project/Figures/"

spp.short.names = c("C. punctatum","P. senegalus","A. ruthenus Ca 1","A. ruthenus Ca 2","P. progenetica","D. rerio","O. mykiss","O. latipes","P. annectens","M. musculus","L. africana")

  # Runs in around 10 minuts
  inicio = Sys.time()
  source(paste0(spp.code.folder,"Data extraction functions.R"))
  source(paste0(spp.code.folder,"Sequence manipulation functions.R"))
  # Generation of original parsed files is in VDJ pipeline.R
  # To load lists from VDJ pipeline step 1:
  for(x in list.files(paste0(spp.data.folder,"220801"))) load(file = paste0(paste0(spp.data.folder,"220801/"),x))
  # To load aggregated lists:
  for(x in list.files(paste0(spp.data.folder,"RepSeq_lists"))) load(file = paste0(paste0(spp.data.folder,"RepSeq_lists/"),x))
  for(x in list.files(paste0(spp.data.folder,"RSS/"))) load(file = paste0(paste0(spp.data.folder,"RSS/"),x))
  source(paste0(spp.code.folder,"IT functions.R"))
  source(paste0(spp.code.folder,"Repertoire data frames generation.R"))
  source(paste0(spp.code.folder,"Plot functions.R"))
  source(paste0(spp.code.folder,"Phylogeny functions.R"))
  source(paste0(spp.code.folder,"Phylogenetic tree generation.R"))
  
  
  load(paste0(spp.data.folder,"RSS/selected.Js"))
  print("Tiempo total")
  print(difftime(Sys.time(),inicio))
}

source("Pipeline/Figures code/Fig_1.R")
source("Pipeline/Figures code/Fig_2.R")
source("Pipeline/Figures code/Fig_3.R")
source("Pipeline/Figures code/Fig_4.R")
source("Pipeline/Figures code/Fig_Sup.R")
