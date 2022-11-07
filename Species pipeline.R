spp.root.folder = "/Users/giorgetti/Desktop/Species_project"
{
  spp.code.folder = file.path(spp.root.folder,"Pipeline/")
  spp.data.folder = file.path(spp.root.folder,"Databases/")
  spp.maps.folder = file.path(spp.root.folder,"Genome maps/")
  spp.tables.folder = file.path(spp.root.folder,"Tables/")
  spp.fig.folder = file.path(spp.root.folder,"Figures/")

  spp.short.names = c("C. punctatum","P. senegalus","A. ruthenus Ca 1","A. ruthenus Ca 2","P. progenetica","D. rerio","O. mykiss","O. latipes","P. annectens","M. musculus","L. africana")

  # Runs in around 10 minutes
  inicio = Sys.time()
  source(paste0(spp.code.folder,"Data extraction functions.R"))
  source(paste0(spp.code.folder,"Sequence manipulation functions.R"))
  # Generation of original parsed files is in VDJ pipeline.R
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

source(file.path(spp.code.folder,"/Figures code/Fig_1.R"))
source(file.path(spp.code.folder,"/Figures code/Fig_2.R"))
source(file.path(spp.code.folder,"/Figures code/Fig_3.R"))
source(file.path(spp.code.folder,"/Figures code/Fig_4.R"))
source(file.path(spp.code.folder,"/Figures code/Fig_Sup"))