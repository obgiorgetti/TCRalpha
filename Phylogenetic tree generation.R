# Load text files from NCBI

{
  eukaryotes = read.delim(paste0(spp.maps.folder,"eukaryotes.txt"))
  eukaryotes[order(as.numeric(eukaryotes$Size..Mb.),decreasing = T),][1:10,]
  eukaryotes[is.na(as.numeric(eukaryotes$Size..Mb.)),]
  eukaryotes.animals = eukaryotes[eukaryotes$Group %in% "Animals" & !(eukaryotes$SubGroup %in% c("Other Animals","Flatworms","Insects","Roundworms")),]
  
  overview = read.delim(paste0(spp.maps.folder,"overview.txt"))
  overview$Size..Mb. = as.numeric(overview$Size..Mb.)
  overview$Chrs = as.numeric(overview$Chrs)
  overview[is.na(overview$Size..Mb.),]
  overview.animals = overview[!is.na(overview$Chrs) & overview$Group == "Animals" & !(overview$SubGroup %in% c("Other Animals","Flatworms","Insects","Roundworms")),]
  names(overview.animals)[1] = 'Organism_name'
}


taxa <- tnrs_match_names(names = unique(as.character(sapply(spp.classification,function(x)x[x$rank=='species',1]))))
taxon_map <- structure(taxa$search_string, names = taxa$unique_name)
tree <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])
taxa.red <- tnrs_match_names(names = as.character(sapply(spp.classification,function(x)x[x$rank=='species',1])[grep("bichir|Zebrafi|Brownbanded bamb|Rat|cobra|Chicken|West af|Leishan spiny toad",names(spp.classification),value = F)]))
taxon_map.red <- structure(taxa.red$search_string, names = taxa.red$unique_name)
tree.red <- tol_induced_subtree(ott_id(taxa.red)[is_in_tree(ott_id(taxa.red))])
spp.cl.translate = sapply(spp.classification,function(x)x[x[,2]=="species",1])
tree.tip.translate = sub(" \\(.*","",gsub("_"," ",sub("_ott[0-9].*","",tree$tip.label)))
tree.tip.translate[!tree.tip.translate%in%spp.cl.translate]
tree$tip.label[!tree.tip.translate%in%spp.cl.translate]
tree.v2 = tree
tree.v2$tip.label = tree.tip.translate
pruned_tree = drop.tip(tree.v2, tree.v2$tip.label[tree.v2$tip.label %in% spp.df$spp[is.na(spp.df$J_frame)]|!tree.v2$tip.label %in% spp.df$spp])
pruned.tree.tip.translate = sub(" \\(.*","",gsub("_"," ",sub("_ott[0-9].*","",pruned_tree$tip.label)))
col.pt = rep("black",nrow(pruned_tree$edge))
for (i in 1:length(levels(spp.df$group)))col.pt[which.edge(pruned_tree,sort(match(spp.df$spp[spp.df$group == levels(spp.df$group)[i]],pruned.tree.tip.translate)))] = spp.col.xt[[i]]
