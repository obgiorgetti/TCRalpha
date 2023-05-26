# This code was used to extract the J sequences, it will not be run in the pipeline described in the readme file, since output was included in Databases.
# To reproduce this code, the complete list of C alpha loci described in Supplementary table 4 are needed.

source(paste0(spp.maps.folder,'C alpha database.R'))

Alpha.files = list.files(path= spp.maps.folder ,pattern="alpha",recursive = TRUE)
Alpha.files.ss = gsub("_"," ",gsub(".[0-9]|.*/|.alpha.*|_H[DNS].*","",Alpha.files)) # file names to match with Calphas
spp.classification = spp.classification[names(spp.classification) %in% Alpha.files.ss]
sequences.V = do.call(c,unname(lapply(c("MF","ZF","RT","MM","LA"),function(x)subseq(dictionary.extractor(spp.genome.data.list[[x]][["a"]],spp.TCR.dict[[x]][["a"]]$V.GR.dict)[[1]],-120))))
sequences.J = do.call(c,unname(lapply(c("MF","ZF","RT","MM","LA"),function(x)get.J.Seq(x,"a"))))

spp.gene.map = list()
spp.VJ.sequences = list()
spp.taxonomy = list()
spp.gene.chart = list()
spp.gene.map.width = list()
vertebrates = which(match(Alpha.files.ss,names(which(spp.classification.filter("Vertebrata"))),nomatch = 0)>0)
genome.indices = vertebrates

gmap = sapply(genome.indices,function(i)list(i))

{
  inicio = Sys.time()
  gmap[which(sapply(gmap,length) <= 1)]=mclapply(genome.indices[which(sapply(gmap,length) <= 1)],function(spp){
    {
      print(spp)
      Ja.sequences = c()  
      genome.data = NA
      my.file = Alpha.files[spp]
      taxonomy = spp.classification[[gsub("_"," ",gsub(".*/|_alpha.*|_[0-9]_alpha.*|_H[DNS].*","",my.file))]]
      print(Alpha.files[spp])
      print(taxonomy[nrow(taxonomy),])
      genome.data = readDNAStringSet(paste0("/home/giorgetti/Species Project/Species paper/Genome maps/Genome data/",my.file))
      print(names(genome.data))
      check.protein.patterns = T
      
      comienzo = Sys.time()
      Calpha.index = c()
      VJ.df = c()
      canonical_CJC = FALSE
      gpe.VJ.nt = genomic.VJ.pattern.extractor(genome.data,sequences.V,sequences.J,verbose = F,thr.V = .6,thr.J = .5,extend.V = 0,V.end.length = 200)
      
      if(sum(grepl("Ja",gpe.VJ.nt$VJ.df.summary$gene))==0){ cat("\nNo J found for",Alpha.files[spp],"\n")
        return(NULL)}
      
      #Cs    
      
      Calpha.index = which(sapply(names(Calphas.ss),function(x)grepl(sub("\\..*","",x),my.file)))
      gpl.Calpha = lapply(Calpha.index,function(Ca)genomic.pattern.locator.mk2(genome.data,sub("\\*.*","",Calphas.ss[Ca]),0,0)[[1]][c(2,4)]) # added a sub to remove stop codon in the Grayling genome
      Ca.df=lapply(gpl.Calpha,function(x)data.frame(t(c("C alpha",unlist(x),0))))
      if(length(Calpha.index>0))if(sum(sapply(1:length(Calpha.index),function(i)nrow(gpl.Calpha[[i]])>1)))Ca.df=lapply(gpl.Calpha,function(x)data.frame(cbind("C alpha",x,0)))
      Ca.df = do.call(rbind,Ca.df)
      if(!is.null(Ca.df))colnames(Ca.df) = colnames(gpe.VJ.nt$VJ.df)
      
      Cdelta.index = which(sapply(names(Cdeltas.ss),function(x)grepl(sub("delta","alpha",sub("\\..*","",x)),my.file)))
      gpl.Cdelta = lapply(Cdelta.index,function(Cd)genomic.pattern.locator.mk2(genome.data,sub("\\*.*","",Cdeltas.ss[Cd]),0,0)[[1]][c(2,4)])
      Cd.df=lapply(gpl.Cdelta,function(x)data.frame(t(c("C delta",unlist(x),0))))
      if(length(Cdelta.index>0))if(sum(sapply(1:length(Cdelta.index),function(i)nrow(gpl.Cdelta[[i]])>1)))Cd.df=lapply(gpl.Cdelta,function(x)data.frame(cbind("C delta",x,0)))
      Cd.df = do.call(rbind,Cd.df)
      if(!is.null(Cd.df))colnames(Cd.df) = colnames(gpe.VJ.nt$VJ.df)
      
      Cigh.index = which(sapply(names(Cigh.ss),function(x)grepl(sub("igh","alpha",sub("\\..*","",x)),my.file)))
      gpl.Cigh = if(length(Cigh.index)>0)genomic.pattern.locator.mk2(genome.data,paste(Cigh.ss[Cigh.index],collapse = "|"),0,0)[[1]][c(2,4)]
      Cigh.df = data.frame((cbind("C Heavy",gpl.Cigh,0)))
      if(!is.null(gpl.Cigh))colnames(Cigh.df) = colnames(gpe.VJ.nt$VJ.df)
      VJC.df = rbind(gpe.VJ.nt$VJ.df,Ca.df,Cd.df)
      
      if(check.protein.patterns ==TRUE){
        Va.pattern = "\\w{63}D[STA]A\\wY\\wC|\\w{63}DSGTY\\wC"
        Ja.pattern = "\\w{3}FG\\wGT\\w[LV]\\w[VI]\\w{2}"
        gpl.Ja = genomic.pattern.locator.mk2(genome.data,Ja.pattern,155,45)
        gpl.Va = genomic.pattern.locator.mk2(genome.data,Va.pattern,0,70*3)
        VJC.df.aa= c()
        Va.aa.F = c()
        Va.aa.seq.F = c()
        Va.aa.R = c()
        Va.aa.seq.R = c()
        Ja.aa.F = c()
        Ja.aa.seq.F = c()
        Ja.aa.R = c()
        Ja.aa.seq.R = c()
        
        if(!is.na(gpl.Va)[1]){ 
          VJC.df.aa = cbind("Va",gpl.Va[[1]][[2]],gpl.Va[[1]][[4]],0)[!gpl.Va[[1]][[2]]%in%VJC.df$position,] # doesn't take into account the unlikely possiblilty both J and V are found at same position
          
          Va.aa.F = gpl.Va[[1]][[2]][!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="F"]
          if(length(Va.aa.F)>0)gpe.VJ.nt$VJ.index.list$V.indexes.f =c(gpe.VJ.nt$VJ.index.list$V.indexes.f, Va.aa.F)
          Va.aa.R = gpl.Va[[1]][[2]][!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="R"]
          if(length(Va.aa.R)>0)gpe.VJ.nt$VJ.index.list$V.indexes.r = c(gpe.VJ.nt$VJ.index.list$V.indexes.r,Va.aa.R) 
          
          Va.aa.seq.F = (gpl.Va[[2]][!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="F"] )
          names(Va.aa.seq.F) = Va.aa.F
          if(length(Va.aa.seq.F)>0){gpe.VJ.nt$VJ.sequences.list$V.sequences.f = DNAStringSet(c(gpe.VJ.nt$VJ.sequences.list$V.sequences.f, Va.aa.seq.F))
          Va.aa.end.seq.F = do.call(c,sapply(which(!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="F"),function(i)subseq(genome.data,start = as.numeric(gpl.Va[[1]][[2]][i])+210 -3,width = 200)))
          names(Va.aa.end.seq.F) = Va.aa.F
          gpe.VJ.nt$VJ.sequences.list$V.end.sequences.f = DNAStringSet(c(DNAStringSet(gpe.VJ.nt$VJ.sequences.list$V.end.sequences.f), Va.aa.end.seq.F))
          }  
          
          Va.aa.seq.R = (gpl.Va[[2]][!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="R"])
          names(Va.aa.seq.R) = Va.aa.R
          if(length(Va.aa.seq.R)>0){
            gpe.VJ.nt$VJ.sequences.list$V.sequences.r = DNAStringSet(c(gpe.VJ.nt$VJ.sequences.list$V.sequences.r, Va.aa.seq.R))
            Va.aa.end.seq.R = reverseComplement(do.call(c,sapply(which(!gpl.Va[[1]][[2]]%in%VJC.df$position & gpl.Va[[1]][[4]]=="R"),function(i)subseq(genome.data,end = as.numeric(gpl.Va[[1]][[2]][i]) +2,width = 200))))
            names(Va.aa.end.seq.R) = Va.aa.R
            gpe.VJ.nt$VJ.sequences.list$V.end.sequences.r = DNAStringSet(c(DNAStringSet(gpe.VJ.nt$VJ.sequences.list$V.end.sequences.r), Va.aa.end.seq.R))
            
          }
        }
        
        if(!is.na(gpl.Ja)[1]){
          VJC.df.aa = rbind(VJC.df.aa,cbind("Ja",gpl.Ja[[1]][[2]],gpl.Ja[[1]][[4]],0)[!gpl.Ja[[1]][[2]]%in%VJC.df$position,])
          
          Ja.aa.F = gpl.Ja[[1]][[2]][!gpl.Ja[[1]][[2]]%in%VJC.df$position & gpl.Ja[[1]][[4]]=="F"]
          if(length(Ja.aa.F)>0)gpe.VJ.nt$VJ.index.list$J.indexes.f =c(gpe.VJ.nt$VJ.index.list$J.indexes.f, Ja.aa.F)
          Ja.aa.R = gpl.Ja[[1]][[2]][!gpl.Ja[[1]][[2]]%in%VJC.df$position & gpl.Ja[[1]][[4]]=="R"]
          if(length(Ja.aa.R)>0)gpe.VJ.nt$VJ.index.list$J.indexes.r = c(gpe.VJ.nt$VJ.index.list$J.indexes.r,Ja.aa.R)  
          
          Ja.aa.seq.F = (gpl.Ja[[2]][!gpl.Ja[[1]][[2]]%in%VJC.df$position & gpl.Ja[[1]][[4]]=="F"] )
          names(Ja.aa.seq.F) = Ja.aa.F
          if(length(Ja.aa.seq.F)>0)gpe.VJ.nt$VJ.sequences.list$J.sequences.f = DNAStringSet(c(gpe.VJ.nt$VJ.sequences.list$J.sequences.f, Ja.aa.seq.F))
          Ja.aa.seq.R = (gpl.Ja[[2]][!gpl.Ja[[1]][[2]]%in%VJC.df$position & gpl.Ja[[1]][[4]]=="R"])
          names(Ja.aa.seq.R) = Ja.aa.R
          if(length(Ja.aa.seq.R)>0)gpe.VJ.nt$VJ.sequences.list$J.sequences.r = DNAStringSet(c(gpe.VJ.nt$VJ.sequences.list$J.sequences.r, Ja.aa.seq.R))
        }
        
        if(!is.null(dim(VJC.df.aa)))colnames(VJC.df.aa) = colnames(VJC.df)
        VJC.df = rbind(VJC.df,VJC.df.aa)
        VJ.df.s = VJ.df.summarizer(VJC.df,nchar(genome.data))
        VJ.sc = VJ.score(gpe.VJ.nt,VJ.df.s$VJ.df ,sequences.V, sequences.J)  
        VJ.df.s.2 = VJ.df.summarizer(VJ.df.s$VJ.df[VJ.sc$in_frame,],nchar(genome.data))
        VJ.seqs.2 = lapply(gpe.VJ.nt$VJ.sequences.list,function(x)x[names(x)%in%VJ.df.s.2$VJ.df$position])
      }
      
      print(VJ.df.s.2$VJ.df.summary)
      print(round(difftime(Sys.time(),comienzo),1))
      
      cat("\n")
    }
    return(list(VJ.df.s.2,VJ.seqs.2,taxonomy,nchar(genome.data),names(genome.data)))
  },mc.cores = 100)
  print(difftime(Sys.time(),inicio))
}

which(sapply(gmap,length) <= 1)
{
  spp.gene.map = lapply(gmap[sapply(gmap,length)==5],function(x)x[[1]]$VJ.df.summary)
  names(spp.gene.map) = Alpha.files[genome.indices[sapply(gmap,length)==5]]
  spp.gene.chart = lapply(gmap[sapply(gmap,length)==5],function(x)x[[1]]$VJ.df)
  names(spp.gene.chart) = Alpha.files[genome.indices[sapply(gmap,length)==5]]
  
  spp.VJ.sequences = lapply(gmap[sapply(gmap,length)==5],function(x)x[[2]])
  names(spp.VJ.sequences)= Alpha.files[genome.indices[sapply(gmap,length)==5]]
  spp.taxonomy= lapply(gmap[sapply(gmap,length)==5],function(x)x[[3]])
  names(spp.taxonomy)= Alpha.files[genome.indices[sapply(gmap,length)==5]]
  
  spp.gene.map.width = lapply(gmap[sapply(gmap,length)==5],function(x)x[[4]])
  names(spp.gene.map.width)= Alpha.files[genome.indices[sapply(gmap,length)==5]]
  spp.gene.names = lapply(gmap[sapply(gmap,length)==5],function(x)x[[5]])
  names(spp.gene.names)= Alpha.files[genome.indices[sapply(gmap,length)==5]]
  
}

spp.gm.translate = gsub("_"," ",gsub(".*/|_alpha.*|_[0-9]","",names(spp.gene.chart))) # move to after for loop 1

spp.Ca = sapply(spp.gene.chart,function(x)grep("C alpha",x$segment))
spp.Cd = sapply(spp.gene.chart,function(x)grep("C delta",x$segment))
selected.Js = list()
C.ranges = list()
for(spp in names(spp.gene.chart)){
  {
    CaF = which(spp.gene.chart[[spp]]$segment=="C alpha" & spp.gene.chart[[spp]]$sense == "F")
    CdF = which(spp.gene.chart[[spp]]$segment=="C delta" & spp.gene.chart[[spp]]$sense == "F")
    CaR = which(spp.gene.chart[[spp]]$segment=="C alpha" & spp.gene.chart[[spp]]$sense == "R")
    CdR = which(spp.gene.chart[[spp]]$segment=="C delta" & spp.gene.chart[[spp]]$sense == "R")
    
    selected.Js.F = c()
    selected.Js.R = c()
    C.range = c()
    
    if(length(CaF) >0 & length(CdF)>0){
      CadF.indices = lapply(CaF,function(x)if(sum(CdF<x)>0)c(as.numeric(spp.gene.chart[[spp]][x,2]),as.numeric(spp.gene.chart[[spp]][max(CdF[CdF<x]),2])))
      if(length(CadF.indices)>1){CadF.sel = tapply(do.call(rbind,CadF.indices)[,1],do.call(rbind,CadF.indices)[,2],min)
      CadF.indices = lapply(1:length(CadF.sel),function(x)as.numeric(c(CadF.sel[x],names(CadF.sel[x]))))}
      
      selected.Js.F = lapply(CadF.indices,function(x)if (length(x)> 0) spp.VJ.sequences[[spp]]$J.sequences.f[as.numeric(names(spp.VJ.sequences[[spp]]$J.sequences.f)) %in% x[[1]]:x[[2]]])
      selected.Js.F = selected.Js.F[sapply(selected.Js.F,length)>0]
      selected.Js.F = lapply(selected.Js.F,function(x)x[order(as.numeric(names(x)))])
      C.range = c(C.range,CadF.indices)
    }
    
    if(length(CaR) >0 & length(CdR)>0){
      CadR.indices = lapply(CaR,function(x)if(sum(CdR>x)>0)c(as.numeric(spp.gene.chart[[spp]][x,2]),as.numeric(spp.gene.chart[[spp]][min(CdR[CdR>x]),2])))
      if(length(CadR.indices)>1){CadR.sel = tapply(do.call(rbind,CadR.indices)[,1],do.call(rbind,CadR.indices)[,2],max)
      CadR.indices = lapply(1:length(CadR.sel),function(x)as.numeric(c(CadR.sel[x],names(CadR.sel[x]))))}
      selected.Js.R = lapply(CadR.indices,function(x)if (length(x)> 0) spp.VJ.sequences[[spp]]$J.sequences.r[as.numeric(names(spp.VJ.sequences[[spp]]$J.sequences.r)) %in% x[[1]]:x[[2]]])
      selected.Js.R = selected.Js.R[sapply(selected.Js.R,length)>0]
      selected.Js.R = lapply(selected.Js.R,function(x)x[order(as.numeric(names(x)))])
      C.range = c(C.range,CadR.indices)
      
    }
    
    selected.Js[[spp]] = append(selected.Js.F,selected.Js.R)
    C.ranges[[spp]] = C.range
    
  }
}



C.ranges.df = sapply(C.ranges,function(x)gsub('-','',paste(sapply(x,diff),collapse = ',')))

selected.Js.2 = list()
for(spp in names(spp.gene.chart)){
  {
    CaF = which(spp.gene.chart[[spp]]$segment=="C alpha" & spp.gene.chart[[spp]]$sense == "F")
    CdF = which(spp.gene.chart[[spp]]$segment=="C delta" & spp.gene.chart[[spp]]$sense == "F")
    CaR = which(spp.gene.chart[[spp]]$segment=="C alpha" & spp.gene.chart[[spp]]$sense == "R")
    CdR = which(spp.gene.chart[[spp]]$segment=="C delta" & spp.gene.chart[[spp]]$sense == "R")
    
    selected.Js.F = c()
    selected.Js.R = c()
    if(length(CaF) >0){
      CadF.indices = lapply(CaF,function(x)c(as.numeric(spp.gene.chart[[spp]][x,2])))
      selected.Js.F = lapply(CadF.indices,function(x)if (length(x)> 0) spp.VJ.sequences[[spp]]$J.sequences.f[as.numeric(names(spp.VJ.sequences[[spp]]$J.sequences.f)) %in% (x[[1]]- 1000000):(x[[1]])])
      selected.Js.F = selected.Js.F[sapply(selected.Js.F,length)>0]
      selected.Js.F = lapply(selected.Js.F,function(x)x[order(as.numeric(names(x)))])
      selected.Js.F = unique(do.call(c,selected.Js.F))
    }
    
    if(length(CaR) >0){
      CadR.indices = lapply(CaR,function(x)c(as.numeric(spp.gene.chart[[spp]][x,2])))
      selected.Js.R = lapply(CadR.indices,function(x)if (length(x)> 0) spp.VJ.sequences[[spp]]$J.sequences.r[as.numeric(names(spp.VJ.sequences[[spp]]$J.sequences.r)) %in% x[[1]]:(x[[1]] +1000000)])
      selected.Js.R = selected.Js.R[sapply(selected.Js.R,length)>0]
      selected.Js.R = lapply(selected.Js.R,function(x)x[order(as.numeric(names(x)))])
      selected.Js.R = unique(do.call(c,selected.Js.R))
      
    }
    
    selected.Js.2[[spp]] = append(selected.Js.F,selected.Js.R)
  }
}



if(find.RSS == T){
  
  {
    spp.RSS = list()
    spp.RSS.xt = list()
    spp.RSS.MM = list()
    spp.RSS.2 = list()
    spp.RSS.xt.2 = list()
    spp.RSS.MM.2 = list()
  }
  inicio = Sys.time()
  spp.RSS[sapply(selected.Js,length)>0] = mclapply(which(sapply(selected.Js,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("ZF","a"),selected.Js[[i]][[1]],iter.max = 5),mc.cores = 100)
  spp.RSS.xt[sapply(selected.Js,length)>0] = mclapply(which(sapply(selected.Js,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("ZF","a",extend = 5),selected.Js[[i]][[1]],iter.max = 5),mc.cores = 100)
  spp.RSS.MM[sapply(selected.Js,length)>0] = mclapply(which(sapply(selected.Js,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("MM","a"),selected.Js[[i]][[1]],iter.max = 5),mc.cores = 100)
  spp.RSS.2[sapply(selected.Js.2,length)>0] = mclapply(which(sapply(selected.Js.2,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("ZF","a"),selected.Js.2[[i]],iter.max = 5),mc.cores = 100)
  spp.RSS.xt.2[sapply(selected.Js.2,length)>0] = mclapply(which(sapply(selected.Js.2,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("ZF","a",extend = 5),selected.Js.2[[i]],iter.max = 5),mc.cores = 100)
  spp.RSS.MM.2[sapply(selected.Js.2,length)>0] = mclapply(which(sapply(selected.Js.2,length)>0),function(i)RSS.scorer.iterative(get.J.RSS.matrix("MM","a"),selected.Js.2[[i]],iter.max = 5),mc.cores = 100)
  names(spp.RSS) = names(sapply(selected.Js,length)>0)
  names(spp.RSS.xt) = names(sapply(selected.Js,length)>0)
  
  print(difftime(Sys.time(),inicio))
  
  save(spp.RSS, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS')
  save(spp.RSS.xt, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.xt')
  save(spp.RSS.MM, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.MM')
  save(spp.RSS.2, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.2')
  save(spp.RSS.xt.2, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.xt.2')
  save(spp.RSS.MM.2, file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.MM.2')
  
  
} else {
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS')
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.xt')
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.MM')
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.2')
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.xt.2')
  load( file = '/data/boehm/group/giorgetti/Species_project/RSS/spp.RSS.MM.2')
}
spp.cl.translate = sapply(spp.classification,function(x)x[x[,2]=="species",1])

{
  spp.df = data.frame(spp = spp.cl.translate,common.name = names(spp.classification))
  spp.df$Ca_nr = NA
  spp.df$Cd_nr = NA
  spp.df$J_nr = NA
  spp.df$J_range = NA
  spp.df$C.ranges = NA
  spp.df$J_frame = NA
  spp.df$J_ec = NA
  spp.df$group = NA
  spp.df$genome.size = NA
  spp.df$read.width = NA
  
  spp.df$J_nr.2 = NA
  spp.df$J_ec.2 = NA
  spp.df$J_frame.2 = NA
  spp.df$J_frame.MM = NA
  spp.df$J_frame.xt = NA
  spp.df$J_frame.2 = NA
  spp.df$J_frame.MM.2 = NA
  spp.df$J_frame.xt.2 = NA
  spp.df$file.name = NA
  spp.df$file.folder = NA
  
  spp.df$J_nr.2[match(Alpha.files.ss[match(names(selected.Js.2),Alpha.files)],spp.df$common.name)] = sapply(selected.Js.2,length)
  
  
  selected.Js.ec = sapply(spp.RSS.xt,function(x)if(length(x)>0)sum(entropy.content(consensusMatrix(x$ss.RSS)[c("A","C","G","T"),29:32]))else NA)
  spp.df$J_ec[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = selected.Js.ec
  
  selected.Js.ec.2 = sapply(spp.RSS.xt,function(x)if(length(x)>0)sum(entropy.content(consensusMatrix(x$ss.RSS)[c("A","C","G","T"),]))else NA)
  spp.df$J_ec.2[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = selected.Js.ec.2
  
  
  spp.df$J_frame[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS),function(i)mean(spp.RSS[[i]]$RSS.index%%3==0))
  spp.df$J_frame.2[match(Alpha.files.ss[match(names(selected.Js.2),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS.2),function(i)mean(spp.RSS.2[[i]]$RSS.index%%3==0))
  
  spp.df$J_frame.MM[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS.MM),function(i)mean(spp.RSS.MM[[i]]$RSS.index%%3==0))
  spp.df$J_frame.xt[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS.xt),function(i)mean(spp.RSS.xt[[i]]$RSS.index%%3==0))
  
  spp.df$J_frame.2[match(Alpha.files.ss[match(names(selected.Js.2),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS.2),function(i)mean(spp.RSS.2[[i]]$RSS.index%%3==0))
  
  spp.df$J_frame.MM.2[match(Alpha.files.ss[match(names(selected.Js.2),Alpha.files)],spp.df$common.name)][sapply(spp.RSS.MM.2,length)==4] = sapply(which(sapply(spp.RSS.MM.2,length)==4),function(i)mean(spp.RSS.MM.2[[i]]$RSS.index%%3==0))
  spp.df$J_frame.xt.2[match(Alpha.files.ss[match(names(selected.Js.2),Alpha.files)],spp.df$common.name)] = sapply(1:length(spp.RSS.xt.2),function(i)mean(spp.RSS.xt.2[[i]]$RSS.index%%3==0))
  
  
  spp.df$group = factor(spp.df$group, levels = c("Cartilaginous fishes","Basal fishes","Teleosts","Transition","Amphibians","Reptiles","Birds","Mammals") )
  spp.df$group[spp.classification.filter("Teleostei")] = "Teleosts"
  spp.df$group[spp.classification.filter("Actinopterygii") & !spp.classification.filter("Teleostei")] = "Basal fishes"
  spp.df$group[spp.classification.filter("Chondrichthyes")] = "Cartilaginous fishes"
  spp.df$group[spp.classification.filter("Mammalia")] = "Mammals"
  spp.df$group[spp.classification.filter("Sauropsida")& !spp.classification.filter("Aves")] = "Reptiles"
  spp.df$group[spp.classification.filter("Aves")] = "Birds"
  spp.df$group[spp.classification.filter("Amphibia")] = "Amphibians"
  spp.df$group[spp.classification.filter("Coelacanth|Dipnomorpha")] = "Transition"
  spp.df$genome.size = as.numeric(eukaryotes$Size..Mb.[match(spp.df$spp,eukaryotes$X.Organism.Name)])
  spp.df$read.width[match(Alpha.files.ss[match(names(spp.gene.map.width),Alpha.files)],spp.df$common.name)] = as.numeric(spp.gene.map.width)
  spp.df$file.name[match(Alpha.files.ss[match(names(spp.gene.names),Alpha.files)],spp.df$common.name)] = as.character(spp.gene.names)
  spp.df$file.folder[match(Alpha.files.ss[match(names(spp.gene.names),Alpha.files)],spp.df$common.name)] = names(spp.gene.names)
  spp.df$Ca_nr = sapply( spp.gene.map,function(x)sum(grepl("C alpha",x$gene)))[spp.df$file.folder]
  spp.df$Cd_nr = sapply( spp.gene.map,function(x)sum(grepl("C delta",x$gene)))[spp.df$file.folder]
  spp.df$J_nr[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = sapply(selected.Js,function(x)length(do.call(c,unname(x))))
  spp.df$J_range[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = sapply(1:length(selected.Js),function(i)if(length( selected.Js[[i]])>0)diff(range(as.numeric(names(selected.Js[[i]][[1]]))))else NA)
  
  spp.df$C.ranges[match(Alpha.files.ss[match(names(selected.Js),Alpha.files)],spp.df$common.name)] = C.ranges.df
  spp.df$accession = eukaryotes.animals[match(spp.df$spp,eukaryotes.animals$X.Organism.Name),]$Assembly.Accession
  
  spp.accession.nas = sapply(spp.df[is.na(spp.df$accession),][,1],function(x)grep(x,eukaryotes.animals$X.Organism.Name,value = T))
  spp.df$accession[is.na(spp.df$accession )][ sapply(spp.accession.nas,length)==1] = eukaryotes.animals$Assembly.Accession[match(spp.accession.nas[ sapply(spp.accession.nas,length)==1],eukaryotes.animals$X.Organism.Name)]
  mm.acc = match(lapply(strsplit(spp.df$spp[is.na(spp.df$accession )]," "),function(x)paste(x[1],x[2],x[2])),eukaryotes.animals$X.Organism.Name,nomatch = 0)
  spp.df$accession[is.na(spp.df$accession )][mm.acc > 0 ] = eukaryotes.animals$Assembly.Accession[mm.acc]
  spp.df$accession[spp.df$common.name == "Deer mouse"] = 'GCA_003704035.3'
  spp.df$accession[spp.df$common.name == "Minifish"] = 'JABUMV000000000'
  spp.df$accession[spp.df$common.name == "Arctic char"] = 'GCA_002910315.2' # Alpinus missing in eukariotes.txt, mistake on their part
  spp.df$accession[spp.df$common.name == "African ostrich"] = 'GCA_000698965.1' # For mismatch Struthio camelus / Struthio camelus australis
  spp.df$accession[spp.df$common.name == "Golden eagle"] = 'GCA_900496995.2' # Aquila chrysaetos / Aquila chrysaetos chrysaetos 
  spp.df$accession[spp.df$common.name == "Arctic char"] = 'GCA_002910315.2' # Alpinus missing in eukariotes.txt, mistake on their part
  spp.df$accession[spp.df$common.name == "Prairie rattlesnake"] = 'GCA_003400415.2' # Crotalus viridis is Crotalus viridis viridis 
  spp.df$accession[spp.df$common.name == "American flamingo"] = 'GCA_009819775.1 ' # Phoenicopterus ruber / Phoenicopterus ruber ruber
  spp.df$accession[spp.df$common.name == "Japanese medaka"] = 'GCA_002234695.1'
  spp.df$accession[spp.df$common.name == "German shepherd"] = 'GCF_011100685.1'
}

#spp.classification$`Asian redtail catfish` = classification("Hemibagrus wyckioides", db = 'ncbi')[[1]]


