library(rotl)
library(taxize)
library(ape)

spp.classification.filter = function(taxonomy)sapply(spp.classification,function(x)if(!is.null(x))sum(grepl(taxonomy,x$name))>0 else F) 

VJ.df.summarizer = function(VJ.df,genome.length){
  
  VJ.df = VJ.df[order(as.numeric(VJ.df$position)),]
  VJ.df$gap = c(diff(as.numeric(VJ.df$position)),genome.length-max(as.numeric(VJ.df$position)))
  rle.df = rle(paste(VJ.df$segment,VJ.df$sense))
  VJ.df.summary = data.frame(rle.df$lengths,rle.df$values)
  VJ.df.sp = split(VJ.df,rep(1:length(rle.df$values),rle.df$length))
  median.gap = round(sapply(1:length(VJ.df.sp),function(i)if(nrow(VJ.df.sp[[i]])>1)median(VJ.df.sp[[i]]$gap[-length(VJ.df.sp[[i]]$gap)])else 0))
  max.gap = round(sapply(1:length(VJ.df.sp),function(i)if(nrow(VJ.df.sp[[i]])>1)max(VJ.df.sp[[i]]$gap[-length(VJ.df.sp[[i]]$gap)])else 0))
  start.at = sapply(1:length(VJ.df.sp),function(i)VJ.df.sp[[i]]$position[[1]])
  rango = sapply(1:length(VJ.df.sp),function(i)diff(range(as.numeric(VJ.df.sp[[i]]$position))))
  end.at = sapply(1:length(VJ.df.sp),function(i)max(as.numeric(VJ.df.sp[[i]]$position)))
  gap.next = sapply(1:length(VJ.df.sp),function(i)VJ.df.sp[[i]]$gap[length(VJ.df.sp[[i]]$gap)])
  
  VJ.df.summary = cbind(start.at,rango,end.at,gap.next,VJ.df.summary,median.gap,max.gap)
  colnames(VJ.df.summary) = c("start","range","end","gap to next", "segments","gene","median gap","max gap")
  rownames(VJ.df) = NULL
  return(list(VJ.df.summary=VJ.df.summary,VJ.df=VJ.df))
}

genomic.VJ.pattern.extractor = function(genome.data, secuencias.V, secuencias.J, thr.V=0.3, thr.J=0.3, gene = "a", V.length = 210, J.length = 200,V.end.length = 100,extend.V = 0,verbose = T){
  if(verbose ==T)cat(names(genome.data),"\n","\n")
  
  cS.V = consensus.regex(secuencias.V, thr =  thr.V, verbose = verbose, plot.logo = F)
  if(verbose ==T)cat("\n")
  cS.J = consensus.regex(secuencias.J, thr = thr.J, verbose = verbose, plot.logo = F)
  if(verbose ==T)cat("\n")
  inicio = Sys.time()
  
  J.indexes.f = c()
  J.indexes.r = c()
  J.sequences.f = c()
  J.sequences.r = c()
  V.indexes.f = c()
  V.indexes.r = c()
  V.sequences.f = c()
  V.sequences.r = c()
  V.end.sequences.f = c()
  V.end.sequences.r = c()
  
  V.gg.f = genomic.gregexpr(genome.data,cS.V[1],selected = 1)
  V.gg.r = genomic.gregexpr(genome.data,cS.V[2],selected = 1)
  J.gg.f = genomic.gregexpr(genome.data,cS.J[1],selected = 1)
  J.gg.r = genomic.gregexpr(genome.data,cS.J[2],selected = 1)
  
  if(sum(V.gg.f$hits.chart)>0){
    V.indexes.f = as.numeric(V.gg.f$gregexpr.hits[[1]][[1]]+nchar(secuencias.V[1])-cS.V$starts.at+1)-V.length
    V.sequences.f = do.call(c,sapply(V.indexes.f,function(x)subseq(genome.data,start = x,width = V.length+extend.V)))
    V.end.sequences.f = do.call(c,sapply(V.indexes.f+V.length-3,function(x)subseq(genome.data,start = x,width =V.end.length)))
  }
  if(sum(V.gg.r$hits.chart)>0){
    V.indexes.r = as.numeric(V.gg.r$gregexpr.hits[[1]][[1]]-1-(nchar(secuencias.V[1])-(cS.V$starts.at+cS.V$pattern.length)))
    V.indexes.r = V.indexes.r[V.indexes.r-V.end.length+3>0] # This removes Vs that are in the edge
    V.sequences.r = reverseComplement(do.call(c,sapply(V.indexes.r,function(x)subseq(genome.data,start = x-extend.V, width = V.length+extend.V))))
    V.end.sequences.r = reverseComplement(do.call(c,sapply(V.indexes.r-V.end.length+3,function(x)subseq(genome.data,start = x, width = V.end.length))))
    
  }
  
  
  
  if(sum(J.gg.f$hits.chart)>0){
    J.indexes.f = as.numeric(J.gg.f$gregexpr.hits[[1]][[1]]+nchar(secuencias.J[1])-cS.J$starts.at+1)-J.length
    J.sequences.f = do.call(c,sapply(J.indexes.f,function(x)subseq(genome.data,start = x,width = J.length)))
  }
  if(sum(J.gg.r$hits.chart)>0){
    J.indexes.r = J.gg.r$gregexpr.hits[[1]][[1]]-1-(nchar(secuencias.J[1])-(cS.J$starts.at+cS.J$pattern.length))
    J.sequences.r = reverseComplement(do.call(c,sapply(J.indexes.r,function(x)subseq(genome.data,start = x,width = J.length))))
  }
  
  VJ.index.list = list(V.indexes.f=V.indexes.f,V.indexes.r=V.indexes.r,J.indexes.f=J.indexes.f,J.indexes.r=J.indexes.r)
  VJ.sequences.list = list(V.sequences.f = V.sequences.f,V.sequences.r = V.sequences.r,J.sequences.f = J.sequences.f,J.sequences.r = J.sequences.r,V.end.sequences.f = V.end.sequences.f,V.end.sequences.r = V.end.sequences.r)
  #VJ.index.list.C_GT = list(V.indexes.f+V.length-3,V.indexes.r+2,J.indexes.f+J.length,J.indexes.r)
  
  df.index = as.numeric(unlist(VJ.index.list))
  df.orientation = rep(c("F","R","F","R"),sapply(VJ.index.list,length))
  
  df.gen = rep(paste0(c("V","V","J","J"),gene),sapply(VJ.index.list,length))
  VJ.df = data.frame(cbind(df.gen,df.index,df.orientation))
  colnames(VJ.df) = c("segment","position","sense")
  rownames(VJ.df) = NULL
  if (length(VJ.df[[1]])==0)return (NULL)
  VJ.df.summary = VJ.df.summarizer(VJ.df,nchar(genome.data))
  
  if(verbose == T)print(difftime(Sys.time(),inicio))
  
  for (i in which(sapply(VJ.index.list,length)>0)){names(VJ.sequences.list[[i]]) = VJ.index.list[[i]]
  if(i %in% 1:2)names(VJ.sequences.list[[i+4]]) = names(VJ.sequences.list[[i]])}
  
  
  return(list(VJ.df.summary = VJ.df.summary[[1]],VJ.df = VJ.df.summary[[2]],VJ.index.list = VJ.index.list, VJ.sequences.list = VJ.sequences.list))
}

consensus.regex = function(sequences,thr,verbose = F,plot.logo = F){
  if (length(unique(nchar(sequences)))>1) stop("sequences must be of equal length")
  ec.seq =entropy.content(consensusMatrix(sequences))
  e.c.filter=ec.seq<=thr
  ss.cS = strsplit(consensusString(consensusMatrix(sequences)),"")[[1]]
  string.cS = rep(".",nchar(sequences)[1])
  string.cS[e.c.filter] =ss.cS[e.c.filter]
  cS.string = sub("\\.*$","",sub("\\.*","",paste(string.cS,collapse = "")))
  cS.string.rc = gsub("N",".",reverseComplement(DNAString(gsub("\\.","N",cS.string))))
  starts.at = as.numeric(which.max(table(factor(regexpr(cS.string,sequences),levels = 1:nchar(sequences[[1]])))))
  pattern.length = nchar(cS.string)
  if (verbose==T){cat("pattern entropy",2*sum(e.c.filter),"\n")
    cat("total entropy loss",sum(ec.seq[e.c.filter]),"\n")
    cat("starts at",starts.at,"--- pattern length",pattern.length,"\n")
    cat("consensus string",cS.string,"\n")}
  if (plot.logo==T)plot(ggseqlogo(consensusMatrix(subseq(sequences,starts.at,width = nchar(cS.string[1])))[1:4,]))
  return(list(fwd = cS.string, rev = cS.string.rc, both.rf = paste(cS.string,cS.string.rc,sep = "|"), pattern.entropy = 2*sum(e.c.filter),entropy.loss = sum(ec.seq[e.c.filter]),starts.at = starts.at,pattern.length = pattern.length))}

genomic.gregexpr = function(genome.data,pattern,selected,min.dist = 100000){
  gg.genome.data = lapply(selected,function(i)gregexpr(pattern,genome.data[[i]]))
  hits.chart = sapply(1:length(gg.genome.data),function(i)sort(as.numeric(table(cumsum(diff(gg.genome.data[[i]][[1]])>min.dist))),decreasing = T)[1:10])
  hits.chart[is.na(hits.chart)] = 0
  hits.order = order(colSums(hits.chart),decreasing = T)
  if(length(selected)==1)return(list(hits.chart = hits.chart,gregexpr.hits = gg.genome.data))
  df = data.frame(hits.chart[,hits.order][,1:min(10,ncol(hits.chart))])
  colnames(df) = selected[hits.order][1:min(10,ncol(hits.chart))]
  df = df[rowSums(df==0)!=ncol(df),]
  return(list(gregexpr.hits = gg.genome.data,hits.chart = df))
}

genomic.pattern.locator.mk2 = function(genomic.data,pattern,DNA.pre,DNA.length,verbose = F){
  genomic.data = paste(paste(rep("N",DNA.pre+DNA.length),collapse = ""),genomic.data,paste(rep("N",DNA.pre+DNA.length),collapse = ""),sep = "")
  all.6.rf = tr.DNA6(genomic.data)
  hits = lapply(1:6,function(j)sapply(1:length(pattern),function(i)gregexpr(pattern[i],all.6.rf[[j]])))
  if(sum(sapply(hits,function(x)sum(unlist(x)==-1)))==6){
    if(verbose == T)cat("Not found!")
    return(NA)}
  #if(sum(sapply(hits,function(i)(unlist(i)==-1)))==6)return("Not found!")
  indices = lapply(1:6,function(i)which(unlist(sapply(hits[[i]],function(x)x[1]!=-1))))
  indices.non0 = which(sapply(indices,sum)!=0)
  lista = list()
  largo=nchar(genomic.data)
  for (i in indices.non0){
    if (i <=3){
      lista[[i]]=list(pattern[indices[[i]]],sapply(hits[[i]][indices[[i]]],function(x)as.numeric(x*3+i-3)-DNA.pre),sapply(hits[[i]][indices[[i]]],function(x)sapply(x,function(y)as.character(subseq(genomic.data,y*3+i-3-DNA.pre,width = DNA.length+DNA.pre)))),"F")}else{
        lista[[i]]=list(pattern[indices[[i]]],sapply(hits[[i]][indices[[i]]],function(x)(largo+3)-as.numeric(x*3)-(i-4)-(DNA.length-1)),sapply(hits[[i]][indices[[i]]],function(x)sapply(x,function(y)as.character(subseq(genomic.data,(largo+3)-as.numeric(y*3)-(i-4)-(DNA.length-1),width = DNA.length+DNA.pre)))),"R")
      }
  }
  
  resultado = list()
  for (rf in indices.non0){
    df = data.frame()
    for (k in 1:length(lista[[rf]][[2]])){
      if(length(lista[[rf]][[1]])==1 & length(lista[[rf]][[2]])>1)lista[[rf]][[1]] = rep(lista[[rf]][[1]],length(lista[[rf]][[2]]))
      
      temp.df =data.frame(lista[[rf]][[1]][k],lista[[rf]][[2]][k],lista[[rf]][[3]][k],lista[[rf]][[4]])
      names(temp.df) = c("translated seq","index 5'","germline seq","orientation")
      df <- rbind(df,temp.df)
    }
    resultado[[rf]] = df
  }
  final = do.call(rbind,resultado)
  final2=final[!duplicated(final[,2]),]
  res=final2[order(final2[,2]),]
  
  secuencias = c()
  for (i in 1:nrow(res))if (res[i,][4] == "R")secuencias[i]=(as.character(reverseComplement(DNAStringSet(res[,3])[i])))else(secuencias[i]=as.character(DNAStringSet(res[,3])[i]))
  res$`index 5'` = res$`index 5'` - (DNA.pre+DNA.length)
  
  
  return(list(res,secuencias))
  
}

VJ.score = function(gpe.VJ,VJ.df,secuencias.V,secuencias.J,from.gpe = T){
  Vs.fr.indexes = NULL
  Js.fr.indexes = NULL
  
  if(from.gpe){gpe.VJ.sl = gpe.VJ$VJ.sequences.list} else gpe.VJ.sl = gpe.VJ
  gpe.filter = which(sapply(gpe.VJ.sl,length)>0)
  
  if(sum((1:2) %in% gpe.filter)>0){
    Vs.fr.indexes = (1:2)[1:2%in%gpe.filter]
    V.length = nchar(gpe.VJ.sl[[Vs.fr.indexes[1]]])[1]
    V.cut.length = min(nchar(secuencias.V[1]),V.length)
    cM.Vs= consensusMatrix(subseq(secuencias.V,-V.cut.length))[1:4,]
    for (fr in Vs.fr.indexes)
      
      
      PWM.V.scores = lapply(Vs.fr.indexes,function(fr)sapply(1:length(gpe.VJ.sl[[fr]]),function(seq)PWMscoreStartingAt(cM.Vs,subseq(gpe.VJ.sl[[fr]],-V.cut.length)[[seq]],1)))
    V.in.frame = lapply(Vs.fr.indexes,function(fr)ns.filter(translate(gpe.VJ.sl[[fr]],if.fuzzy.codon = "solve")))
  }else{PWM.V.scores = NA
  V.in.frame = NA
  Vend.pattern = NA
  }
  if(sum(3:4 %in% gpe.filter)>0){
    
    Js.fr.indexes = (3:4)[3:4%in%gpe.filter]
    J.cut.length = nchar(secuencias.J)[1] - (as.numeric(names(which.max(table(regexpr("TT.GG....GG",secuencias.J)))))-9) + 1
    
    cM.Js= consensusMatrix(subseq(secuencias.J,-J.cut.length))[1:4,]
    PWM.J.scores = lapply(Js.fr.indexes,function(fr)sapply(1:length(gpe.VJ.sl[[fr]]),function(seq)PWMscoreStartingAt(cM.Js,subseq(gpe.VJ.sl[[fr]],-J.cut.length)[[seq]],1)))
    J.in.frame = lapply(Js.fr.indexes,function(fr)ns.filter(translate(subseq(gpe.VJ.sl[[fr]],-J.cut.length),if.fuzzy.codon = "solve")))
    
    FGXG.pattern = lapply(Js.fr.indexes,function(fr)grepl("FG.G",translate(subseq(gpe.VJ.sl[[fr]],-J.cut.length+9,width = 12),if.fuzzy.codon = "solve")))
    WGXG.pattern = lapply(Js.fr.indexes,function(fr)grepl("WG.G",translate(subseq(gpe.VJ.sl[[fr]],-J.cut.length+9,width = 12),if.fuzzy.codon = "solve")))
  }else{PWM.J.scores = NA
  J.in.frame = NA
  FGXG.pattern = NA
  }
  VJ.df$in_frame = NA
  VJ.df$PWM.score = NA
  VJ.df$FGxG = NA
  VJ.df$WGxG = NA
  if(!is.null(Vs.fr.indexes)){for(i in Vs.fr.indexes)VJ.df$in_frame[grepl("V",VJ.df$segment)& VJ.df$sense==c("F","R")[i]] = V.in.frame[[which(Vs.fr.indexes==i)]]
  for(i in Vs.fr.indexes)VJ.df$PWM.score[grepl("V",VJ.df$segment)& VJ.df$sense==c("F","R")[i]] = PWM.V.scores[[which(Vs.fr.indexes==i)]]
  }
  if(!is.null(Js.fr.indexes)){for(i in Js.fr.indexes)VJ.df$in_frame[grepl("J",VJ.df$segment)& VJ.df$sense==c("F","R")[i-2]] = J.in.frame[[which(Js.fr.indexes==i)]]
  for(i in Js.fr.indexes)VJ.df$PWM.score[grepl("J",VJ.df$segment)& VJ.df$sense==c("F","R")[i-2]] = PWM.J.scores[[which(Js.fr.indexes==i)]]
  for(i in Js.fr.indexes)VJ.df$FGxG[grepl("J",VJ.df$segment)& VJ.df$sense==c("F","R")[i-2]] = FGXG.pattern[[which(Js.fr.indexes==i)]]
  for(i in Js.fr.indexes)VJ.df$WGxG[grepl("J",VJ.df$segment)& VJ.df$sense==c("F","R")[i-2]] = WGXG.pattern[[which(Js.fr.indexes==i)]]
  }
  VJ.df$in_frame[grepl("C",VJ.df$segment)] = T
  return(VJ.df)
}

RSS.scorer.iterative = function(cM.RSS,VJ.sequences,check.to = 160,iter.max = 1){  
  cM.RSS.0 = cM.RSS
  cM.diff = 1
  iter= 0
  while(iter<iter.max & cM.diff!=0){
    iter = iter + 1
    cat("iteration ",iter,"\n")
    if((dim(cM.RSS.0)[2]+check.to-1)>nchar(VJ.sequences)[[1]]){check.to = nchar(VJ.sequences[1])-dim(cM.RSS.0)[2]+1
    cat("check.to moved to fit J.sequences width. New position ->",check.to,"\n")
    }
    VJ.sequences = paste(VJ.sequences,paste(rep("N",dim(cM.RSS.0)[2]),collapse = ""),sep = "") #extend in case some sequence has a match at the end to avoid problems with cutting of the sequences for making the new consensus
    score.RSS = sapply(1:length(VJ.sequences),function(j)sapply(1:check.to,function(i)PWMscoreStartingAt(cM.RSS.0,VJ.sequences[[j]],i)))
    score.RSS.max = apply(score.RSS,2,which.max)  
    ss.RSS= subseq(VJ.sequences,score.RSS.max,width = dim(cM.RSS)[2])
    cM.RSS.1 =consensusMatrix(ss.RSS)[rownames(consensusMatrix(ss.RSS))%in%(c("A","C","G","T")),]
    cM.diff = sum(abs(cM.RSS.1-cM.RSS.0))
    cM.RSS.0 = cM.RSS.1
    cat("1-norm distance to previous iteration:",cM.diff,"\n")
  }
  
  return(list(score.RSS = score.RSS,RSS.index = score.RSS.max,ss.RSS = ss.RSS, cM.RSS =cM.RSS.1)) 
}
