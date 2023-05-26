# Wrapper functions used in downstream analysis

fastq.loader =function(archivos,carpetas,s){
  s.Fastq.list = list()
  for(i in 1:length(archivos))
    s.Fastq.list[[i]] = mapply(function(x,y)assign(x,readFastq(y)),archivos[[i]][(s*2-1):(s*2)],carpetas[[i]][(s*2-1):(s*2)])  #load the fastq files (needs libraries)
  return(s.Fastq.list)
}

spp.repseq.mc = function(s,spp,n.Cs=2,umi.thrs = 2,thrs.reads.list = 2,TCR.dicc = spp.TCR.dict,genomes.list = spp.genome.data.list, pattern.list = spp.patterns,V.mm.L = 180){f.l_n.s = lapply(1:1,function(i)fastq.loader(f.o_n$archivos,f.o_n$carpetas,s))
Rs = lapply(f.l_n.s,function(x)lapply(1:2,function(i)do.call(c,lapply(x,function(y)sread(y[[i]])))))[[1]]
return(repseq.pipeline.mk1.3(spp,1:n.Cs,Rs = Rs,umi.thr = umi.thrs, thr.reads.list = thrs.reads.list,TCR.dict = TCR.dicc,genomes.data.list = genomes.list,patterns = pattern.list,V.mm.l = V.mm.L))}

repseq.pipeline.mk1.3 = function(sp,sp.indexes,Rs,umi.thr = 2, thr.reads.list=2,umi.exclude = "TTTTTTTT|TTTT.TTTT", TCR.dict = spp.TCR.dict, genomes.data.list = spp.genome.data.list, patterns = spp.patterns, V.mm.l = 180){
  parsed = list()
  ec.list = list()
  #if(!sp %in%c('BS','GB','ST','MF','ZF','RT','MK','LF','MM','LA'))stop("Unknown species")
  if (length(umi.thr)==1 & length(sp.indexes)>1) umi.thr = rep(umi.thr,length(sp.indexes))
  if (length(thr.reads.list)==1 & length(sp.indexes)>1) thr.reads.list = rep(thr.reads.list,length(sp.indexes))
  {
    V.df.list = lapply(sp.indexes,function(i)TCR.dict[[sp]][[i]]$V.df)
    V.GR.dict.list = lapply(sp.indexes,function(i)TCR.dict[[sp]][[i]]$V.GR.dict)
    J.GR.dict.list = lapply(sp.indexes,function(i)TCR.dict[[sp]][[i]]$J.GR.dict)
    V.alleles.list = lapply(sp.indexes,function(i)TCR.dict[[sp]][[i]]$V.alleles)
    J.alleles.list = lapply(sp.indexes,function(i)TCR.dict[[sp]][[i]]$J.alleles)
    genome.data.list = genomes.data.list[[sp]]
    rte = repseq.table.extractor.mk2.1(Rs[[1]],Rs[[2]],C.rc.pattern.list = patterns[[sp]]$C.rc.patterns[sp.indexes],start.of.triplet.list = patterns[[sp]]$sot.list[sp.indexes],CDR3.consensus.list = patterns[[sp]]$CDR3.consensus[sp.indexes],J.consensus.list = patterns[[sp]]$J.consensus[sp.indexes],CDR3.consensus.length.list = patterns[[sp]]$CDR3.consensus.lengths[sp.indexes],V.dict.list = lapply(V.df.list,function(x)x$nt),umi_cdr3_threshold = umi.thr, UMI.threshold = min(umi.thr),V.mm.largo = V.mm.l )
  }
  if (is.null(rte)) return (NULL)
  for (i in 1:length(sp.indexes)){
    if(is.null(rte$clone.table.list[[i]])){ parsed[[i]] = NA
    ec.list[[i]] = NA
    next}
    rte.bind = do.call(rbind,rte$clone.table.list[[i]])
    rte.bind = rte.bind[!grepl(umi.exclude,rte.bind$UMI) & rte.bind$Freq>=thr.reads.list[[i]],]
    rte.bind$V.segment = V.translate(V.segment = rte.bind$V.segment,hit.list = V.df.list[[i]]$candidates.hit,V.GR.dict = V.GR.dict.list[[i]])
    rte.bind.ec = error.correction.mk6.clustered(rte.bind)
    ec.list[[i]] = rte.bind.ec
    V.ends = c(as.character(getSeq(genome.data.list[[i]],flank(V.GR.dict.list[[i]]-3,width = 100,start = F))),toupper(V.alleles.list[[i]]))
    J.nt = subseq(c(as.character(getSeq(genome.data.list[[i]],resize(J.GR.dict.list[[i]],200,fix = "end"))),toupper(J.alleles.list[[i]])),1,-4-as.numeric(grepl('a',names(rte$clone.table.list[i]))*3))
    parsed[[i]] = SP.TCR.data.parser.mk2(rte.bind.ec$corrected.list,subseq(V.ends,1,30),J.nt,21,100,1,60,19,V.dict = names(V.GR.dict.list[[i]]))
  }
  names(parsed) = names(TCR.dict[[sp]])[sp.indexes]
  names(ec.list) = names(TCR.dict[[sp]])[sp.indexes]
  return(list(parsed = parsed,rte = rte,ec = ec.list))
}

repseq.table.extractor.mk2.1 = function(R1,R2,C.rc.pattern.list,start.of.triplet.list,CDR3.consensus.list,J.consensus.list,CDR3.consensus.length.list,V.dict.list,UMI.threshold = 2,umi_cdr3_threshold =2,V.largos = c(180,170,150,140,120,80,60),V.mm.largo = 180){
  {
    inicio1 = Sys.time()
    inicio = Sys.time()
    nCs = length(C.rc.pattern.list)
    
    s.umi = UMI.extractor.3(R1,R2,umi.threshold = UMI.threshold)
    if (s.umi$stats$UMI.threshold ==0) return(NULL)
    print ("Umi.processing")
    print (difftime(Sys.time(),inicio))
    inicio = Sys.time()
    
    s.CDR3.e = lapply(1:nCs,function(i)CDR3.extractor(s.umi$R1,C.rc.pattern.list[[i]],start.of.triplet.list[[i]]))
    
    if (nCs>1){
      double.Cs=rowSums(do.call(cbind,lapply(1:nCs,function(i)s.CDR3.e[[i]]$filter.CDR3.extractor)))
      print(paste(sum(double.Cs>1),"double Cs"))
      if (sum(double.Cs>1)>0)for (i in 1:nCs){
        s.CDR3.e[[i]]$translated.reads =  s.CDR3.e[[i]]$translated.reads[!which(s.CDR3.e[[i]]$filter.CDR3.extractor) %in%which(double.Cs>1)]
        s.CDR3.e[[i]]$reading.frames =  s.CDR3.e[[i]]$reading.frames[!which(s.CDR3.e[[i]]$filter.CDR3.extractor) %in%which(double.Cs>1)]
        s.CDR3.e[[i]]$filter.CDR3.extractor[double.Cs>1] = FALSE
        
      }
    }
    
    print ("CDR3 extraction")
    print (difftime(Sys.time(),inicio))
    inicio = Sys.time()
    
    # Canonical V end trimming  
    s.CDR3.tr = lapply(1:nCs,function(i)Canonical.CDR3.j.trimmer.Mk2(s.CDR3.e[[i]]$translated.reads,reverseComplement(s.umi$R1[s.CDR3.e[[i]]$filter.CDR3.extractor]),s.CDR3.e[[i]]$reading.frames, CDR3.consensus.list[[i]],J.consensus.list[[i]],pre.C = CDR3.consensus.length.list[[i]]*3-5))
    
    CDR3.tr_cdr3s = lapply(1:nCs,function(i)do.call(c,s.CDR3.tr[[i]]$CDR3s))
    CDR3.tr_v.ends = lapply(1:nCs,function(i)do.call(c,s.CDR3.tr[[i]]$V.end))
    CDR3.tr_fil = lapply(1:nCs,function(i)do.call(c,s.CDR3.tr[[i]]$CDR3.filter))
    UMI_CDR3 = lapply(1:nCs,function(i)paste(s.umi$UMI[s.CDR3.e[[i]]$filter.CDR3.extractor][CDR3.tr_fil[[i]]],CDR3.tr_cdr3s[[i]]))
    filter_UMI_CDR3 = lapply(1:nCs,function(i)CDR3.tr_fil[[i]][as.logical(table(UMI_CDR3[[i]])[UMI_CDR3[[i]]]  >= umi_cdr3_threshold[[i]])])
    
    print ("Canonical trimming")
    print (difftime(Sys.time(),inicio))
    inicio = Sys.time()
    
    secuencias.r2 = lapply(1:nCs,function(i)s.umi$R2[s.CDR3.e[[i]]$filter.CDR3.extractor][filter_UMI_CDR3[[i]]])
    secuencias.r1 = lapply(1:nCs,function(i)CDR3.tr_v.ends[[i]][CDR3.tr_fil[[i]]%in%filter_UMI_CDR3[[i]]])
    
  }
  s.Vs.list = list()
  for(i in 1:nCs){
    s.Vs.list[[i]] = R1.V.call(secuencias.r1[[i]],V.dict.list[[i]],rotular = F,largos = V.largos)
    filter.solve = s.Vs.list[[i]]$v_n!=1 | s.Vs.list[[i]]$l <=100
    vcPD.Vs = vcountPDict(PDict(DNAStringSet(subseq(V.dict.list[[i]],-min(V.mm.largo,nchar(V.dict.list[[i]])))),tb.start = min(V.mm.largo,nchar(V.dict.list[[i]]))),(secuencias.r1[[i]][filter.solve]),max.mismatch = 1)
    filter.vcPD = colSums(vcPD.Vs) == 1
    vcPD.Vs.w = which(vcPD.Vs[,filter.vcPD]>0,arr.ind = T)
    if(length(vcPD.Vs.w)==1) vcPD.Vs.w = as.matrix(vcPD.Vs.w) # new line to correct dimensions problem when length of vcPD.Vs.w is 1
    s.Vs.list[[i]]$v[filter.solve][filter.vcPD] = vcPD.Vs.w[,1]
    s.Vs.list[[i]]$l[filter.solve][filter.vcPD] = paste0(as.character(V.mm.largo),".1")
    s.Vs.list[[i]]$v_n[filter.solve][filter.vcPD] = 1
  }
  
  #vcPD.Vs.paste=tapply(vcPD.Vs.w[,1],vcPD.Vs.w[,2],paste,collapse = ',')
  
  s.Vs = lapply(s.Vs.list,function(x)x$v)
  s.Vs.n = lapply(s.Vs.list,function(x)table(x$v_n))
  s.Vs.l = lapply(s.Vs.list,function(x)table(x$l))
  
  
  
  print(" ")
  print ("V calling")
  print (difftime(Sys.time(),inicio))
  inicio = Sys.time()
  s.ftp.list = lapply(1:nCs,function(i)if(sum(s.Vs[[i]]!=0)>100)fast.table.parser.v1(s.umi,s.CDR3.e[[i]],s.CDR3.tr[[i]],filter_UMI_CDR3[[i]],s.Vs[[i]],threshold = 2))
  s.stats.list = list(s.umi$stats,C.filter = sapply(s.CDR3.e,function(x)sum(x[[2]])),Canonical.filter = sapply(CDR3.tr_fil,length),V.n = s.Vs.n,V.l =s.Vs.l)
  cc.90 = lapply(1:nCs,function(i)char.cumsum(secuencias.r1[[i]],0.90))
  excluded.Vs = lapply(1:nCs,function(i)sort(table(subseq(secuencias.r1[[i]][s.Vs[[i]]==0][nchar(secuencias.r1[[i]][s.Vs[[i]]==0])>=cc.90[[i]]],-cc.90[[i]])),decreasing = T)[1:100])
  names(excluded.Vs)= names(C.rc.pattern.list)
  names(s.ftp.list) = names(C.rc.pattern.list)
  
  print ("Table parsing")
  print (difftime(Sys.time(),inicio))
  inicio = Sys.time()
  print ("Total")
  print (difftime(Sys.time(),inicio1))
  return(list(clone.table.list = s.ftp.list, stats = s.stats.list,excluded.Vs = excluded.Vs))
  
}

UMI.extractor.3 = function(sequenceR1,sequenceR2,umi.threshold = 1,mm=0,short_UPM="GCAGAGT\\w{4}T\\w{4}T\\w{4}T",UMI.start= 6,UMI.length=16,R2.trim = 80, umi.exclude = "TTTTTTTT|TTTT.TTTT|N"){
  inicio = Sys.time()
  largo = length(sequenceR2)
  rr.UMI.start = gregexpr(short_UPM,subseq(sequenceR2,1,R2.trim))
  count_UPM = elementNROWS(rr.UMI.start)
  count_UPM[count_UPM==1][rr.UMI.start[count_UPM==1]==-1] = 0  
  find_UPM = unlist(rr.UMI.start[count_UPM==1])
  umis = as.character(subseq(sequenceR2[count_UPM==1],find_UPM+UMI.start,width = UMI.length))
  filto.umi.exclude = grep(umi.exclude,umis)
  filtro = (1:largo)[count_UPM==1][subseq(umis,1,1)=="T"&subseq(umis,6,6)=="T"&subseq(umis,11,11)=="T"&subseq(umis,16,16)=="T"]
  filtro = filtro[-filto.umi.exclude]
  umis = umis[-filto.umi.exclude]
  umis = subseq(umis,2,-2)
  
  umi.filter = as.numeric(table(umis)[umis])>=umi.threshold
  filtro = filtro[umi.filter]
  return = list(R1 = sequenceR1[filtro],
                R2 = sequenceR2[filtro],
                UMI = umis[umi.filter],
                filter = filtro,
                stats = list(Total = largo,UPM_count=table(count_UPM),UMI_structure=length(umi.filter), UMI.threshold = sum (umi.filter),Time = round(difftime(Sys.time(),inicio),2)))
}

CDR3.extractor <- function(step1.1,C.rc.pattern = "GGAGCT[AG]TGAC",start.of.triplet = 10){
  cat("\n",sep="")
  cat("CDR3.extractor",sep="")
  inicio=Sys.time()
  pattern_match = regexpr(C.rc.pattern,step1.1)
  rfs.1 = as.numeric(nchar(step1.1)-pattern_match)+start.of.triplet+1
  rfs.2 = rfs.1%%3 + 1
  resultados = c()
  
  filter.CDR3.extractor = pattern_match != -1
  ancho = (floor(nchar(step1.1[1])/3)-1)*3 # I'm assuming all sequences of equal length here
  
  RFj = suppressWarnings(translate(subseq(reverseComplement(step1.1[filter.CDR3.extractor]),rfs.2[filter.CDR3.extractor],width = ancho),if.fuzzy.codon = "solve" , no.init.codon = TRUE))
  cat ("  --  Sequences with C pattern: ",sum(filter.CDR3.extractor),"  (",round(sum(filter.CDR3.extractor)*100/length(pattern_match),1),"%)\n" , sep = "")
  print(round(difftime(Sys.time(),inicio,units = "auto"),digits = 2))
  return (list(translated.reads=RFj,filter.CDR3.extractor=filter.CDR3.extractor,reading.frames=rfs.2[filter.CDR3.extractor]))
}

Canonical.CDR3.j.trimmer.Mk2 = function(filtered.translated.reads, filtered.r1.reads , reading.frames.table, CDR3.consensus, J.consensus, post.J = 27,pre.C = 13){
  cat("\nCDR3 trimmer", sep ="")
  rf.orders = as.numeric(names(sort(table(reading.frames.table),decreasing = T)))
  J.match = regexpr(J.consensus,filtered.translated.reads)
  filtered.r1.reads.aa =tr.DNA3(filtered.r1.reads[J.match!=-1])
  CDR3.pos.1 = regexpr(CDR3.consensus,filtered.r1.reads.aa[[rf.orders[1]]])
  CDR3.pos.1 [CDR3.pos.1 >  J.match[J.match!=-1]] = -1
  CDR3.pos.2 = regexpr(CDR3.consensus,filtered.r1.reads.aa[[rf.orders[2]]][CDR3.pos.1==-1])
  CDR3.pos.2 [CDR3.pos.2 >  J.match[J.match!=-1][CDR3.pos.1==-1]] = -1
  CDR3.pos.3 = regexpr(CDR3.consensus,filtered.r1.reads.aa[[rf.orders[3]]][CDR3.pos.1==-1][CDR3.pos.2==-1])
  CDR3.pos.3 [CDR3.pos.3 >  J.match[J.match!=-1][CDR3.pos.1==-1][CDR3.pos.2==-1]] = -1
  CDR3.position.list = list(CDR3.pos.1,CDR3.pos.2,CDR3.pos.3)
  
  filter.list = list(which(J.match!=-1),which(J.match!=-1)[CDR3.pos.1==-1],which(J.match!=-1)[CDR3.pos.1==-1][CDR3.pos.2==-1])
  
  rf1= lapply(1:3,function(i)filter.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders-1)%%3+1)[i]])
  rf2= lapply(1:3,function(i)filter.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders+0)%%3+1)[i]])
  rf3= lapply(1:3,function(i)filter.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders+1)%%3+1)[i]])
  {
    #each rf gives the position of CDR3 consensus in seq (as aminoacid) whithout the need to add filters
    #for example if we translate all reading frames for tr.DNA3(filtered.r1.reads[rf1[[1]]])
    #we will see the sequences where the consensus is in the reading frame 1 and also the J and C are in reading frame 1
    #in all the other lists only the CDR3 consensus is in the right frame (in rf2[[3]] for example, 
    #CDR3 is in the right reading frame in the third position of the list, while C and J are not,
    #but are 1 nucleotide away from each other in all 3 elements of the list because of this command(rf.orders%%3+1))
  }
  cdr1 = lapply(1:3,function(i)CDR3.position.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders-1)%%3+1)[i]])
  cdr2 = lapply(1:3,function(i)CDR3.position.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders+0)%%3+1)[i]])
  cdr3 = lapply(1:3,function(i)CDR3.position.list[[i]][ CDR3.position.list[[i]]  !=-1 & reading.frames.table[filter.list[[i]]] == ((rf.orders+1)%%3+1)[i]])
  
  CDR3s.1 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf1[[i]]],rf.orders[i]),cdr1[[i]]*3+pre.C,(J.match[rf1[[i]]]*3+post.J)))
  CDR3s.2 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf2[[i]]],rf.orders[i]),cdr2[[i]]*3+pre.C,(J.match[rf2[[i]]]*3+((rf.orders+0)%%3+1-rf.orders)[i]+post.J)))
  CDR3s.3 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf3[[i]]],rf.orders[i]),cdr3[[i]]*3+pre.C,(J.match[rf3[[i]]]*3+((rf.orders+1)%%3+1-rf.orders)[i]+post.J)))
  
  
  pre.CDR3s.1 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf1[[i]]],rf.orders[i]),1,cdr1[[i]]*3+pre.C+2))
  pre.CDR3s.2 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf2[[i]]],rf.orders[i]),1,cdr2[[i]]*3+pre.C+2))
  pre.CDR3s.3 = lapply(1:3,function(i)subseq(subseq(filtered.r1.reads[rf3[[i]]],rf.orders[i]),1,cdr3[[i]]*3+pre.C+2))
  
  canonical.filter.list = list(sort(unlist(rf1)),sort(unlist(rf2)),sort(unlist(rf3)))
  CDR3s.list = list(do.call("c",CDR3s.1)[order(unlist(rf1))],do.call("c",CDR3s.2)[order(unlist(rf2))],do.call("c",CDR3s.3)[order(unlist(rf3))])
  pre.CDR3s.list = list(do.call("c",pre.CDR3s.1)[order(unlist(rf1))],do.call("c",pre.CDR3s.2)[order(unlist(rf2))],do.call("c",pre.CDR3s.3)[order(unlist(rf3))])
  largo = length(filtered.translated.reads)
  kept = length(unlist(canonical.filter.list))
  cat ("  --  Found canonical CDR3s ",kept," (",round(100*kept/largo,1), "%)  --  In frame ", round(100*length(canonical.filter.list[[1]])/kept,1), "%\n" ,sep = "")
  return(list(CDR3s=CDR3s.list,CDR3.filter=canonical.filter.list,V.end = pre.CDR3s.list))
}

R1.V.call = function(CDR3.ends,V.dict,largos = c(),nombres = c(),largos.min = 24,rotular = F){
  if(rotular & length(nombres) == 0)nombres = names(V.dict)
  CDR3.n = nchar(CDR3.ends)
  if(length(largos)==0)largos = sort(unique(CDR3.n),decreasing = T)
  largos = largos[largos<=min(nchar(V.dict)) & largos >= largos.min]
  pb <- txtProgressBar(min = 1, length(largos), style = 3)
  largos.conflict = lapply(largos,function(x)as.numeric(table(as.character(subseq(V.dict,-x)))[as.character(subseq(V.dict,-x))]))
  largos.conflict.g = lapply(largos,function(x)graph.split(stringDist(as.character(subseq(V.dict,-x))),distancia = 0))
  largos.conflict.1 = lapply(largos.conflict.g,function(x)sapply(x,function(y)y[[1]]))
  largos.conflict.n = lapply(largos.conflict.g,function(x)sapply(x,length))
  if(rotular)largos.conflict.p = lapply(largos.conflict.g,function(x)sapply(x,function(y)paste(nombres[y],collapse = ','))) else largos.conflict.p = lapply(largos.conflict.g,function(x)sapply(x,function(y)paste(y,collapse = ',')))
  CDR3.v = rep(0,length(CDR3.n))
  CDR3.vn = rep(0,length(CDR3.n))
  CDR3.l = rep(0,length(CDR3.n))
  
  for(i in 1:length(largos)){
    
    
    setTxtProgressBar(pb, i)  
    l = largos[i]
    if(sum(CDR3.n>=l)==0)next
    mm=match(subseq(CDR3.ends[CDR3.vn!=1 & CDR3.n>=l],-l),subseq(V.dict,-l),nomatch = 0)
    mm.n = CDR3.vn[CDR3.vn!=1 & CDR3.n>=l]
    if(sum(mm!=0)==0)next
    mm.1 = match(mm[mm!=0],largos.conflict.1[[i]])
    #
    mm.new = rep(0,sum(mm!=0))
    mm.seq = c()
    mm.new[is.na(mm.1)] = 1
    if(sum(!is.na(mm.1))>0) mm.new[!is.na(mm.1)] = largos.conflict.n[[i]][mm.1[!is.na(mm.1)]]
    mm.old = CDR3.vn[CDR3.vn!=1 & CDR3.n>=l][mm!=0]
    
    
    mm.seq[is.na(mm.1)] = mm[mm!=0][is.na(mm.1)]
    if(sum(!is.na(mm.1))>0)mm.seq[!is.na(mm.1)] = largos.conflict.p[[i]][mm.1[!is.na(mm.1)]]
    
    
    CDR3.l[CDR3.vn!=1 & CDR3.n>=l][mm!=0][mm.new<mm.old | mm.old ==0] = largos[i]
    CDR3.v[CDR3.vn!=1 & CDR3.n>=l][mm!=0][mm.new<mm.old | mm.old ==0] = mm.seq[mm.new<mm.old | mm.old ==0]
    CDR3.vn[CDR3.vn!=1 & CDR3.n>=l][mm!=0][mm.new<mm.old | mm.old ==0] = mm.new[mm.new<mm.old | mm.old ==0]
  }
  if(rotular)CDR3.v[CDR3.vn==1] = nombres[as.numeric(CDR3.v[CDR3.vn==1])]
  return(list(v = CDR3.v,v_n = CDR3.vn,l= CDR3.l))}

fast.table.parser.v1 = function(step1.x,CDR3e.x,Ctr.x,filter_UMI_CDR3_x,Vxs,threshold = 2){
  # threshold is the minimum count at which UMI-CDR3s-Vs combinations are kept (should be selected based on read depth and sample size)
  #
  # threshold.single is the threshold below which single UMIs or CDR3s are eliminated to make the function run faster
  
  # Make an index for sample length
  
  # Make data frames for each sample
  frames.x.UMIs = as.character(step1.x$UMI[CDR3e.x$filter.CDR3.extractor][filter_UMI_CDR3_x])
  frames.x.CDR3s = as.character(do.call("c",Ctr.x$CDR3s)[do.call("c",Ctr.x$CDR3.filter)  %in% filter_UMI_CDR3_x])
  frames.x.Vs = as.character(Vxs)
  # Apply and combine threshold.single filter
  
  # Combine and remove columns where the count is less than 2 for either
  frames.umi.cdr3.x.df = split(data.frame(UMIs = frames.x.UMIs, CDR3s = frames.x.CDR3s), frames.x.Vs)
  
  # Count UMI-CDR3-Vs combinations and keep only above threshold (which should be selected based on read depth)
  frames.umi.count.x = lapply(frames.umi.cdr3.x.df[-1],function(x)separate(data.frame(table(paste(x[,1],x[,2],sep = ","))),1,into = c("UMI","CDR3")))
  frames.umi.count.x.2 = frames.umi.count.x[sapply(frames.umi.count.x,nrow)>0]
  frames.umi.count.x.3 = lapply(frames.umi.count.x.2,function(x)x[x[,3]>=threshold,])
  frames.umi.count.x.4 = frames.umi.count.x.3[sapply(frames.umi.count.x.3,nrow)>0]
  frames.umi.count.x.5 = lapply(1:length(frames.umi.count.x.4),function(i)cbind(frames.umi.count.x.4[[i]],V.segment=names(frames.umi.count.x.4)[i]))
  cat("Duplicated UMI+CDR3s per sample: " , sum(table(unlist(lapply(frames.umi.count.x.5,function(x)paste(paste(x[,1],x[,2],sep = ","),x[,4],sep = ","))))!=1))
  
  return (frames.umi.count.x.5)
}

V.translate = function(V.segment,hit.list,V.GR.dict){
  V.indexes = lapply(strsplit(unique(V.segment),','),function(x)as.numeric(x))
  V.df.keys = hit.list
  
  if(sum(grepl(',',hit.list))>0)V.df.keys = lapply(strsplit(hit.list,','),function(x)as.numeric(unique(x)))
  V.df.keys[!is.na(hit.list)] = lapply(V.df.keys[!is.na(hit.list)],function(x)paste(sort(unique(names(V.GR.dict)[x])),collapse = ','))
  if(sum(is.na(hit.list))>0)V.df.keys[is.na(hit.list)] = paste("nf",formatC(1:sum(is.na(hit.list)),width = floor(log10(sum(is.na(V.df.keys)))+1),flag = 0),sep = '.')
  V.index.list = lapply(V.indexes,function(i)paste(unique(V.df.keys[i]),collapse = ','))
  
  return(do.call(c,V.index.list[match(V.segment,unique(V.segment))]))
}

error.correction.mk6.clustered = function(secuencias){
  # adj is an adjacency matrix 
  adj = lapply(1:length(unique(secuencias[,4])),function(i)which(strsplit(unique(secuencias[,4]),",")%in%strsplit(unique(secuencias[,4]),",")[[i]]))
  UMI.CDR3.paste = lapply(split(secuencias,secuencias[,4]),function(x)apply(x[,1:2],1,paste0,collapse = "."))
  UMI_CDR3.matrix = sapply(1:length(UMI.CDR3.paste),function(i)sapply(1:length(UMI.CDR3.paste),function(j)sum(UMI.CDR3.paste[[i]]%in%UMI.CDR3.paste[[j]])))>0
  UMI_CDR3.gb = graph.builder(which(matrix.builder(UMI_CDR3.matrix),arr.ind = T))
  adj2 = adj
  if (length(UMI_CDR3.gb)>0){
    
    for (i in 1:length(UMI_CDR3.gb)) adj2[[as.numeric(UMI_CDR3.gb[[i]][[1]])]] = unique(c(adj2[[as.numeric(UMI_CDR3.gb[[i]][[1]])]] ,as.numeric(UMI_CDR3.gb[[i]][[-1]]) ))
  }
  g2 <- graph_from_adj_list(adj2)
  
  V.clusters = lapply(1:clusters(g2)$no,function(i)unique(secuencias[,4])[clusters(g2)$membership==i])
  by.cluster = lapply(1:clusters(g2)$no,function(i)secuencias[secuencias[,4]%in%V.clusters[[i]],])
  secuencias.ec = lapply(by.cluster,error.correction.mk5)
  
  corrected.same.length = do.call(rbind,lapply(secuencias.ec,function(x)x[[1]]))
  corrected.same.length.sp=split(corrected.same.length,corrected.same.length$UMI)
  corrected.same.length.g = lapply(corrected.same.length.sp[sapply(corrected.same.length.sp,nrow)>1],function(x)unlist(graph.builder(which(matrix.builder(stringDist(x$CDR3))<=3,arr.ind = T))))
  to.exclude = unlist(lapply(corrected.same.length.sp[sapply(corrected.same.length.sp,nrow)>1][sapply(corrected.same.length.g,length)>0],function(x)rownames(x)[-which.max(x$Freq)]))
  if(!is.null(to.exclude)){corrected.clone.list = corrected.same.length[-match(to.exclude,rownames(corrected.same.length)),]
  errors.2 = corrected.same.length.sp[sapply(corrected.same.length.sp,nrow)>1][sapply(corrected.same.length.g,length)>0]
  errors.3 = lapply(1:length(errors.2),function(i)errors.2[[i]][corrected.same.length.g[sapply(corrected.same.length.g,length)>0][[i]],])}else {
    corrected.clone.list = corrected.same.length
    errors.3 = NA}
  
  return(lista = list(corrected.list= corrected.clone.list, errors = do.call(c,lapply(secuencias.ec,function(x)x[[2]])), errors.length = errors.3))
}

error.correction.mk5 =  function(secuencias) {
  retain = c()
  discard= list()
  if (nrow(secuencias)==1)retain = rbind(retain,secuencias) else{
    nf = nodes.finder.UMI(secuencias[,1])
    if(sum(nf)==0) retain = rbind(retain,secuencias) else{
      target.umi = DNAStringSet(secuencias[,1])[nf] # this assumes that there won't be many errors unconnected, but it might not be a safe assumption for CDR3s (it probably is for UMIs)
      vD.umi=vcountPDict(PDict(target.umi,tb.start = 5,tb.end = 5),target.umi,max.mismatch = 1)
      target.cdr3 = DNAStringSet(secuencias[,2])[nf]
      vD.cdr3 =vcountPDict(PDict(target.cdr3 ,tb.start = 1,tb.end = 1),target.cdr3 ,max.mismatch = 2)
      vD.joint = vD.umi+vD.cdr3
      res = which(vD.joint==2,arr.ind = TRUE)
      edges = res[res[,1]!=res[,2],]
      if (length(edges)==0) retain = rbind(retain,secuencias) else{
        gr = graph(as.vector(t(edges)))
        cls=clusters(gr)
        cls.list = lapply(seq_along(cls$csize)[cls$csize > 1], function(x) V(gr)[cls$membership %in% x])
        pre.sel = lapply(cls.list,function(clu)secuencias[(nf),][clu,])
        pre.sel = lapply(pre.sel,function(x)x[order(x[,3],decreasing = T),])
        sel = do.call(rbind,lapply(pre.sel,function(x)x[which.max(x$Freq),]))
        discard = append(discard,pre.sel)
        retain = rbind(retain,sel)
        retain = rbind(retain,secuencias[!(nf),])
        retain = rbind(retain,secuencias[nf,][!(1:nrow(secuencias[nf,])%in%(unlist(cls.list))),])
      }
    }
  }
  return (list(corrected.list = retain,errors = discard))
}

SP.TCR.data.parser.mk2 = function(corrected.clone.list,V.end.nt,J.nt,J.l.min,J.l.max,J.clone.subseq,J.corr.length,J.start.corr,V.dict){
  df.final = corrected.clone.list
  groupColumns = c("CDR3","V.segment")
  dataColumns = c("Freq")
  df.final = ddply(df.final, groupColumns, list(UMI.count = nrow, Total = function(x) colSums(x[dataColumns])))
  #Make Vs end list
  zf_V.sp = split(df.final,df.final$V.segment) #this split is automatically ordered by CDR3 (since it is the first column)
  V.end = SP.V.end_calling.mk2(df.final,SP.V.end.nt = V.end.nt,V.dict.names = V.dict)
  df.final = do.call(rbind,lapply(1:length(V.end),function(i)cbind(zf_V.sp[[i]],V.end=V.end[[i]])))
  # Make J list
  Js.df = SP.J_call(clone.list = as.character(df.final$CDR3),SP.J.Ns = J.nt,l.min = J.l.min,l.max = J.l.max,clone.subseq = J.clone.subseq,corr.length = J.corr.length)
  df.J_names = unlist(lapply(strsplit(as.character(Js.df$J.matches.final),","),function(x)paste0(names(J.nt)[as.numeric(x)],collapse = ",")))
  
  df.final$J.segment = df.J_names
  df.final$J.start = nchar(as.character(df.final$CDR3))-Js.df$J.start-J.start.corr
  
  #add a few more columns
  df.final$CDR3.C.to.F = as.character(subseq(DNAStringSet(df.final$CDR3,1,-28)))
  df.final$CDR3.aa= as.character(translate(subseq(DNAStringSet(df.final$CDR3,1,-28)),if.fuzzy.codon = "solve"))
  df.umi.sum = sum(df.final$UMI.count)
  df.Freq.sum = sum(df.final$Total.Freq)
  df.for.tcR = data.frame("Umi.count" =df.final$UMI.count , "Umi.proportion" =df.final$UMI.count/df.umi.sum , "Read.count" = df.final$Total.Freq , "Read.proportion" =  df.final$Total.Freq/df.Freq.sum , "CDR3.nucleotide.sequence" = df.final$CDR3.C.to.F , "CDR3.amino.acid.sequence" = df.final$CDR3.aa , "V.gene" = as.character(df.final$V.segment) , "J.gene" = as.character(df.final$J.segment), "D.gene" = rep("", length(df.final$UMI)) , "V.end" = df.final$V.end , "J.start" = df.final$J.start, "D5.end" = rep(-1, length(df.final$UMI)), "D3.end" = rep(-1, length(df.final$UMI)) , "VD.insertions" = rep(-1, length(df.final$UMI)) , "DJ.insertions" = rep(-1, length(df.final$UMI)) , "Total.insertions" = as.numeric(unlist(df.final$J.start -df.final$V.end -1)) ,stringsAsFactors = FALSE)  
  return(df.for.tcR)
}

SP.V.end_calling.mk2 = function(df.final,SP.V.end.nt,V.dict.names){
  cdr3_V.sp = split(df.final$CDR3,df.final$V.segment)
  
  #now as general for loop
  V.start.i = list()
  inicio = Sys.time()
  nombres = strsplit(names(cdr3_V.sp),",") #names of the Vs found from each split (sometimes more than one)
  
  for(i in 1:length(cdr3_V.sp)){
    if(!names(cdr3_V.sp)[[i]]%in%V.dict.names){ V.start.i[[i]] = rep(NA,length(cdr3_V.sp[[i]])) 
    next}
    V.start.i[[i]] = rep(0,length(cdr3_V.sp[[i]]))
    patron.1 = lapply(nombres[[i]], function(x)as.character(unique(SP.V.end.nt[names(SP.V.end.nt) %in% x])))
    patron = unique(unlist(patron.1))
    
    # here I check length of patron and add Ns to the short one
    maximo = max(nchar(patron))
    for (p in 1:length(patron)){
      if (nchar(patron[p])< max(nchar(patron)))patron[p] = paste0(patron[p],paste0(rep("N",maximo-nchar(patron[p])),collapse = ""),collapse = "")}
    
    Vs.end.matrix.list = lapply(1:maximo,function(j)match(subseq(as.character(cdr3_V.sp[[i]]),1,j),subseq(patron,1,j),nomatch = 0))
    for (k in length(Vs.end.matrix.list):1){
      V.start.i[[i]][V.start.i[[i]] == 0][Vs.end.matrix.list[[k]][V.start.i[[i]] == 0] %in% 1:length(patron)] = k
    }
  }
  return(V.start.i)
}

SP.J_call = function(clone.list, SP.J.Ns ,l.min = 21,l.max = 60,clone.subseq = 1, corr.length = 60){
  short.CDR3s = clone.list[nchar(clone.list)<l.max]
  correctivo = nchar(short.CDR3s)
  correctivo2 = nchar(clone.list)<l.max
  #I prolongue the shorter CDR3s to be able to perform a match search, with a DIFFERENT character than I extended Js "X" vs "N"
  X.CDR3s = unlist(lapply(short.CDR3s,function(x)paste0(paste0(c(rep("X",(l.max-nchar(x)))),collapse = ""),x,collapse = ""))) 
  #And replace
  clone.list[nchar(clone.list)<l.max] = X.CDR3s
  
  #J.matches = lapply(l.min:l.max,function(i)match(subseq(clone.list,-i,-clone.subseq),subseq(SP.J.Ns,(l.max+1-i)),nomatch = 0)) # replaced by line below, 25.03.21
  J.matches = lapply(l.min:l.max,function(i)match(subseq(clone.list,-i,-clone.subseq),subseq(SP.J.Ns,end = -1,width = i),nomatch = 0))
  
  J.start = rep(0,length(clone.list))
  J.matches.final = rep(0,length(clone.list))
  for (k in length(J.matches):1){
    J.start[J.start == 0][J.matches[[k]][J.start == 0] != 0 ] = k
    J.matches.final[J.matches.final == 0] = J.matches[[k]][J.matches.final==0]
  }
  #J.start[correctivo2] = J.start[correctivo2]-(60-correctivo)
  
  #correction for ambiguous calling
  sD.Ja = lapply(21:corr.length,function(i)stringDist(subseq(SP.J.Ns,((corr.length+1)-i)),method = "hamming"))
  conflicting.clusters = lapply(1:(corr.length-20),function(i)graph.builder(which(matrix.builder(sD.Ja[[i]])==0,arr.ind = TRUE)))
  df.J = data.frame(J.matches.final,J.start)
  for (i in 1:length(conflicting.clusters)){
    if (length(conflicting.clusters[[i]])>0)
      for (j in 1:length(conflicting.clusters[[i]])){
        df.J$J.matches.final[df.J$J.start == i][df.J[df.J$J.start == i,]$J.matches.final %in% as.numeric(conflicting.clusters[[i]][[j]])] = paste0(as.numeric(conflicting.clusters[[i]][[j]]),collapse = ",")
      }
  }
  
  return(df.J)
}

df.list.aggregate.UMIs= function(parsed.list,gene = "a",species = "",parsed.list.names = '',df.ncol = 17,df.col.exclude = c(-2,-4),na.rem = F,TCR.dict = spp.TCR.dict,data.list = spp.genome.data.list,correct.V.end = F){
  if(identical(parsed.list.names,''))parsed.list.names = names(parsed.list)
  pl.filter = sapply(1:length(parsed.list),function(i)length(parsed.list[[i]]$parsed[[gene]])>0)
  parsed.list = parsed.list[pl.filter]
  parsed.list.names = parsed.list.names[pl.filter]
  parsed.x= lapply(1:length(parsed.list),function(i)cbind(parsed.list[[i]]$parsed[[gene]], Sample.bio.name =  parsed.list.names[[i]]))
  parsed.x.rbind = do.call(rbind,parsed.x[sapply(parsed.x,length)==df.ncol])
  if(!na.rem)parsed.x.rbind[is.na(parsed.x.rbind)] = -1
  parsed.x.agg = aggregate(cbind(Umi.count , Read.count) ~ ., data = parsed.x.rbind[,c(df.col.exclude)], FUN = sum, na.rm = TRUE)
  df.ncol.2 = df.ncol-length(df.col.exclude)
  parsed.x.agg = parsed.x.agg[,c((df.ncol.2-1):df.ncol.2,1:(df.ncol.2-2))]
  parsed.x.agg$Umi.proportion = unlist(with(parsed.x.agg,tapply(Umi.count,Sample.bio.name,prop.table)))
  parsed.x.agg$Read.proportion = unlist(with(parsed.x.agg,tapply(Read.count,Sample.bio.name,prop.table)))
  parsed.x.agg$pc.nt = pc.generator.nt(parsed.x.agg)
  pc.nt.tb = table(unlist(tapply(parsed.x.agg$pc.nt,parsed.x.agg$Sample.bio.name,unique)))
  parsed.x.agg$pub = pc.nt.tb[parsed.x.agg$pc.nt]
  parsed.x.agg$RSS.V = TCR.dict[[species]][[gene]]$V.GR.dict$heptamer.start[match(parsed.x.agg$V.gene, names(TCR.dict[[species]][[gene]]$V.GR.dict))]
  parsed.x.agg$RSS.J = TCR.dict[[species]][[gene]]$J.GR.dict$heptamer.end[match(parsed.x.agg$J.gene, names(TCR.dict[[species]][[gene]]$J.GR.dict))]
  parsed.x.agg$L = nchar(parsed.x.agg$CDR3.nucleotide.sequence)
  parsed.x.agg$RSS.V[is.na(parsed.x.agg$RSS.V)] = NA
  parsed.x.agg$RSS.J[is.na(parsed.x.agg$RSS.J)] = NA
  parsed.x.agg$Total.insertions[is.na(parsed.x.agg$RSS.V)] = NA
  parsed.x.agg$V.end[parsed.x.agg$RSS.V == 0] = NA
  parsed.x.agg$J.length = parsed.x.agg$L - parsed.x.agg$J.start +1
  if(correct.V.end)parsed.x.agg = df.V.end.adjuster(parsed.x.agg,v.end.alleles.list =as.character(c( get.V.end.Seq(species,gene,50,data.list = data.list,TCR.dict = TCR.dict)[names(TCR.dict[[species]][[gene]]$V.GR.dict)%in%names(TCR.dict[[species]][[gene]]$V.alleles)],TCR.dict[[species]][[gene]]$V.alleles)))
  #parsed.x.agg$RSS.V[is.na(parsed.x.agg$RSS.V)] = 0
  #parsed.x.agg$RSS.J[is.na(parsed.x.agg$RSS.J)] = 0
  #parsed.x.agg.df.2 = df.2.publicity.calculator(lapply(split(parsed.x.agg,parsed.x.agg$Sample.bio.name),function(x)x[!duplicated(x$pc.nt),]))
  return(parsed.x.agg)}

df.V.end.adjuster = function(df,v.end.alleles.list,ss = 25){
  df.original = df
  df = df[df$V.gene %in% names(v.end.alleles.list),]
  df$V.segment  = df$V.gene
  df$CDR3 = df$CDR3.nucleotide.sequence
  if(sum(df$L < ss)>0)df$CDR3[df$L < ss] = paste0(df$CDR3[df$L < ss],sapply(df$L[df$L < ss],function(x)paste(rep("X",x),collapse = "")))
  V.ends = SP.V.end_calling.mk2(df[df$V.gene %in% names(v.end.alleles.list),],SP.V.end.nt = subseq(v.end.alleles.list,1,ss),V.dict.names = unique(names(v.end.alleles.list)))
  V.sp = sort(unique(df$V.gene[df$V.gene%in%names(v.end.alleles.list)]))
  for(i in 1:length(V.sp))df[df$V.segment%in%V.sp[[i]],]$V.end = V.ends[[i]]
  df$V.segment = NULL
  df$CDR3 = NULL
  df$Total.insertions = df$J.start-(df$V.end+1)
  df.original[df.original$V.gene %in% names(v.end.alleles.list),] = df
  return(df.original)
}

df.D.adjuster = function(df,D.dictionary,min.D.l = 5,D.gene.names = c()){
  if (length(D.gene.names) != length(D.dictionary)) D.gene.names = 1:length(D.dictionary)
  if (length(D.dictionary)==1){
    D.s.p = D.segment.parser.mc(D.dictionary,df$CDR3.nucleotide.sequence,min.D.length = min.D.l)
    D.start = D.s.p$D.start
    D.end = D.s.p$D.end
    D.gene = df$D.gene
    D.gene[D.start != -1] = D.gene.names
    D.gene[D.start == -1] = -1
  }else{
    
    D.s.p = lapply(D.dictionary,function(d)D.segment.parser.mc(d,df$CDR3.nucleotide.sequence,min.D.length = min.D.l))
    D.s.p.cbind = do.call(cbind,lapply(D.s.p,function(x)x$D.length))
    D.s.p.max = apply(D.s.p.cbind[rowSums(D.s.p.cbind==-1)<2,],1,max)
    D.s.p.w.max = D.s.p.cbind[rowSums(D.s.p.cbind==-1)<2,]==D.s.p.max
    D.gene = rep(-1,length(D.s.p[[1]]$D.length))
    D.indices.1 = apply(D.s.p.w.max[rowSums(D.s.p.w.max)==1,],1,which)
    D.gene[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max)==1] = D.gene.names[D.indices.1]
    D.gene[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max) >1]  = apply(D.s.p.w.max[rowSums(D.s.p.w.max) >1,],1,function(x)paste(D.gene.names[x],collapse = ';'))
    
    D.start = rep(-1,length(D.s.p[[1]]$D.length))
    D.end = rep(-1,length(D.s.p[[1]]$D.length))
    
    for(i in unique(D.indices.1))D.start[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max)==1][D.indices.1==i] = D.s.p[[i]]$D.start[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max)==1][D.indices.1==i]
    for(i in unique(D.indices.1))D.end[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max)==1][D.indices.1==i] = D.s.p[[i]]$D.end[rowSums(D.s.p.cbind==-1)<2][rowSums(D.s.p.w.max)==1][D.indices.1==i]
  }
  df$D5.end = D.start
  df$D3.end = D.end
  df$VD.insertions[df$D5.end>0] = (df$D5.end - df$V.end)[df$D5.end>0]-1
  df$DJ.insertions[df$D5.end>0] = (df$J.start - df$D3.end)[df$D5.end>0]-1
  df$D.gene = D.gene
  
  return(df)
}

D.segment.parser.mc = function(spp.D.seq,cdr3s,min.D.length = 6,n.cores = 100){
  Db.kmers = sapply(min.D.length:nchar(spp.D.seq),function(i)kmer.generator(spp.D.seq,largo = i))
  #Db.matches = sapply(length(Db.kmers):1,function(i)sapply(1:length(Db.kmers[[i]]),function(j)regexpr(Db.kmers[[i]][j],cdr3s)))
  Db.matches = mclapply(length(Db.kmers):1,function(i)sapply(1:length(Db.kmers[[i]]),function(j)regexpr(Db.kmers[[i]][j],cdr3s)),mc.cores = n.cores)
  Db.indexes = mclapply(Db.matches,function(x)apply(x,1,function(y)if(sum(y!= -1) ==1) y[y!= -1]else 0),mc.cores = 100)
  
  #Db.indexes = lapply(Db.matches,function(x)apply(x,1,function(y)if(sum(y!= -1) ==1) y[y!= -1]else 0))
  names(Db.indexes) = nchar(spp.D.seq):min.D.length
  
  D.start = rep(-1,length(cdr3s))
  D.length = D.start
  D.end = D.start
  for (i in 1:(length(Db.indexes)))D.start[which(Db.indexes[[i]]!=0)] = Db.indexes[[i]][which(Db.indexes[[i]]!=0)]
  for (i in 1:(length(Db.indexes)))D.length[which(Db.indexes[[i]]!=0)] = as.numeric(names(Db.indexes[i]))
  D.end = D.length+D.start-1
  D.end[D.end == -3] = -1
  return(list(D.start = D.start,D.length = D.length,D.end = D.end))
}

rules.generation = function(TCR.a.df,sp,gene = "a"){
  V.seqs = get.V.end.Seq(sp,gene)
  V.cuts = spp.TCR.dict[[sp]][[gene]]$V.GR.dict$heptamer.start
  names(V.cuts) = names(V.seqs)
  J.seqs = get.J.Seq(sp,gene)[spp.TCR.dict[[sp]][[gene]]$J.GR.dict$heptamer.end < -50]
  J.cuts = spp.TCR.dict[[sp]][[gene]]$J.GR.dict[names(J.seqs)]$heptamer.end
  names(J.cuts)= names(J.seqs)
  
  
  rule.1.Vs = subseq(V.seqs[subseq(V.seqs,8,9)=="TG"],1,9)
  rule.1.Js = subseq(J.seqs,J.cuts+3,-34)[subseq(J.seqs,J.cuts+1,width = 2)=="TG"]
  rule.1 = unname(unlist(lapply(rule.1.Vs,function(x)paste0(x,rule.1.Js))))
  rule.1.V = paste(rule.1,rep(names(rule.1.Vs),each = length(rule.1.Js)),sep = ';')
  rule.1.VJ = paste(rule.1.V,rep(names(rule.1.Js),length(rule.1.Vs)),sep = ';')
  
  rule.2.Vs = subseq(V.seqs[subseq(V.seqs,12,12)=="G"],1,12)
  rule.2.Js = subseq(J.seqs,J.cuts+3,-34)[subseq(J.seqs,J.cuts+2,width = 1)=="G"]
  rule.2 = unname(unlist(lapply(rule.2.Vs,function(x)paste0(x,rule.2.Js))))
  rule.2.V = paste(rule.2,rep(names(rule.2.Vs),each = length(rule.2.Js)),sep = ';')
  rule.2.VJ = paste(rule.2.V,rep(names(rule.2.Js),length(rule.2.Vs)),sep = ';')
  
  rule.3.Vs = subseq(V.seqs[subseq(V.seqs,14,14)=="C"],1,14)
  rule.3.Js = subseq(J.seqs,J.cuts+5,-34)[subseq(J.seqs,J.cuts+4,width = 1)=="C"]
  rule.3 = unname(unlist(lapply(rule.3.Vs,function(x)paste0(x,rule.3.Js))))
  rule.3.V = paste(rule.3,rep(names(rule.3.Vs),each = length(rule.3.Js)),sep = ';')
  rule.3.VJ = paste(rule.3.V,rep(names(rule.3.Js),length(rule.3.Vs)),sep = ';')
  
  rule.4.Vs = subseq(V.seqs[subseq(V.seqs,14,14)=="C"],1,14)
  rule.4.Js = subseq(J.seqs,J.cuts+2,-34)[subseq(J.seqs,J.cuts+2,width = 1)=="G"]
  rule.4 = unname(unlist(lapply(rule.4.Vs,function(x)paste0(x,rule.4.Js))))
  rule.4.V = paste(rule.4,rep(names(rule.4.Vs),each = length(rule.4.Js)),sep = ';')
  rule.4.VJ = paste(rule.4.V,rep(names(rule.4.Js),length(rule.4.Vs)),sep = ';')
  
  rules = rep(-1,nrow(TCR.a.df))
  rules[TCR.a.df$pc.nt %in% rule.1.VJ] = 1
  rules[TCR.a.df$pc.nt %in% rule.2.VJ] = 2
  rules[TCR.a.df$pc.nt %in% rule.3.VJ] = 3
  rules[TCR.a.df$pc.nt %in% rule.4.VJ] = 4
  return (list(rules,list(rule.1.VJ,rule.2.VJ,rule.3.VJ,rule.4.VJ)))
}

D.locator.2 = function(df,D.names,D.seq,kmer.lengths = c(2),exclude = " "){
  res = rep(-1,nrow(df))
  
  for (Di in D.names){
    df.sel = df[with(df,D.gene== Di & D3.end > D5.end),]
    Db.sel = D.seq[[Di]]
    kmer.sel = lapply(kmer.lengths,function(i)kmer.generator(Db.sel,i)[[1]])
    kmer.sel = lapply(kmer.sel,function(x)x[!grepl(exclude,x)])
    kmer.sel.u = lapply(1:length(kmer.lengths),function(i)names(which(table(kmer.sel[[i]])==1)))
    kmer.pos.u = lapply(1:length(kmer.lengths),function(i)match(kmer.sel.u[[i]],kmer.sel[[i]]))
    Dseqs = with(df.sel,subseq(CDR3.nucleotide.sequence,D5.end,D3.end))
    D.c = rep(-1,length(Dseqs))
    for(i in 1:1)for(j in 1:length(kmer.sel.u[[i]])){
      rr.ij = regexpr(kmer.sel.u[[i]][[j]],Dseqs[D.c==-1])
      D.c[D.c==-1][rr.ij != -1] = rr.ij[rr.ij != -1] - kmer.pos.u[[i]][[j]] + 4 + df.sel[D.c==-1,][rr.ij != -1,]$D5.end 
    }
    res[df$D.gene== Di] = D.c
  }
  return(res)
}
