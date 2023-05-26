rules.generation = function(TCR.a.df,sp,gene = "a"){
  V.seqs = get.V.end.Seq(sp,gene)
  V.cuts = spp.TCR.dict[[sp]][[gene]]$V.GR.dict$heptamer.start
  names(V.cuts) = names(V.seqs)
  J.seqs = get.J.Seq(sp,gene)[spp.TCR.dict[[sp]][[gene]]$J.GR.dict$heptamer.end < -50]
  J.cuts = spp.TCR.dict[[sp]][[gene]]$J.GR.dict[names(J.seqs)]$heptamer.end
  names(J.cuts)= names(J.seqs)
  
  
  
  rule.1.Vs = subseq(V.seqs[subseq(V.seqs,8,9)=="TG"],1,9)
  if(!is.null(spp.TCR.dict[[sp]][[gene]]$V.alleles)){
    rule.1.Vs = c(rule.1.Vs,DNAStringSet(subseq(spp.TCR.dict[[sp]][[gene]]$V.alleles,1,9)))
  }
  rule.1.Js = subseq(J.seqs,J.cuts+3,-34)[subseq(J.seqs,J.cuts+1,width = 2)=="TG"]
  rule.1 = unname(unlist(lapply(rule.1.Vs,function(x)paste0(x,rule.1.Js))))
  rule.1.V = paste(rule.1,rep(names(rule.1.Vs),each = length(rule.1.Js)),sep = ';')
  rule.1.VJ = paste(rule.1.V,rep(names(rule.1.Js),length(rule.1.Vs)),sep = ';')
  
  rule.2.Vs = subseq(V.seqs[subseq(V.seqs,12,12)=="G"],1,12)
  if(!is.null(spp.TCR.dict[[sp]][[gene]]$V.alleles)){
    rule.2.Vs = c(rule.2.Vs,DNAStringSet(subseq(spp.TCR.dict[[sp]][[gene]]$V.alleles,1,12)))
  }
  rule.2.Js = subseq(J.seqs,J.cuts+3,-34)[subseq(J.seqs,J.cuts+2,width = 1)=="G"]
  rule.2 = unname(unlist(lapply(rule.2.Vs,function(x)paste0(x,rule.2.Js))))
  rule.2.V = paste(rule.2,rep(names(rule.2.Vs),each = length(rule.2.Js)),sep = ';')
  rule.2.VJ = paste(rule.2.V,rep(names(rule.2.Js),length(rule.2.Vs)),sep = ';')
  
  rule.3.Vs = subseq(V.seqs[subseq(V.seqs,14,14)=="C"],1,14)
  if(!is.null(spp.TCR.dict[[sp]][[gene]]$V.alleles)){
    rule.3.Vs = c(rule.3.Vs,DNAStringSet(subseq(spp.TCR.dict[[sp]][[gene]]$V.alleles,1,14)))
  }
  rule.3.Js = subseq(J.seqs,J.cuts+5,-34)[subseq(J.seqs,J.cuts+4,width = 1)=="C"]
  rule.3 = unname(unlist(lapply(rule.3.Vs,function(x)paste0(x,rule.3.Js))))
  rule.3.V = paste(rule.3,rep(names(rule.3.Vs),each = length(rule.3.Js)),sep = ';')
  rule.3.VJ = paste(rule.3.V,rep(names(rule.3.Js),length(rule.3.Vs)),sep = ';')
  
  rule.4.Vs = subseq(V.seqs[subseq(V.seqs,14,14)=="C"],1,14)
  if(!is.null(spp.TCR.dict[[sp]][[gene]]$V.alleles)){
    rule.4.Vs = c(rule.4.Vs,DNAStringSet(subseq(spp.TCR.dict[[sp]][[gene]]$V.alleles,1,14)))
  }
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


# Data frame extension

# Nucleotide removal - nibbling count 
{
  for (i in 1:11)TCR.a.spp[[i]]$V.nt.rem = TCR.a.spp[[i]]$RSS.V-TCR.a.spp[[i]]$V.end-1
  for (i in 1:11)TCR.a.spp[[i]]$V.nt.rem = TCR.a.spp[[i]]$V.nt.rem*as.numeric(TCR.a.spp[[i]]$V.nt.rem>=0)
  for (i in 1:11)TCR.a.spp[[i]]$J.nt.rem =-1*(TCR.a.spp[[i]]$RSS.J+TCR.a.spp[[i]]$J.length+34)
  for (i in 1:11)TCR.a.spp[[i]]$J.nt.rem = TCR.a.spp[[i]]$J.nt.rem*as.numeric(TCR.a.spp[[i]]$J.nt.rem>=0)
  
  for (i in 1:17)TCR.b.spp[[i]]$V.nt.rem = TCR.b.spp[[i]]$RSS.V-TCR.b.spp[[i]]$V.end-1
  for (i in 1:17)TCR.b.spp[[i]]$V.nt.rem = TCR.b.spp[[i]]$V.nt.rem*as.numeric(TCR.b.spp[[i]]$V.nt.rem>=0)
  for (i in 1:17)TCR.b.spp[[i]]$J.nt.rem = -1*(TCR.b.spp[[i]]$RSS.J+TCR.b.spp[[i]]$J.length+31)
  for (i in 1:17)TCR.b.spp[[i]]$J.nt.rem = TCR.b.spp[[i]]$J.nt.rem*as.numeric(TCR.b.spp[[i]]$J.nt.rem>=0)
  
  TCR.a.spp.VJ = lapply(TCR.a.spp,function(x)x[!is.na(with(x,RSS.V + RSS.J)),])
  TCR.b.spp.VJ = lapply(TCR.b.spp,function(x)x[!is.na(with(x,RSS.V + RSS.J)),])
}
TCR.b.spp.VJ.2 = TCR.b.spp.VJ
for(sp in names(TCR.b.spp.VJ.2)){
  if (sum(with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start)) >0 ){
    TCR.b.spp.VJ.2[[sp]][with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start),]$D.gene = -1
    TCR.b.spp.VJ.2[[sp]][with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start),]$VD.insertions = -1
    TCR.b.spp.VJ.2[[sp]][with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start),]$DJ.insertions = -1
    TCR.b.spp.VJ.2[[sp]][with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start),]$D3.end = -1
    TCR.b.spp.VJ.2[[sp]][with(TCR.b.spp.VJ.2[[sp]], D5.end > J.start),]$D5.end = -1
  }}
TCR.b.spp.VJ = TCR.b.spp.VJ.2
rm(TCR.b.spp.VJ.2)
# Recombination categories count

for(sp in 1:11){
  TCR.a.spp.VJ[[sp]]$rc.cat = 0
  TCR.a.spp.VJ[[sp]]$rc.cat[TCR.a.spp.VJ[[sp]]$V.nt.rem ==0 & TCR.a.spp.VJ[[sp]]$J.nt.rem ==0 &  TCR.a.spp.VJ[[sp]]$Total.insertions ==0] = 1
  
  TCR.a.spp.VJ[[sp]]$rc.cat[(TCR.a.spp.VJ[[sp]]$V.nt.rem >0 | TCR.a.spp.VJ[[sp]]$J.nt.rem >0) &  TCR.a.spp.VJ[[sp]]$Total.insertions <= 0] = 2
  
  TCR.a.spp.VJ[[sp]]$rc.cat[TCR.a.spp.VJ[[sp]]$Total.insertions < 0] = 2
  TCR.a.spp.VJ[[sp]]$rc.cat[TCR.a.spp.VJ[[sp]]$V.nt.rem ==0 & TCR.a.spp.VJ[[sp]]$J.nt.rem ==0 &  TCR.a.spp.VJ[[sp]]$Total.insertions > 0] = 3
  TCR.a.spp.VJ[[sp]]$rc.cat[(TCR.a.spp.VJ[[sp]]$V.nt.rem >0 | TCR.a.spp.VJ[[sp]]$J.nt.rem >0) &  TCR.a.spp.VJ[[sp]]$Total.insertions > 0] = 4
}

for(sp in 1:11) TCR.a.spp.VJ[[sp]]$rules = rules.generation(TCR.a.spp.VJ[[sp]],gsub("[^A-Z]","",names(TCR.a.spp.VJ)[[sp]]),c(rep("a",6),"a.1","a.2",rep("a",3))[[sp]])[[1]]



vj.a.table.sp = lapply(names(TCR.a.spp),function(x)with(TCR.a.spp[[x]][with(TCR.a.spp[[x]], !is.na(RSS.V) & !is.na(RSS.J) & !duplicated(pc.nt)),],table(paste(V.gene,J.gene))))
names(vj.a.table.sp) = names(TCR.a.spp)

vj.b.table.sp = lapply(names(TCR.b.spp),function(x)with(TCR.b.spp[[x]][with(TCR.b.spp[[x]], !is.na(RSS.V) & !is.na(RSS.J) & !duplicated(pc.nt)),],table(paste(V.gene,J.gene))))
names(vj.b.table.sp) = names(TCR.b.spp)


# Microhomology
{
  MH.rule1 = do.call(rbind,lapply(names(TCR.a.spp.VJ),function(sp){(with(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$L > 20 ,][with(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$L > 20 ,], Total.insertions <= -2 & subseq(CDR3.nucleotide.sequence,8,9) == "TG" & -33 - RSS.J - J.length == 1),],table(factor(rep(Total.insertions,Umi.count),levels = -2:-4)))/ sum(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$L > 20 ,]$Umi.count))}))
  MH.rule1 = MH.rule1[1:11,]
  MH.rule1.df.1 = data.frame(MH.rule1)
  MH.rule1.df = cbind(sp = names(TCR.a.spp.VJ),MH.rule1.df.1)
  colnames(MH.rule1.df) = c("sp",-2:-4)
  MH.rule1.gather = gather(MH.rule1.df,"MH","Umi",2:4)
  MH.rule1.gather$sp = factor(MH.rule1.gather$sp,levels =names(TCR.a.spp.VJ)[c(6,10,7,8,1:3,9,11,4,5)])
}

# Repertoire data frame: spp.repertoire.df  
{
  
  spp.repertoire.a.df = data.frame("UMIs" = sapply(TCR.a.spp.VJ,function(x)sum(x$Umi.count)),
                                   "clones" = sapply(TCR.a.spp.VJ,function(x)length(unique(x$pc.nt))),
                                   "VJ combinations" = sapply(TCR.a.spp.VJ,function(sp)with(sp,length(unique(paste(V.gene,J.gene))))),
                                   "Identified.V" = round(sapply(TCR.a.spp,function(sp)with(sp,mean(rep(!is.na(RSS.V),Umi.count)))),3),
                                   "Identified.J" = round(sapply(TCR.a.spp,function(sp)with(sp,mean(rep(!is.na(RSS.J),Umi.count)))),3),
                                   "Identified.VJ" = round(sapply(TCR.a.spp,function(sp)with(sp,mean(rep(!is.na(RSS.V + RSS.J),Umi.count)))),3),
                                   "Mean length" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(L,Umi.count)))),1),
                                   "No insertions" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(Total.insertions <= 0 ,Umi.count)))),3),
                                   "MH" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(Total.insertions <= -1 ,Umi.count)))),3),
                                   "Entropy" = t(sapply(names(TCR.a.spp),function(x)round(c(unlist(h.spp.a[[x]]$statistics[1:2])),1))),
                                   "H.V" = do.call(rbind,ger.spp.a)[,2],
                                   "H.J" = do.call(rbind,ger.spp.a)[,3],
                                   "RSS.V.RF3" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(RSS.V%%3==0,Umi.count)))),3),
                                   "RSS.J.RF1" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(RSS.J%%3==0,Umi.count)))),3),
                                   "In.frame" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(L%%3==0,Umi.count)))),3),
                                   "V.nt.rem" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(V.nt.rem,Umi.count)))),2),
                                   "J.nt.rem" = round(sapply(TCR.a.spp.VJ,function(sp)with(sp,mean(rep(J.nt.rem,Umi.count)))),2))
  
  spp.repertoire.a.df = cbind(sp = rownames(spp.repertoire.a.df),spp.repertoire.a.df)
  spp.repertoire.a.df$H.N = spp.repertoire.a.df$Entropy.entropy.nt.VorJ - spp.repertoire.a.df$H.V - spp.repertoire.a.df$H.J - do.call(rbind,ger.spp.a)[,1] + do.call(rbind,ger.spp.a)[,4]
  rownames(spp.repertoire.a.df) = NULL
  spp.repertoire.a.df = spp.repertoire.a.df[c(6,10,7,8,1:3,9,11,4,5),]
  spp.repertoire.a.df$class = c(rep("sharks",1),rep("basal fishes",3),rep("teleosts",4),rep("lobe-finned fishes",1),rep("mammals",2))
  
  
  
  spp.repertoire.b.df = data.frame("UMIs" = sapply(TCR.b.spp.VJ,function(x)sum(x$Umi.count)),
                                   "clones" = sapply(TCR.b.spp.VJ,function(x)length(unique(x$pc.nt))),
                                   "VJ combinations" = sapply(TCR.b.spp.VJ,function(sp)with(sp,length(unique(paste(V.gene,J.gene))))),
                                   "Identified.V" = round(sapply(TCR.b.spp,function(sp)with(sp,mean(rep(!is.na(RSS.V),Umi.count)))),3),
                                   "Identified.J" = round(sapply(TCR.b.spp,function(sp)with(sp,mean(rep(!is.na(RSS.J),Umi.count)))),3),
                                   "Identified.VJ" = round(sapply(TCR.b.spp,function(sp)with(sp,mean(rep(!is.na(RSS.V + RSS.J),Umi.count)))),3),
                                   "Mean length" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(L,Umi.count)))),1),
                                   "No VD insertions" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp[sp$D.gene %in% 1:2,],mean(rep(VD.insertions <= 0 ,Umi.count)))),3),
                                   "No DJ insertions" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp[sp$D.gene %in% 1:2,],mean(rep(DJ.insertions <= 0 ,Umi.count)))),3),
                                   "No VDJ insertions" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp[sp$D.gene %in% 1:2,],mean(rep(VD.insertions <= 0 & DJ.insertions <= 0,Umi.count)))),3),
                                   "Entropy" = rbind(t(sapply(1:12,function(x)round(c(unlist(h.spp.b[[x]]$statistics[1:2])),1))),c(0,0),t(sapply(14:17,function(x)round(c(unlist(h.spp.b[[x]]$statistics[1:2])),1)))),
                                   "H.V" = do.call(rbind,ger.spp.b)[,2],
                                   "H.J" = do.call(rbind,ger.spp.b)[,3],
                                   "RSS.V.RF3" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(RSS.V%%3==0,Umi.count)))),3),
                                   "RSS.J.RF1" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(RSS.J%%3==0,Umi.count)))),3),
                                   "In.frame" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(L%%3==0,Umi.count)))),3),
                                   "V.nt.rem" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(V.nt.rem,Umi.count)))),2),
                                   "J.nt.rem" = round(sapply(TCR.b.spp.VJ,function(sp)with(sp,mean(rep(J.nt.rem,Umi.count)))),2))
  
  spp.repertoire.b.df = cbind(sp = rownames(spp.repertoire.b.df),spp.repertoire.b.df)
  spp.repertoire.b.df$H.N = spp.repertoire.b.df$Entropy.entropy.nt.VorJ - spp.repertoire.b.df$H.V - spp.repertoire.b.df$H.J - do.call(rbind,ger.spp.b)[,1] + do.call(rbind,ger.spp.b)[,4]
  spp.repertoire.b.df$H.N [spp.repertoire.b.df$H.N <0] = 0
  rownames(spp.repertoire.b.df) = NULL
  spp.repertoire.b.df = spp.repertoire.b.df[c(6:7,11:15,8:9,1:3,10,16:17,4,5),]
  spp.repertoire.b.df$class = c(rep("sharks",2),rep("basal fishes",7),rep("teleosts",4),rep("lobe-finned fishes",2),rep("mammals",2))
  
  
}

write.csv2(spp.repertoire.a.df,paste0(spp.tables.folder,"/spp.repertoire.a.df.csv"))
write.csv2(spp.repertoire.b.df,paste0(spp.tables.folder,"/spp.repertoire.b.df.csv"))



# Mutants data frames

# ZF.Ca.mut.df.list = list(ZF.a2.e7.dfl,ZF.a3.e7.dfl,ZF.a.mut.a3.dfl.MUT,ZF.a.mut.a23.dfl.MUT,ZF.a.mut.a35.dfl.MUT)
#  ZF.Ca.wt.df.list = list(ZF.a1.e7.dfl,ZF.a.mut.a3.dfl.WT,ZF.a.mut.a23.dfl.WT,ZF.a.mut.a35.dfl.WT)

{
  ZF.Ca.df.mut = data.frame(s = c("E7_frame_2","E7_frame_3","A3_Mut","A23_Mut","A35_Mut"),in.frame = sapply(ZF.Ca.mut.df.list,function(x)mean(x$L%%3==0)),genotype = "MUT",sp = "ZF")
  ZF.Ca.df.wt = data.frame(s = c("E7_frame_1","A3_WT","A23_WT","A35_WT"),in.frame = sapply(ZF.Ca.wt.df.list,function(x)mean(x$L%%3==0)),genotype = "WT",sp = "ZF")
  MM.Ca.df.mut = data.frame(s = unique(MM.a.mut.a32.dfl.MUT$Sample.bio.name),in.frame = tapply(MM.a.mut.a32.dfl.MUT$L,MM.a.mut.a32.dfl.MUT$Sample.bio.name,function(x)mean(x%%3==0)),genotype = "MUT",sp = "MM")
  MM.Ca.df.wt = data.frame(s = unique(TCR.a.spp$MM$Sample.bio.name),in.frame = tapply(TCR.a.spp$MM$L,TCR.a.spp$MM$Sample.bio.name,function(x)mean(x%%3==0)),genotype = "WT",sp = "MM")
  Ca.df = rbind(ZF.Ca.df.mut,ZF.Ca.df.wt,MM.Ca.df.mut,MM.Ca.df.wt)
  Ca.df$genotype = factor(Ca.df$genotype, levels = c("WT","MUT"))
  Ca.df$in.frame = round(Ca.df$in.frame*100,2)
  
  ZF.Ca.df.mut.umi = data.frame(s = c("E7_frame_2","E7_frame_3","A3_Mut","A23_Mut","A35_Mut"),in.frame = sapply(ZF.Ca.mut.df.list,function(x)mean(rep(x$L,x$Umi.count)%%3==0)),genotype = "MUT",sp = "ZF")
  ZF.Ca.df.wt.umi = data.frame(s = c("E7_frame_1","A3_WT","A23_WT","A35_WT"),in.frame = sapply(ZF.Ca.wt.df.list,function(x)mean(rep(x$L,x$Umi.count)%%3==0)),genotype = "WT",sp = "ZF")
  MM.Ca.df.mut.umi = data.frame(s = unique(MM.a.mut.a32.dfl.MUT$Sample.bio.name),in.frame = tapply(rep(MM.a.mut.a32.dfl.MUT$L,MM.a.mut.a32.dfl.MUT$Umi.count),rep(MM.a.mut.a32.dfl.MUT$Sample.bio.name,MM.a.mut.a32.dfl.MUT$Umi.count),function(x)mean(x%%3==0)),genotype = "MUT",sp = "MM")
  MM.Ca.df.wt.umi = data.frame(s = unique(TCR.a.spp$MM$Sample.bio.name),in.frame = tapply(rep(TCR.a.spp$MM$L,TCR.a.spp$MM$Umi.count),rep(TCR.a.spp$MM$Sample.bio.name,TCR.a.spp$MM$Umi.count),function(x)mean(x%%3==0)),genotype = "WT",sp = "MM")
  Ca.df.umi = rbind(ZF.Ca.df.mut.umi,ZF.Ca.df.wt.umi,MM.Ca.df.mut.umi,MM.Ca.df.wt.umi)
  Ca.df.umi$genotype = factor(Ca.df$genotype, levels = c("WT","MUT"))
  Ca.df.umi$in.frame = round(Ca.df.umi$in.frame*100,2)
  df.Ca.mut = do.call(rbind,ZF.Ca.mut.df.list)[with(do.call(rbind,ZF.Ca.mut.df.list),!is.na(RSS.V) & !is.na(RSS.J)),]
  df.Ca.wt = do.call(rbind,ZF.Ca.wt.df.list)[with(do.call(rbind,ZF.Ca.wt.df.list),!is.na(RSS.V) & !is.na(RSS.J)),]
  df.Ca.mut$V.nt.rem = df.Ca.mut$RSS.V-df.Ca.mut$V.end-1
  df.Ca.mut$V.nt.rem = df.Ca.mut$V.nt.rem*as.numeric(df.Ca.mut$V.nt.rem >=0)
  df.Ca.mut$J.nt.rem =-1*(df.Ca.mut$RSS.J+df.Ca.mut$J.length+34)  
  df.Ca.mut$J.nt.rem = df.Ca.mut$J.nt.rem*as.numeric(df.Ca.mut$J.nt.rem >=0)
  
  df.Ca.mut$rc.cat = 0
  df.Ca.mut$rc.cat[ df.Ca.mut$V.nt.rem ==0 &  df.Ca.mut$J.nt.rem ==0 &   df.Ca.mut$Total.insertions ==0] = 1
  
  df.Ca.mut$rc.cat[( df.Ca.mut$V.nt.rem >0 |  df.Ca.mut$J.nt.rem >0) &   df.Ca.mut$Total.insertions <= 0] = 2
  
  df.Ca.mut$rc.cat[ df.Ca.mut$Total.insertions < 0] = 2
  df.Ca.mut$rc.cat[ df.Ca.mut$V.nt.rem ==0 &  df.Ca.mut$J.nt.rem ==0 &   df.Ca.mut$Total.insertions > 0] = 3
  df.Ca.mut$rc.cat[( df.Ca.mut$V.nt.rem >0 |  df.Ca.mut$J.nt.rem >0) &   df.Ca.mut$Total.insertions > 0] = 4
  
  
}



{
  df.ZF.Va.RSS = data.frame(round(prop.table(table(factor(with(TCR.a.spp.VJ$ZF,rep(RSS.V,Umi.count)),levels = 12:21))),3))
  names(df.ZF.Va.RSS)[1] = c("RSS")
  df.ZF.Va.RSS.germline = data.frame(table(factor(spp.TCR.dict$ZF$a$V.GR.dict$heptamer.start, levels = 12:21)))
  names(df.ZF.Va.RSS.germline) = names(df.ZF.Va.RSS)
  df.ZF.Ja.RSS = data.frame(round(prop.table(table(factor(rep(TCR.a.spp$ZF$RSS.J[!is.na(TCR.a.spp$ZF$RSS.J)],TCR.a.spp$ZF$Umi.count[!is.na(TCR.a.spp$ZF$RSS.J)]),levels = -57:-69))),3))
  names(df.ZF.Ja.RSS)[1] = c("RSS")
  df.ZF.Ja.RSS.germline = data.frame(table(factor(spp.TCR.dict$ZF$a$J.GR.dict$heptamer.end,levels = -57:-69)))
  names(df.ZF.Ja.RSS.germline) = names(df.ZF.Ja.RSS)
  
  df.ZF.RSS.list = list(df.ZF.Va.RSS,df.ZF.Va.RSS.germline,df.ZF.Ja.RSS,df.ZF.Ja.RSS.germline)
}


{
  sp = "ZF"
  df= TCR.a.spp.VJ[[sp]][with(TCR.a.spp.VJ[[sp]], L>=27 & L<=57 & RSS.J %in% seq(-57,-69,by=-3)),]
  df.sel = df[,with(df,c("Umi.count","RSS.V","RSS.J","L"))]
  temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
  temp.agg$RSS.J = factor(temp.agg$RSS.J,levels = c(seq(-57,-69,by = -3)))
}

# Rules dataframes
TCR.a.spp.rules = TCR.a.spp.VJ

rules.explainable = sapply(names(TCR.a.spp.VJ),function(sp)mean(rep(TCR.a.spp.VJ[[sp]]$Total.insertions,TCR.a.spp.VJ[[sp]]$Umi.count)<=0))
rules.spp = lapply(names(TCR.a.spp),function(sp)prop.table(tapply(TCR.a.spp.rules[[sp]][with(TCR.a.spp.rules[[sp]], !is.na(RSS.V) & !is.na(RSS.J)),]$Umi.count,(TCR.a.spp.rules[[sp]][with(TCR.a.spp.rules[[sp]], !is.na(RSS.V) & !is.na(RSS.J)),]$rules),sum)))
for(i in 1:11) TCR.a.spp.rules[[i]]$rules[with(TCR.a.spp.rules[[i]],rules %in% 1:3 & Total.insertions==0)] = 0
rules.tb.1 = lapply(1:11,function(i) sort(prop.table(table(rep(TCR.a.spp.rules[[i]]$rules[!is.na(TCR.a.spp.rules[[i]]$RSS.V) & !is.na(TCR.a.spp.rules[[i]]$RSS.J)],TCR.a.spp.rules[[i]]$Umi.count[!is.na(TCR.a.spp.rules[[i]]$RSS.V) & !is.na(TCR.a.spp.rules[[i]]$RSS.J)]))),decreasing = T))
rules.tb =sapply(1:11,function(i)rules.tb.1[[i]][as.character(c(1,2,3,4))])  
colnames(rules.tb) = names(TCR.a.spp)
rules.tb.cumsum = apply(rules.tb,2,cumsum)

rules.df = cbind(gather(data.frame(rules.tb),key = "sp"),rules = 1:4)
rules.df = cbind(rules.df,groups = spp.repertoire.a.df$class[match(rules.df$sp,spp.repertoire.a.df$sp)])
rules.df$groups = factor(rules.df$groups,levels = unique(rules.df$groups)[c(3,4,1,5,2)])
rules.df$explainable = as.numeric(rules.df$value/rules.explainable[rules.df$sp])

rules.df.c = cbind(gather(data.frame(rules.tb.cumsum),key = "sp"),rules = 1:4)
rules.df.c = cbind(rules.df.c,groups = spp.repertoire.a.df$class[match(rules.df$sp,spp.repertoire.a.df$sp)])
rules.df.c$groups = factor(rules.df.c$groups,levels = unique(rules.df.c$groups)[c(3,4,1,5,2)])

rules.ZF.df = data.frame(rule = c(paste0("rule",1:4),'unkonwn'), freq = c(rules.tb[,2],1-sum(rules.tb[,2])))

# In frameness V
{
  df.V.frame.df = data.frame(sort(tapply(with(TCR.a.spp$ZF[!is.na(TCR.a.spp$ZF$RSS.V),],rep(L,Umi.count) %% 3==0),with(TCR.a.spp$ZF[!is.na(TCR.a.spp$ZF$RSS.V),],rep(V.gene,Umi.count)),mean)))
  df.V.frame.df$Vs = rownames(df.V.frame.df)
  rownames(df.V.frame.df) = NULL
  colnames(df.V.frame.df)[1] = "in.frame"
  df.V.frame.df = df.V.frame.df[,2:1]
  df.V.frame.df = cbind(df.V.frame.df, RSS = spp.TCR.dict$ZF$a$V.GR.dict[df.V.frame.df$Vs]$heptamer.start)
}

# In frameness J  
{
  df.J.frame.df = data.frame(sort(tapply(with(TCR.a.spp$ZF[!is.na(TCR.a.spp$ZF$RSS.J),],rep(L,Umi.count) %% 3==0),with(TCR.a.spp$ZF[!is.na(TCR.a.spp$ZF$RSS.J),],rep(J.gene,Umi.count)),mean)))
  df.J.frame.df$Js = rownames(df.J.frame.df)
  rownames(df.J.frame.df) = NULL
  colnames(df.J.frame.df)[1] = "in.frame"
  df.J.frame.df = df.J.frame.df[,2:1]
  df.J.frame.df = cbind(df.J.frame.df, RSS = spp.TCR.dict$ZF$a$J.GR.dict[df.J.frame.df$Js]$heptamer.end)
}


{
  J5.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$J.gene == "TRAJ5_01",],rep(L,Umi.count)),levels = 30:55))
  V52.1.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "TRAV52_1_01",],rep(L,Umi.count)),levels = 30:55))
  V11.6.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "TRAV11_6_01",],rep(L,Umi.count)),levels = 30:55))
  nV1.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "new.Va_01",],rep(L,Umi.count)),levels = 30:55))
  V65.2.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "TRAV65_2_01",],rep(L,Umi.count)),levels = 30:55))
  J5.V11.6.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "TRAV11_6_01" & TCR.a.spp$ZF$J.gene == "TRAJ5_01",],rep(L,Umi.count)),levels = 30:55))
  J5.V52.1.L.table = table(factor(with(TCR.a.spp$ZF[TCR.a.spp$ZF$V.gene == "TRAV52_1_01" & TCR.a.spp$ZF$J.gene == "TRAJ5_01",],rep(L,Umi.count)),levels = 30:55))
}


{
  spp.V.frame.df = data.frame(sapply(1:11,function(i)prop.table(table(factor(rep(TCR.a.spp[[i]]$RSS.V-1,TCR.a.spp[[i]]$Umi.count)%%3,levels = 0:2)))))
  colnames(spp.V.frame.df) = names(TCR.a.spp)
  spp.V.frame.df = cbind(gather(spp.V.frame.df,key = "sp"),frame = 0:2)
  spp.V.frame.df$sp = factor(spp.V.frame.df$sp, levels = unique(spp.V.frame.df$sp)[c(6,10,7,8,1:3,9,11,4,5)])
  spp.V.frame.df$frame = factor(spp.V.frame.df$frame)
}
{
  spp.J.frame.df = data.frame(sapply(1:11,function(i)prop.table(table(factor(rep(TCR.a.spp[[i]]$RSS.J,TCR.a.spp[[i]]$Umi.count)%%3,levels = 0:2)))))
  colnames(spp.J.frame.df) = names(TCR.a.spp)
  spp.J.frame.df = cbind(gather(spp.J.frame.df,key = "sp"),frame = 0:2)
  spp.J.frame.df$sp = factor(spp.J.frame.df$sp, levels = unique(spp.J.frame.df$sp)[c(6,10,7,8,1:3,9,11,4,5)])
  spp.J.frame.df$frame = factor(spp.J.frame.df$frame)
}

spp.categories = spp.repertoire.a.df$class[c(1,2,6,9,10)]

{
  spp.RSS.V.tb = lapply(1:11,function(i)data.frame(round(prop.table(with(TCR.a.spp.VJ[[i]],table(factor(rep(RSS.V,Umi.count),levels = 1:25)))),3),names(TCR.a.spp.VJ)[i]))
  for (i in 1:11)colnames(spp.RSS.V.tb[[i]]) = c("RSS","fraction","sp")
  spp.RSS.V.df = do.call(rbind,spp.RSS.V.tb)
  spp.RSS.V.df = cbind(spp.RSS.V.df,class = spp.repertoire.a.df$class[match(spp.RSS.V.df$sp,spp.repertoire.a.df$sp)])
  spp.RSS.V.df.by.class <- spp.RSS.V.df %>% 
    group_by(class, RSS) %>% 
    dplyr::summarise(fraction = sum(fraction))
}
{
  spp.RSS.J.tb = lapply(1:11,function(i)data.frame(round(prop.table(with(TCR.a.spp.VJ[[i]],table(factor(rep(RSS.J,Umi.count),levels = -78:-50)))),3),names(TCR.a.spp.VJ)[i]))
  for (i in 1:11)colnames(spp.RSS.J.tb[[i]]) = c("RSS","fraction","sp")
  spp.RSS.J.df = do.call(rbind,spp.RSS.J.tb)
  spp.RSS.J.df = cbind(spp.RSS.J.df,class = spp.repertoire.a.df$class[match(spp.RSS.J.df$sp,spp.repertoire.a.df$sp)])
  spp.RSS.J.df.by.class <- spp.RSS.J.df %>% 
    group_by(class, RSS) %>% 
    dplyr::summarise(fraction = sum(fraction))
}

