# Fig A
names.pairs.list.col = lapply(1:11,function(i)list(names(TCR.a.spp.VJ)[[i]],c(rep("a",6),"a.1","a.2",rep("a",3))[[i]],c(3,3,3,5,5,1,2,2,3,2,4)[[i]]))[c(6,10,7,8,1:3,9,11,4,5)]

weighted.RSS.V.umi.col = function(sp,gene,col.n,posicion = 6, largo = 14){
  spx = gsub("[^A-Z]","",sp)
  filtro = rep(TCR.a.spp.VJ[[sp]]$V.gene,TCR.a.spp.VJ[[sp]]$Umi.count)[rep(TCR.a.spp.VJ[[sp]]$V.gene,TCR.a.spp.VJ[[sp]]$Umi.count) %in% names(get.V.end.Seq(spx,gene))]
  
  colores = c('darkseagreen 2','cadetblue2',"royal blue",'dark orchid',"red")
  mcs = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(colores[as.numeric(col.n)],4))  
  p = ggseqlogo(consensusMatrix(subseq(get.V.end.Seq(spx,gene)[filtro],1,width = largo))[1:4,] , col_scheme = mcs) + 
    theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(), axis.title.y = element_blank(),axis.text.y=element_blank()) 
  #ggtitle(sp)
  return(p)}


RSS.V.umi.grobs = lapply(c(1:11),function(i)do.call(weighted.RSS.V.umi.col,names.pairs.list.col[[i]]))

pdf(paste0(spp.fig.folder,"Fig3.A_V.pdf"),width = 3,height = 6)
grid.arrange(grobs = RSS.V.umi.grobs,ncol = 1)
dev.off()

weighted.RSS.J.umi.col = function(sp,gene,col.n,posicion = -1+3, largo = 5+3){
  spx = gsub("[^A-Z]","",sp)
  filtro = rep(TCR.a.spp[[sp]]$J.gene,TCR.a.spp[[sp]]$Umi.count)[rep(TCR.a.spp[[sp]]$J.gene,TCR.a.spp[[sp]]$Umi.count) %in% names(get.J.Seq(spx,gene))]
  colores = c('darkseagreen 2','cadetblue2',"royal blue",'dark orchid',"red")
  mcs = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(colores[as.numeric(col.n)],4))  
  p = ggseqlogo(consensusMatrix(subseq(get.J.Seq(spx,gene)[filtro],spp.TCR.dict[[spx]][[gene]]$J.GR.dict[filtro]$heptamer.end-posicion,width = largo))[1:4,], col_scheme = mcs) + 
    theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(), axis.title.y = element_blank(),axis.text.y=element_blank()) 
  
  return(p)}



RSS.J.umi.grobs = lapply(c(1:11),function(i)do.call(weighted.RSS.J.umi.col,names.pairs.list.col[[i]]))

pdf(paste0(spp.fig.folder,"Fig3.A_J.pdf"),width = 2,height = 6)
grid.arrange(grobs = RSS.J.umi.grobs,ncol = 1)
dev.off()



# Fig 3B

MH.rule1.gather.2 = MH.rule1.gather
MH.rule1.gather.2$MH = factor(-as.numeric(MH.rule1.gather.2$MH))

p.rule1 = ggplot(MH.rule1.gather.2) + geom_col(mapping = aes(x= sp,y = Umi,fill = sp), color = 'black') + 
  
  scale_fill_manual(values = rep(spp.col,c(1,3,4,1,2)),labels = c("4 nucleotides", "3 nucleotides", "2 nucleotides")) +
  scale_x_discrete(labels = c("C. punctatum","P. senegalus","A. ruthenus Ca 1","A. ruthenus Ca 2","P. progenetica","D. rerio","O. mykiss","O. latipes","P. annectens","M. musculus","L. africana"),) +
  labs(y = "Fraction of the repertoire", x = "Species", fill = "Microhomology length") + 
  
  theme_minimal() +
  geom_text(aes(label = MH,x= sp,y = Umi), color = "black", size = 6, position = position_stack(vjust = 0.5)) + 
  
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
    axis.text.y = element_text(size = 12))

ggsave(paste0(spp.fig.folder,"Fig3.B.pdf"),p.rule1,width = 9,height = 7)


# Fig 3C


p.HvsI = ggplot(spp.repertoire.a.df,aes(x = as.numeric(No.insertions), y = as.numeric(H.N))) + 
  scale_y_continuous(limits = c(-4,25)) +
  scale_color_manual(values = spp.col[c(2,4,5,1,3)]) +
  scale_fill_manual(values = spp.col[c(2,4,5,1,3)]) +
  theme_linedraw() +
  stat_summary(fun.data=mean_cl_normal) + 
  
  geom_smooth(method='lm', formula= y~x) +
  geom_label_repel(mapping = aes( fill = class,label = spp.short.names),size =8,fontface = 'italic') +
  labs(x = "Repertoire fraction with no insertions", y = "CDR3 V and J independent entropy (bits)") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=25),
        legend.position = "none") 

p.HvsI
ggsave(paste0(spp.fig.folder,"Fig3.C.pdf"),p.HvsI,width = 9,height = 9)


# Fig 3D
spp.categories = unique(spp.repertoire.a.df$class)
{
  
  
  spp.RSS.V.df.plots = lapply(1:5,function(i){
    ggplot(spp.RSS.V.df.by.class[spp.RSS.V.df.by.class$class == spp.categories[i],]) + 
      geom_col(mapping = aes(x=RSS,y=fraction/c(1,3,4,1,2)[i]), color = 'black',fill = spp.col[i],position = "dodge") +
      theme_classic() +
      scale_x_discrete(labels = c(rep('',14),15,rep('',10))) +
      #scale_y_continuous(limits = c(0,1.05)) +
      labs(x = element_blank(), y = element_blank()) + 
      theme(axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25)
            
      )
  })
  pdf(paste0(spp.fig.folder,"Fig3.D.Vs.pdf"),width = 24,height = 4)
  grid.arrange(grobs = spp.RSS.V.df.plots,nrow = 1)
  dev.off()
  
  spp.RSS.J.df.by.class.2 = spp.RSS.J.df.by.class
  spp.RSS.J.df.by.class.2$RSS = factor(-as.numeric(as.character(spp.RSS.J.df.by.class.2$RSS)),levels = sort(unique(-as.numeric(as.character(spp.RSS.J.df.by.class.2$RSS))))) 
  spp.RSS.J.df.plots = lapply(1:5,function(i){
    ggplot(spp.RSS.J.df.by.class.2[spp.RSS.J.df.by.class.2$class == spp.categories[i],]) + 
      scale_x_discrete(breaks = seq(0, 76, by = 6)) +
      geom_col(mapping = aes(x=RSS,y=fraction/c(1,3,4,1,2)[i]), color = 'black',fill = spp.col[i],position = "dodge") +
      theme_classic() +
      #scale_y_continuous(limits = c(0,0.6)) +
      labs(x = element_blank(), y = element_blank()) + 
      
      #labs(x = "J RSS position", y = "Fraction of the repertoire") + 
      theme(axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title=element_text(size=18))
  })
  pdf(paste0(spp.fig.folder,"Fig3.D.Js.pdf"),width = 24,height = 4)
  grid.arrange(grobs = spp.RSS.J.df.plots,nrow = 1)
  dev.off()
}

####
{
  spp.names = unique(sub("\\..*","",names(TCR.a.spp.VJ)))
  spp.a.df = lapply(spp.names,function(x)do.call(rbind,TCR.a.spp.VJ[grepl(x,names(TCR.a.spp.VJ))]))
  names(spp.a.df) = spp.names
  spp.a.ins = lapply(spp.a.df,function(x)with(x,rep(Total.insertions,Umi.count)))
  spp.a.xy = rbind(sapply(spp.a.ins,function(x)mean(x[x>0])*mean(x>0)),sapply(spp.a.ins,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]
  rbind(sapply(spp.a.ins,function(x)mean(x[x<0])*mean(x<0)),sapply(spp.a.ins,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]
  
}
####
{
spp.names = unique(sub("\\..*","",names(TCR.b.spp.VJ)))
spp.b.df = lapply(spp.names,function(x)do.call(rbind,TCR.b.spp.VJ[grepl(x,names(TCR.b.spp.VJ))]))
#spp.b.df = lapply(spp.names,function(x)do.call(rbind,TCR.b.spp.VJ.2[grepl(x,names(TCR.b.spp.VJ.2))]))
names(spp.b.df) = spp.names
spp.b.dj = lapply(spp.b.df,function(x)with(x[x$D.gene %in% 1:2,],rep(DJ.insertions,Umi.count)))
spp.b.vd = lapply(spp.b.df,function(x)with(x[x$D.gene %in% 1:2,],rep(VD.insertions,Umi.count)))
spp.b.x = rbind(sapply(spp.b.vd,function(x)mean(x[x>0])*mean(x>0)),sapply(spp.b.dj,function(x)mean(x[x>0])*mean(x>0)))[,c(6,9,7,1,2,3,8,10,4,5)]
spp.b.y = rbind(sapply(spp.b.vd,function(x)mean(x[x<0])*mean(x<0)),sapply(spp.b.dj,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]
}
spp.a.ins.df = data.frame(t(rbind(sapply(spp.a.ins,function(x)mean(x[x>0])*mean(x>0)),-sapply(spp.a.ins,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]))
colnames(spp.a.ins.df) = c("Ins","MH")
spp.a.ins.df$class = factor(spp.repertoire.a.df$class[-4],levels = unique(spp.repertoire.a.df$class))

spp.b.vd.df = data.frame(t(rbind(sapply(spp.b.vd,function(x)mean(x[x>0])*mean(x>0)),-sapply(spp.b.vd,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]))
colnames(spp.b.vd.df) = c("Ins","MH")
spp.b.vd.df$class =spp.a.ins.df$class
spp.b.dj.df = data.frame(t(rbind(sapply(spp.b.dj,function(x)mean(x[x>0])*mean(x>0)),-sapply(spp.b.dj,function(x)mean(x[x<0])*mean(x<0)))[,c(6,9,7,1,2,3,8,10,4,5)]))
colnames(spp.b.dj.df) = c("Ins","MH")
spp.b.dj.df$class =spp.a.ins.df$class

{
  p1=ggplot(spp.a.ins.df) + geom_label_repel(mapping = aes(x = Ins,y = MH, fill = class, label = spp.short.names[-4]),fontface = 'italic') + 
    geom_point(mapping = aes(x = Ins,y = MH)) +
    scale_fill_manual(values = alpha(spp.col, 0.5)) + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x = NULL, y = NULL) 
  
  p2=ggplot(spp.b.vd.df) + geom_label_repel(mapping = aes(x = Ins,y = MH, fill = class, label = sub(" Ca 1","",spp.short.names[-4])),fontface = 'italic') + 
    geom_point(mapping = aes(x = Ins,y = MH)) +
    scale_fill_manual(values = alpha(spp.col, 0.5)) + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x = NULL, y = NULL) 
  
  p3=ggplot(spp.b.dj.df) + geom_label_repel(mapping = aes(x = Ins,y = MH, fill = class, label = sub(" Ca 1","",spp.short.names[-4])),fontface = 'italic') + 
    geom_point(mapping = aes(x = Ins,y = MH)) +
    scale_fill_manual(values = alpha(spp.col, 0.5)) + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x = NULL, y = NULL) 
  
  #labs(x = "Mean insertions per molecule", y = "Mean microhomology residues per molecule")
  
  
}
pdf(paste0(spp.fig.folder,"Fig3.E.pdf"),width = 12,height = 4)
grid.arrange(grobs = list(p1,p2,p3),nrow= 1)
dev.off()
