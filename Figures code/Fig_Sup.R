### Fig S1 ###

{
  pdf(paste0(spp.fig.folder,"Fig_S1_V_germline_position.pdf"),width = 10,height = 5)
  grid.arrange(grobs = lapply(1:2,function(i)ggplot(data=df.ZF.RSS.list[[i]], aes(x=RSS, y=Freq,fill = factor(RSS))) +
                                geom_bar(stat="identity", width=1,colour = "black",alpha = 2.3-i) +
                                scale_fill_manual(values=(c("black","black","black",J2.col,rep("black",4),J3.col,"black"))) +
                                scale_x_discrete(name="nucleotide position of V RSS") + 
                                scale_y_continuous(name=c("Fraction of usage in the repertoire","Number of germline V elements")[i]) +
                                theme_minimal() +
                                theme(legend.position = "none", 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      axis.text = element_text(size = 12),
                                      axis.title = element_text(size = 15)) 
  ),nrow =1 )
  dev.off()
  
  {
    p.V.RSS.1 = ggseqlogo(as.character(get.V.end.Seq("ZF","a",largo = 21)[spp.TCR.dict$ZF$a$V.GR.dict$heptamer.start==15])) +
      scale_y_continuous(limits = c(0,2)) +
      theme(axis.ticks.x = element_blank(),
            #axis.text.x=element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(color = "black",fill = NA,size = 0.3))
    p.V.RSS.2 = ggseqlogo(as.character(get.V.end.Seq("ZF","a",largo = 21)[spp.TCR.dict$ZF$a$V.GR.dict$heptamer.start==15]),col_scheme = cs.black)+
      scale_y_continuous(limits = c(0,2)) +
      theme(axis.ticks.x = element_blank(),
            #axis.text.x=element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(color = "black",fill = NA,size = 0.3))
    p.V.20.RSS.1 = ggseqlogo(as.character(get.V.end.Seq("ZF","a",largo = 26)[spp.TCR.dict$ZF$a$V.GR.dict$heptamer.start==20]))+
      scale_y_continuous(limits = c(0,2)) +
      theme(axis.ticks.x = element_blank(),
            # axis.text.x=element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(color = "black",fill = NA,size = 0.3))
    p.V.20.RSS.2 = ggseqlogo(as.character(get.V.end.Seq("ZF","a",largo = 26)[spp.TCR.dict$ZF$a$V.GR.dict$heptamer.start==20]),col_scheme = cs.black)+
      scale_y_continuous(limits = c(0,2)) +
      theme(axis.ticks.x = element_blank(),
            #axis.text.x=element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_rect(color = "black",fill = NA,size = 0.3))
    pdf(paste0(spp.fig.folder,"Fig_S1_V.15.RSS.pdf"),width = 21/3,height = 7/3)
    grid.arrange(grobs = list(p.V.RSS.1,p.V.RSS.2),ncol = 1)
    dev.off()
    
    pdf(paste0(spp.fig.folder,"Fig_S1_V.20.RSS.pdf"),width = 26/3,height = 7/3)
    grid.arrange(grobs = list(p.V.20.RSS.1,p.V.20.RSS.2),ncol = 1)
    dev.off()
  }
}

### Fig S2 ###
sp ="ZF"
gene = "a"
pdf(file =paste0(spp.fig.folder,"Fig_S2_consensus.pdf"),width = 8,height = 3)
ggseqlogo(get.J.RSS.matrix("ZF","a"))
dev.off()
pdf(file =paste0(spp.fig.folder,"Fig_S2_score.pdf"),width = 12,height = 6)
plot(RSS.si.ZF.a$score.RSS[,25])
dev.off()
pdf(file =paste0(spp.fig.folder,"Fig_S2_sequence.pdf"),width = 80,height = 3)
ggseqlogo(as.character(get.J.Seq("ZF","a")[25]))
dev.off()
RSS.si.ZF.a = RSS.scorer.iterative(get.J.RSS.matrix("ZF","a"),getSeq(spp.genome.data.list[[sp]][[gene]],shift(resize(spp.TCR.dict[[sp]][[gene]]$J.GR.dict,width = 214+27,fix = "end"),14)),check.to = 200,iter.max = 5)
pdf(file = paste0(spp.fig.folder,"Fig_S2_score_matrix.pdf"),width = 8,height = 8)
image(RSS.si.ZF.a$score.RSS[,order(RSS.si.ZF.a$RSS.index)])
dev.off()
pdf(file = paste0(spp.fig.folder,"Fig_S2_score_matrix.ref.pdf"),width = 8,height = 3)
image(matrix(min(RSS.si.ZF.a$score.RSS):max(RSS.si.ZF.a$score.RSS)))
dev.off()

### Fig S3 ###

{  
  
  # Fig J RSS  
  pdf(paste0(spp.fig.folder,"Fig_S3_J_germline_position.pdf"),width = 10,height = 5)
  grid.arrange(grobs = lapply(3:4,function(i)ggplot(data=df.ZF.RSS.list[[i]], aes(x=RSS, y=Freq,fill = factor(RSS))) +
                                scale_x_discrete(breaks = seq(-57,-69,by = -3)) +
                                geom_bar(stat="identity", width=1,colour = "black",alpha = 4.3-i) +
                                scale_fill_manual(values=(c(J1.col,"white","white",J2.col,"white","white",J3.col,"black","black",J4.col,"black","white",J5.col)))+
                                labs (y =c("Fraction of usage in the repertoire","Number of germline V elements")[i-2] , x ="nucleotide position of V RSS") + 
                                
                                theme_minimal() +
                                theme(legend.position = "none", 
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      axis.text = element_text(size = 12),
                                      axis.title = element_text(size = 15))
                              
  ),nrow = 1)
  dev.off()
  
  lapply(seq(-57,-69,by = -3)[1],function(x)subseq(as.character(get.J.Seq("ZF","a",100)[spp.TCR.dict$ZF$a$J.GR.dict$heptamer.end==x]),28))
  
  cM.Js   = lapply(seq(-57,-69,by = -3),function(x)consensusMatrix(subseq(as.character(get.J.Seq("ZF","a",100)[spp.TCR.dict$ZF$a$J.GR.dict$heptamer.end==x]),45+57-7+x)))
  cM.Js.2 = cM.Js
  for(i in 1:4)cM.Js.2[[i]] = cbind(matrix(rep(25,12*(5-i)),ncol = 3*(5-i)),cM.Js.2[[i]])
  
  {
    
    p.J.RSS.1 = lapply(1:5,function(i)ggseqlogo(cM.Js.2[[i]][DNA_BASES,]) +
                         scale_y_continuous(limits = c(0,2)) +
                         theme(axis.ticks.x = element_blank(),
                               #axis.text.x=element_blank(), 
                               axis.title.y = element_blank(),
                               axis.text.y=element_blank(),
                               panel.border = element_rect(color = "black",fill = NA,size = 0.3)))
    
    p.J.RSS.2 = lapply(1:5,function(i)ggseqlogo(cM.Js.2[[i]][DNA_BASES,],col_scheme = cs.black) +
                         scale_y_continuous(limits = c(0,2)) +
                         theme(axis.ticks.x = element_blank(),
                               #axis.text.x=element_blank(), 
                               axis.title.y = element_blank(),
                               axis.text.y=element_blank(),
                               panel.border = element_rect(color = "black",fill = NA,size = 0.3)))
    
    
    
    
    
    
    pdf(paste0(spp.fig.folder,"Fig_S3_Js.RSS.pdf"),width = 18,height = 8)
    grid.arrange(grobs = p.J.RSS.1,ncol = 1)
    dev.off()
    
    pdf(paste0(spp.fig.folder,"Fig_S3_Js.black.RSS.pdf"),width = 18,height = 8)
    grid.arrange(grobs = p.J.RSS.2,ncol = 1)
    dev.off()
    
    
  }
}  

### Fig S4 ###

{
pdf(file =  paste0(spp.fig.folder,"Fig_S4_rank_pub.pdf"),width = 8,height = 8)
with(TCR.a.spp.VJ$ZF[TCR.a.spp.VJ$ZF$Sample.bio.name==2,],plot(sort(Umi.proportion,decreasing = T),log = 'xy',type = 'l',xlab = NA,ylab = NA,las = 1))
sapply(1:6,function(i)with(TCR.a.spp.VJ$ZF[TCR.a.spp.VJ$ZF$Sample.bio.name==i,],lines(sort(Umi.proportion,decreasing = T),col = i)))
legend(10000,0.1,legend = 1:6,fill = 1:6)
dev.off()

{
  sp = "ZF"
  df= TCR.a.spp.VJ[[sp]]
  df.sel = df[,with(df,c("Umi.count","pub","Sample.bio.name"))]
  temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
  temp.agg$cs = unname(unlist(tapply(temp.agg$Umi.count,temp.agg$Sample.bio.name,cumsum)))
  temp.agg$cs.2 = unlist(lapply(tapply(temp.agg$Umi.count,temp.agg$Sample.bio.name,cumsum),function(x)(c(0,x)[-7]+x)/2))
  p.pub.s4 = ggplot(temp.agg) + geom_col(mapping = aes(x = Sample.bio.name,y = Umi.count, fill = pub<2), colour = "black") +
    geom_text(mapping = aes(x = Sample.bio.name, y = cs.2 , label = rep(1:6,6)),size = 7) +
    labs(y = "Umi molecules", x = "Specimen") +
    theme_minimal()
  ggsave(paste0(spp.fig.folder,"Fig_S4_pub.pdf"),p.pub.s4,width = 6,height = 9)
  
  sp = "ZF"
  df= TCR.a.spp.VJ[[sp]]
  df.sel = df[,with(df,c("pc.nt","pub","Sample.bio.name"))]
  temp.agg = aggregate(pc.nt ~ ., data = df.sel, FUN = length)
  temp.agg$cs = unname(unlist(tapply(temp.agg$pc.nt,temp.agg$Sample.bio.name,cumsum)))
  temp.agg$cs.2 = unlist(lapply(tapply(temp.agg$pc.nt,temp.agg$Sample.bio.name,cumsum),function(x)(c(0,x)[-7]+x)/2))
  p.pub.s4.2 = ggplot(temp.agg) + geom_col(mapping = aes(x = Sample.bio.name,y = pc.nt, fill = pub<2), colour = "black") +
    geom_text(mapping = aes(x = Sample.bio.name, y = cs.2 , label = rep(1:6,6)),size = 7) +
    labs(y = "Clonotypes", x = "Specimen") +
    theme_minimal()
  ggsave(paste0(spp.fig.folder,"Fig_S4_pub.2.pdf"),p.pub.s4.2,width = 6,height = 9)
}
}

### Fig S5 ###

 # Use Fig 2B
pdf(file =  paste0(spp.fig.folder,"Fig_S5_consensus.pdf"),width = 8,height = 4)

  grid.arrange(grobs = lapply(1:4,function(i)ggseqlogo(with(TCR.a.spp.VJ$ZF,CDR3.nucleotide.sequence[rules==i & L==36+3*i-3*sum(grepl(i,c(3,4)))]),col_scheme = cs.black) +
                                scale_y_continuous(limits = c(0,2)) +
                                theme(axis.ticks.x = element_blank(),
                                      axis.text.x=element_blank(), 
                                      axis.title.y = element_blank(),
                                      axis.text.y=element_blank(),
                                )),ncol = 1)
  

  
  dev.off()
  
  
  ### Fig S6 ###
  
  {  
    {
      sp = "ZF"
      df= TCR.a.spp[[sp]][with(TCR.a.spp[[sp]], L>=27 & L<=57 & RSS.J %in% seq(-57,-69,by=-3)),]
      df.sel = df[,with(df,c("Umi.count","RSS.V","RSS.J","L"))]
      temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
      temp.agg$RSS.J = factor(temp.agg$RSS.J,levels = c(seq(-57,-69,by = -3)))
    }
    p.J.RSS.cat.sup = ggplot(data=temp.agg, aes(x=L, y=Umi.count/sum(Umi.count),fill = factor(RSS.J))) +
      geom_bar(stat="identity", width=1,colour = "black") +
      scale_fill_manual(values=(c(J1.col,J2.col,J3.col,J4.col,J5.col)))+
      scale_x_continuous(name="CDR3 length (nucleotides)",breaks = seq(27,52, by = 3)) + 
      
      labs (y = "Fraction of the repertoire" , legend = "J RSS position") +
      theme_minimal()  
    ggsave(paste0(spp.fig.folder,"Fig_S6_J.RSS.cat.sup.pdf"),p.J.RSS.cat.sup,width = 6,height = 4,bg = "white")
    
    p.J.RSS.cat.sep = lapply(1:5,function(x)ggplot(data=temp.agg[temp.agg$RSS.J %in%seq(-57,-69,by=-3)[x],], aes(x=factor(L,levels = 25:53), y=Umi.count/sum(Umi.count),fill = factor(RSS.J))) +
                               geom_bar(stat="identity", width=1,colour = "black") +
                               scale_x_discrete(name=element_blank(),breaks = seq(0, 50, by = 3),limits = factor(27:52)) + 
                               scale_y_continuous(name=element_blank(),breaks = seq(0, 0.6, by = 0.2),limits = c(0,0.75)) + 
                               
                               scale_fill_manual(values=c(J1.col,J2.col,J3.col,J4.col,J5.col)[x])+
                               labs (legend = element_blank()) +
                               theme_minimal() +
                               theme(panel.grid.minor = element_blank()))
    pdf(paste0(spp.fig.folder,"Fig_S6_J.RSS.cat.sep.pdf"),width = 9,height = 8)
    grid.arrange(grobs = p.J.RSS.cat.sep,ncol = 1)
    dev.off()
  }
### Fig S7 ###
    rules.prop = prop.table(table(with(TCR.a.spp.VJ[[sp]],rep(rules,Umi.count))))[2:5]
    rules.VJ = lapply(1:4,function(rule)with(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$rules==rule,],sapply(as.numeric(prop.table(table(rep(V.gene,Umi.count)))),function(x)sapply(as.numeric(prop.table(table(rep(J.gene,Umi.count)))),function(y)x*y))))
    rules.VJ.prop = unlist(mapply('*',rules.VJ,rules.prop))
    
    {
      pdf(file = paste0(spp.fig.folder,"Fig_S7_rules.prob.pdf"),width = 5,height = 5)
      plot(pgeom(300000,sort(rules.VJ.prop,decreasing = T)),type = 'l',ylab = "Probability of generation",xlab = "Rank",xaxt = "none",las = 2 )
      axis(side = 1,seq(0,50000,by = 10000),tick = T)
      grid(nx = NULL, ny = NULL,
           lty = 1,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)      # Grid line width
      
      lines(pgeom(300000,sort(rules.VJ.prop,decreasing = T)),col = 'black')
      lines(pgeom(200000,sort(rules.VJ.prop,decreasing = T)),col = 'red')
      dev.off()
    }
  
### Fig S9 ###
    {

    p1s9=VDJ.plot.2(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,V.gene == "TRAV11_6_01" & J.start < L),],VJ.col = "black",alfa = 0.65) + ggtitle("TRAV11_6_01")
    p2s9=VDJ.plot.2(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,V.gene == "TRAV52_1_01" & J.start < L),],VJ.col = "black",alfa = 0.65) + ggtitle("TRAV52_1_01")
    p3s9=VDJ.plot.2(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,J.gene == "TRAJ5_01" & J.start < L),],VJ.col = "black",alfa = 0.65) + ggtitle("TRAJ5_01")
    ggsave(paste0(spp.fig.folder,"Fig_S9_TRAV11_6_01.pdf"),p1s9,width = 10,height = 8)
    ggsave(paste0(spp.fig.folder,"Fig_S9_TRAV52_1_01.pdf"),p2s9,width = 10,height = 8)
    ggsave(paste0(spp.fig.folder,"Fig_S9_TRAJ5_011.pdf"),p3s9,width = 10,height = 8)
    }
    {
    pdf(paste0(spp.fig.folder,"Fig_S9_barplots.pdf"),width = 5,height = 10)
    par(mfrow = c(3,1))
    bp.names = 25:55
    bp.names[bp.names%%3 !=0] = NA
    barplot(prop.table(table(factor(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,J.gene == "TRAJ5_01"),]$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.3), axes = F,main = "TRAJ5_01")
    axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1)
    
    bp.names = 25:55
    bp.names[bp.names%%3 !=0] = NA
    barplot(prop.table(table(factor(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,V.gene == "TRAV11_6_01"),]$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.3), axes = F,main = "TRAV11_6_01")
    axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1)
    
    bp.names = 25:55
    bp.names[bp.names%%3 !=0] = NA
    barplot(prop.table(table(factor(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF,V.gene == "TRAV52_1_01"),]$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.3), axes = F,main = "TRAV52_1_01")
    axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1)
    
    dev.off()    
    }
    
### Fig S10 ###
    
    {  
      # In frame Vs
      
      df.V.frame.df.2 = df.V.frame.df
      df.V.frame.df.2$Vs[!df.V.frame.df.2$Vs%in% c("TRAV11_6_01","TRAV52_1_01")] = ""
      
      p.if.Vs = ggplot(df.V.frame.df.2) +  
        geom_boxplot(mapping = aes(x = factor(RSS),group = factor(RSS),y = in.frame, fill = factor(RSS))) + 
        geom_jitter(mapping = aes(x = factor(RSS),y = in.frame), size = 0.5) +
        geom_text_repel(mapping = aes(x = factor(RSS),y = in.frame,label = Vs)) +
        scale_y_continuous(limits = c(0,1)) +
        theme_linedraw() +
        theme(legend.position="none") +
        labs(x = "V RSS position", y = "Fraction of in-frame rearrangements by V element" ) +
        scale_fill_manual(values=(c(J1.col,J2.col,J3.col,J3.col,J4.col,J5.col))) 
      
      ggsave(paste0(spp.fig.folder,"Fig_S10_Vs_in_frame.pdf"),p.if.Vs,width = 8,height = 6)
      # In frame Js
      
      df.J.frame.df.2 = df.J.frame.df
      df.J.frame.df.2$Js[df.J.frame.df.2$Js !="TRAJ5_01"] = ""
      
      df.J.frame.df.2$RSS = factor(df.J.frame.df.2$RSS,levels = sort(unique(factor(df.J.frame.df.2$RSS)),decreasing = T))
      
      p.if.Js = ggplot(df.J.frame.df.2) +  
        geom_jitter(mapping = aes(x = RSS,y = in.frame), size = 0.5) + #, color = factor(RSS))) + 
        geom_text_repel(mapping = aes(x = RSS,y = in.frame,label = Js)) +
        geom_boxplot(mapping = aes(x = RSS,group = RSS,y = in.frame, fill = RSS)) + 
        scale_y_continuous(limits = c(0,1)) +
        theme_linedraw() +
        labs(x = "J RSS position", y = "Fraction of in-frame rearrangements by J element" ,fill = "J RSS position") +
        scale_fill_manual(values=(c(J1.col,J2.col,J3.col,"black","black",J4.col,J5.col))) +
        theme(legend.position = 'none')
      
      ggsave(paste0(spp.fig.folder,"Fig_S10_Js_in_frame.pdf"),p.if.Js,width = 8,height = 6)
    }
    
    
    ### Fig S11 ###
    
    {  
      # Connor fish barplots
      target_seq_14 = "TGTGCTCTGAGGCC"
      target_seq_15 = "TGTGCT..TGAGGCC"
     
      mut.L.table = table(factor(ZF.23.25.dla.a$L[grepl(target_seq_15,ZF.23.25.dla.a$CDR3.nucleotide.sequence)],levels = 30:55))
      control.L.table = table(factor(ZF.23.25.dla.a$L[grepl(target_seq_14,ZF.23.25.dla.a$CDR3.nucleotide.sequence)],levels = 30:55))
      
      pdf(paste0(spp.fig.folder,"Fig_S11.CF.barplots.wt.pdf"),width = 5,height = 5)
      barplot(prop.table(control.L.table),col = c(in.frame.col,"grey","grey"),yaxt = "n")
      axis(2,seq(0,0.25,by =0.05),las = 2)
      dev.off()
      pdf(paste0(spp.fig.folder,"Fig_S11.CF.barplots.mut.pdf"),width = 5,height = 5)
      barplot(prop.table(mut.L.table),col = c(in.frame.col,"grey","grey"),yaxt = "n")
      axis(2,seq(0,25,by =0.05),las = 2)
      dev.off()
      
      pdf(paste0(spp.fig.folder,"Fig_S11.CF.barplot.if.pdf"),width = 2.6,height = 5)
      
      barplot(matrix(c(prop.table(rev(table(ZF.23.25.dla.a$L[grepl(target_seq_14,ZF.23.25.dla.a$CDR3.nucleotide.sequence)]%%3==0))),prop.table(rev(table(ZF.23.25.dla.a$L[grepl(target_seq_15,ZF.23.25.dla.a$CDR3.nucleotide.sequence)]%%3==0)))),ncol = 2),col = c(in.frame.col,"grey"),las = 2)
      dev.off()
      # Connor fish consensus
      
      p.CF.14 = ggseqlogo(subseq(ZF.23.25.dla.a[grepl(target_seq_14,ZF.23.25.dla.a$CDR3.nucleotide.sequence),]$CDR3.nucleotide.sequence,1,14),col_scheme = cs.black) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
      
      p.CF.15 = ggseqlogo(subseq(ZF.23.25.dla.a[grepl(target_seq_15,ZF.23.25.dla.a$CDR3.nucleotide.sequence),]$CDR3.nucleotide.sequence,1,15),col_scheme = cs.black) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
      
      ggsave(filename = paste0(spp.fig.folder,"Fig_S11.CF.14.pdf"),p.CF.14,width = 16,height = 2)  
      ggsave(filename = paste0(spp.fig.folder,"Fig_S11.CF.15.pdf"),p.CF.15,width = 16,height = 2)  
    }
    
    
    {  
      
      # VDJ plot with rules VJ new color
      rules.plot.list = list()
      MH.col = "cyan"
      for(sp in names(TCR.a.spp)[c(6,10,7,8,1,2,3,9,11,4,5)]){
        {
          df = TCR.a.spp.VJ[[sp]][!is.na(with(TCR.a.spp.VJ[[sp]],RSS.V+RSS.J & L <= 60)),]
          df$D5.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0] + 1
          df$D3.end[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0] - 1
          df$V.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0]
          df$J.start[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0]
          
          rules.col = c('white','sienna4','magenta','dark green','darkorange 1')
          J.RSS.col = c(J1.col,J2.col,J3.col,J4.col,J5.col)
          
          
          Umi.lim = 1
          plot.D = F
          titulo = ''
          estetica = 5
          df = df[df$L %% 3 ==0,]
          df$rules[!df$rules %in% 1:4] = "r0"
          df$rules[df$rules %in% 1:4] = paste0("r",  df$rules[df$rules %in% 1:4])
          df$rules[df$rules == "r0" & df$Total.insertions < 0] = "r5"
          df$rules = factor(df$rules,levels = c("r1","r2","r3","r4","r5","r0"))
          df.sel = df[,with(df,c("Umi.count","V.end","J.start","L","D5.end","D3.end",'rules','Total.insertions'))]
          estetica.vj = F
          temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
          
          temp.agg2 = temp.agg[temp.agg$Umi.count>=Umi.lim,]
          temp.agg3 = temp.agg2[rev(with(temp.agg2,order(rules,V.end+J.start,L))),]
          
          
          v.x1=rep(0,nrow(temp.agg3))
          v.x2=temp.agg3$V.end
          j.x1=(temp.agg3$J.start)-1
          j.x2=(temp.agg3$L)
          vj.filter = temp.agg3$V.end>=temp.agg3$J.start
          vj.x1=temp.agg3$J.start[vj.filter]-1
          vj.x2=temp.agg3$V.end[vj.filter]
          nn.filter = temp.agg3$V.end<temp.agg3$J.start-1
          nn.x1 = temp.agg3$V.end[nn.filter]
          nn.x2 = temp.agg3$J.start[nn.filter]-1
          vj.x1 = rep(-1,length(temp.agg3$D5.end))
          vj.x2 = rep(-1,length(temp.agg3$D5.end))
          vj.x1 = temp.agg3$D5.end-1
          vj.x2 = temp.agg3$D3.end
          vj.filter = vj.x1 >0
          y1=c(0,cumsum(100*prop.table(temp.agg3$Umi.count))[-nrow(temp.agg3)])
          y2=cumsum(100*prop.table(temp.agg3$Umi.count))
          
          rules.x1 = rep(60,nrow(temp.agg3))
          
          rules.x2 = rep(65,nrow(temp.agg3))
          rule = temp.agg3$rules
          d1=data.frame(x1=c(v.x1,j.x1),x2=c(v.x2,j.x2),y1=c(y1,y1),y2=c(y2,y2), t=factor(c(rep(c('V','J'),each=length(v.x1))),levels = c("V","J")))
          d2=rbind(d1,data.frame(x1=nn.x1,x2=nn.x2,y1=y1[nn.filter],y2=y2[nn.filter], t=factor(rep("N",sum(nn.filter)),levels = c("V","J","N"))))
          
          d3 = rbind(d2,data.frame(x1 = vj.x1[vj.filter],x2 = vj.x2[vj.filter],y1=y1[vj.filter],y2=y2[vj.filter],t=factor(rep("VJ",sum(vj.filter)),levels = c("V","J","N","VJ"))))
          colores = c(V.col,J.col,N.col,VJ.col,rules.col,MH.col)
          if (estetica %in% 1:5)d4 = with(temp.agg3,data.frame(L = unique(L), Umi.count = rev(tapply(Umi.count,L,sum))))
          if (estetica %in% 5)d4 = with(temp.agg3,data.frame(L = unique(rules), Umi.count = rev(tapply(Umi.count,rules,sum))))
          d4$y1 = c(0,cumsum(100*prop.table(d4$Umi.count))[-nrow(d4)])
          d4$y2 = c(cumsum(100*prop.table(d4$Umi.count)))
          d5 = rbind(d3,data.frame(x1 = rules.x1,x2 = rules.x2,y1 = y1,y2 = y2,t = factor(rule,levels = c("V","J","N","r0","r1","r2","r3","r4","r5"))))
          labels.x = seq(0, 60, by = 3)
          labels.x[c(F,T)] = ""
          p1= ggplot() + 
            scale_x_continuous(name=element_blank(),breaks = seq(0, 60, by = 3),limits = c(0,80),labels = labels.x) + 
            scale_y_continuous(name=element_blank()) +
            geom_rect(data=d5, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=NA,size = 0.1, alpha=0.7) +
            geom_rect(data=d4, mapping=aes(xmin=0, xmax=65, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.4) +
            geom_rect(data=d4, mapping=aes(xmin=60, xmax=65, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.4, alpha = 1) +
            scale_fill_manual(values=colores) +
            theme_minimal(base_family = "sans") +
            ggtitle(label = (spp.short.names[match(sp,names(TCR.a.spp)[c(6,10,7,8,1,2,3,9,11,4,5)])])) +
            
            theme(panel.grid.major.x = element_line(color = "grey",size = 0.3,linetype = 2),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(size = 12),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank() , 
                  legend.position = 'none',
                  plot.title = element_text(size = 20, face = "italic"))
        }
        
        rules.plot.list[[sp]] = p1
      }
      
      p.spp = grid.arrange(grobs = rules.plot.list,ncol = 3)
      
      p.spp.1 = grid.arrange(grobs = rules.plot.list[1:6],ncol = 3)
      p.spp.2 = grid.arrange(grobs = rules.plot.list[7:11],ncol = 3)
      
      ggsave(filename = paste0(spp.fig.folder,"Fig_S12.pdf"),p.spp,width = 18,height = 25)
      
      
    }

    
### Fig S13 ###
    
    
{ 
    {
      pdf(file = paste0(spp.fig.folder,"Fig_S13.A.pdf"),width = 8,height = 18)
      par(mfrow = c(11,1),mar = c(2,0,2,0))
      for(i in c(6,10,7,8,1,2,3,9,11,4,5)){do.call('entropy.output.VorJ',list( h.spp.a[[i]]$nt[[42]],plot = 4))
        for(l in seq(0,2,by = 0.5))abline(l,0,lty = 1,lwd = 0.3, col = "black")  
        par(new = T)
        do.call('entropy.output.VorJ',list( h.spp.a[[i]]$nt[[42]],plot = 4))}
      dev.off()  
      
      spp.repertoire.df.figS = spp.repertoire.a.df
      
      spp.repertoire.df.figS$sp = factor(c("PP","DR","OM","MM","LA","CP","AR.1","AR.2","OL","PS","PA")[c(6,10,7,8,1,2,3,9,11,4,5)], levels = c("PP","DR","OM","MM","LA","CP","AR.1","AR.2","OL","PS","PA")[c(6,10,7,8,1,2,3,9,11,4,5)])
      spp.repertoire.df.figS$class = factor(spp.repertoire.df.figS$class ,levels = unique(spp.repertoire.df.figS$class))
      levels(spp.repertoire.df.figS$class)[levels(spp.repertoire.df.figS$class) == "transition"] = "lobe-finned fishes"
      
      p.S13.ins =ggplot(spp.repertoire.df.figS) + 
        geom_col(mapping = aes(x = sp, y = as.numeric(No.insertions),fill = class)) + 
        scale_fill_manual(values=c('darkseagreen 2','cadetblue2',"royal blue",'dark orchid',"red"),name = "group") +
        labs(y = "Fraction of the repertoire without non-template nucleotides", x = element_blank()) +
        theme_classic() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12))
      
      ggsave(file = paste0(spp.fig.folder,"Fig_S13.B.pdf"),p.S13.ins,width = 8,height = 5)
      
      p.S13.H = ggplot(spp.repertoire.df.figS) + 
        geom_col(mapping = aes(x = sp, y = as.numeric(H.N),fill = class)) + 
        #scale_color_manual(values=c('cadetblue2',"red",'darkseagreen 2',"royal blue",'dark orchid')) +
        scale_fill_manual(values=c('darkseagreen 2','cadetblue2',"royal blue",'dark orchid',"red"),name = "group") +
        theme_classic() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12))
      
      
      ggsave(file = paste0(spp.fig.folder,"Fig_S13.C.pdf"),p.S13.H,width = 8,height = 5)
      
    } 
    
    MH.nt.tb = sapply(1:11,function(sp)sum(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$Total.insertions < 0 ,]$Umi.count))
    blunt.nt.tb = sapply(1:11,function(sp)sum(TCR.a.spp.VJ[[sp]][TCR.a.spp.VJ[[sp]]$Total.insertions == 0 ,]$Umi.count))
    
    Pv.nt.pc.nt = lapply(1:11,function(sp)lapply(1:5,function(l)with(TCR.a.spp.VJ[[sp]][with(TCR.a.spp.VJ[[sp]],V.end < L & V.end > l & Total.insertions == l),],pc.nt[(V.end == RSS.V-1) & subseq(CDR3.nucleotide.sequence,end = V.end,width = l)==reverseComplement(DNAStringSet(CDR3.nucleotide.sequence,start = V.end+1,width = l))])))
    Pj.nt.pc.nt = lapply(1:11,function(sp)lapply(1:5,function(l)with(TCR.a.spp.VJ[[sp]][with(TCR.a.spp.VJ[[sp]],J.start < L & J.length + RSS.J == -34 & Total.insertions == l),],pc.nt[subseq(CDR3.nucleotide.sequence,end = J.start-1,width = l)==reverseComplement(DNAStringSet(CDR3.nucleotide.sequence,start = J.start,width = l))])))
    
    Pvj.nt.tb = sapply(1:11,function(sp)sapply(1:5,function(l)sum(TCR.a.spp.VJ[[sp]]$Umi.count[TCR.a.spp.VJ[[sp]]$pc.nt %in% unique(c(Pv.nt.pc.nt[[sp]][[l]],Pj.nt.pc.nt[[sp]][[l]]))])))
    Total.nt.tb = sapply(1:11,function(sp)sum(TCR.a.spp.VJ[[sp]]$Umi.count))
    
    {
      No.TdT.tb = rbind(MH.nt.tb,blunt.nt.tb,Pvj.nt.tb)
      rownames(No.TdT.tb) = c("MH","blunt",paste0("P",1:5))
      colnames(No.TdT.tb) = names(TCR.a.spp)
      No.TdT.tb = No.TdT.tb[,c(6,10,7,8,1,2,3,9,11,4,5)]
      Total.nt.tb = Total.nt.tb[c(6,10,7,8,1,2,3,9,11,4,5)]
      No.TdT.tb.2 = apply(No.TdT.tb,1,function(x)x/Total.nt.tb)
      No.TdT.tb.2 = cbind(No.TdT.tb.2,rowSums(No.TdT.tb.2))
      colnames(No.TdT.tb.2)[8] = "Total no TdT"
    }
    
    No.TdT.tb.3 = cbind(No.TdT.tb.2[,1:2],rowSums(No.TdT.tb.2[,3:7]))
    colnames(No.TdT.tb.3)[3] = "P nt"
    pdf(file = paste0(spp.fig.folder,"Fig_S13.D.pdf"),width = 18,height = 10)
    barplot(t(cbind(No.TdT.tb.3,1-rowSums(No.TdT.tb.3))),col = c(1,7,3,'white'))
    dev.off()    
}    
    
### Fig S14 ###
    
    {
      p.J.n_vs_range = ggplot(spp.df[!is.na(spp.df$J_frame),]) + geom_point(mapping = aes(y = J_range,x=J_nr,color = group),size = 0.8) +
        scale_color_manual(values=spp.col.xt) +
        theme_classic() + 
        labs (y= "Range between first and last found J genes (nt)", x = "Number of J elements") +
        scale_y_continuous(limits = c(0,150000))
      ggsave(paste0(spp.fig.folder,"Fig_S14.pdf"),p.J.n_vs_range,width = 6,height = 4)
    }

### Fig S15 ###
    
    
    p.J.ZF_vs_J.MM = ggplot(spp.df) + 
      geom_point(mapping = aes(x = J_frame, y = J_frame.MM)) +
      theme_classic() +
      labs (x = "Initial consensus query D. rerio", y = "Initial consensus query M. musculus")
    ggsave(paste0(spp.fig.folder,"Fig_S15.A.pdf"),p.J.ZF_vs_J.MM,width = 8,height = 8)
    
    {
      p.J.ec_vs_frame = ggplot(spp.df) + 
        geom_point(mapping = aes(x = J_frame,y=J_ec,color = group), size = 0.8) +
        geom_text_repel(mapping = aes(x = J_frame,y=J_ec,color = group,label = sapply(strsplit(spp.df$spp," "),function(x)paste(paste0(subseq(x[[1]],1,1),"."),x[[2]]))), size = 1.5,fontface = 'italic') +
        scale_color_manual(values=spp.col.xt) +
        labs (x= "Conservation of J RSS reading frame", y = "Entropy of J first 5 nucleotides") +

        theme_classic()
    }
    ggsave(paste0(spp.fig.folder,"Fig_S15.B.pdf"),p.J.ec_vs_frame,width = 10,height = 8)
    
    
### Fig S16 ###
    
    spp.short.names = c("C. punctatum","P. senegalus","A. ruthenus Ca 1","A. ruthenus Ca 2","P. progenetica","D. rerio","O. mykiss","O. latipes","P. annectens","M. musculus","L. africana")
    pdf(paste0(spp.fig.folder,"Fig_S16.ref.colors.pdf"),width = 10,height = 5)
    barplot(rep(1,20), col = c(rainbow(47)[c(9:13)],rainbow(20)[c(7:20,1)]),space = 0,names.arg = paste((1:20)*5,"%"))
    dev.off()
    pruned.colors = c(rainbow(47)[c(9:13)],rainbow(20)[c(7:20,1)])[round(as.numeric(spp.df$J_frame)*20)[match(pruned.tree.tip.translate,spp.df$spp)]]
    pruned.tree.tip.translate.2 = paste(paste0("(", round(as.numeric(spp.df$J_frame),2)[match(pruned.tree.tip.translate,spp.df$spp)],")"), sub("[a-z].* ",". ",pruned.tree.tip.translate) , sep =" ")
    pruned.tree.tip.translate.2[76:223] = paste(sub("[a-z].* ",". ",pruned.tree.tip.translate)[76:223],paste0("(", round(as.numeric(spp.df$J_frame),2)[match(pruned.tree.tip.translate,spp.df$spp)][76:223],")") , sep =" ")
    pruned_tree.2 = pruned_tree
    pruned_tree.2$tip.label = pruned.tree.tip.translate.2
  
    pdf(paste0(spp.fig.folder,"Fig_S16.pdf"),width = 9,height = 9)
    par(mfrow = c(1,1))
    plot(pruned_tree, type = 'fan',cex = 0.5, label.offset = .01, no.margin = TRUE, tip.color = pruned.colors, edge.color = col.pt)
    dev.off()