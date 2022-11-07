# Figs A and B
cx.ax = 1.9
{
  pdf(file = paste0(spp.fig.folder,"Fig1.A_and_B.pdf"),width = 12,height = 4)
  par(mfrow=(c(1,2)) , mar = rep(4.5,4))
  bp.names = 25:55
  bp.names[bp.names%%6 !=0] = NA
  barplot(prop.table(table(factor(do.call(rbind,ZF.Ca.wt.df.list)$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.3), axes = F,cex.names=cx.ax)
  axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1,cex.axis=cx.ax)
  
  barplot(prop.table(table(factor(TCR.a.spp$MM$L,levels = 25:55))), col = c(out.of.frame.col,out.of.frame.col,in.frame.col), names = bp.names, ylim = c(0,0.3),axes = F,cex.names=cx.ax)
  axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1,cex.axis=cx.ax)
  barplot(prop.table(table(factor(do.call(rbind,ZF.Ca.mut.df.list)$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.30), axes = F,cex.names=cx.ax)
  axis(2,labels = as.character(c("0","0.10","0.20","0.30")),at = c(0,0.10,0.20,0.30),las = 1,cex.axis=cx.ax)
  barplot(prop.table(table(factor(MM.a.mut.a32.dfl.MUT$L,levels = 25:55))),col = c(out.of.frame.col,out.of.frame.col,in.frame.col),names = bp.names, ylim = c(0,0.12), axes =F,cex.names=cx.ax)
  axis(2,labels = as.character(c("0","0.04","0.08","0.12")),at = c(0,0.04,0.08,0.12),las = 1,cex.axis=cx.ax)
  
  dev.off()
}
# Fig C
{
ggsave(VDJ.plot.2(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF, J.start<=70 & J.start < L),],every.other.x.label = T,VJ.col = "black",alfa = 0.60),filename = paste0(spp.fig.folder,"/Fig1.C.ZF.jpg"),width = 6,height = 8,bg = 'white')
ggsave(VDJ.plot.2(TCR.a.spp.VJ$MM[with(TCR.a.spp.VJ$MM, J.start<=70 & J.start < L),],estetica = 2,every.other.x.label = T,VJ.col = "black",alfa = 0.60),filename = paste0(spp.fig.folder,"/Fig1.C.MM.jpg"),width = 6,height = 8,bg = 'white')

VDJ.plot.2(TCR.a.spp.VJ$ZF[with(TCR.a.spp.VJ$ZF, J.start<=70 & J.start < L),],every.other.x.label = T,VJ.col = "black",alfa = 0.60)
VDJ.plot.2(TCR.a.spp.VJ$MM[with(TCR.a.spp.VJ$MM, J.start<=70 & J.start < L),],estetica = 2,every.other.x.label = T,VJ.col = "black",alfa = 0.60)
}
# Fig D

{
  pdf(file = paste0(spp.fig.folder,"/Fig1.D.pdf"),width = 16,height = 10)
  par(mfrow= c(2,1), mar = rep(2.5,4))
  letras.by.6 = 1:42
  letras.by.6[letras.by.6%%6!=0] = ""
  do.call('entropy.output.VorJ',list( h.spp.a$ZF$nt[[42]],plot = 4,letras = letras.by.6,entropy.limit = 2.0,cex.letras = 1))
  for(l in seq(0,2,by = 0.5))abline(l,0,lty = 1,lwd = 0.3, col = "black")
  #for(i in seq(0,42,by = 6))text((i)*50.5/42-0.4,-0.2,i, xpd=TRUE,cex = 1.2)
  #for(i in seq(12,42,by = 6))text(i*(14.18/12),-0.2,i, xpd=TRUE,cex = 1.2)
  par(new = T)
  do.call('entropy.output.VorJ',list( h.spp.a$ZF$nt[[42]],plot = 4,letras = letras.by.6,entropy.limit = 2.0,cex.letras = 1))
  do.call('entropy.output.VorJ',list( h.spp.a$MM$nt[[42]],plot = 4,letras = letras.by.6,entropy.limit = 2.05))
  for(l in seq(0,2,by = 0.5))abline(l,0,lty = 1,lwd = 0.3, col = "black")  
  par(new = T)
  do.call('entropy.output.VorJ',list( h.spp.a$MM$nt[[42]],plot = 4,letras = letras.by.6,entropy.limit = 2.05))
  
  dev.off()
}


# Fig E

{
  
  pdf(file = paste0(spp.fig.folder,"/Fig1.E.pdf"),width = 5,height = 5)
  {
    plot(as.numeric(prop.table(table(rep(TCR.a.spp.VJ$ZF$pub,TCR.a.spp.VJ$ZF$Umi.count)))), ylim = c(0,0.6),type= 'l', ylab = 'Fraction of the repertoire', xlab = 'Number of individuals',las = 1,cex.axis = 1.2, cex.lab = 1.5)
    points(as.numeric(prop.table(table(rep(TCR.a.spp.VJ$ZF$pub,TCR.a.spp.VJ$ZF$Umi.count)))), pch =15)
    lines(as.numeric(prop.table(table(rep(TCR.b.spp.VJ$ZF$pub,TCR.b.spp.VJ$ZF$Umi.count)))), lty =2)
    points(as.numeric(prop.table(table(rep(TCR.b.spp.VJ$ZF$pub,TCR.b.spp.VJ$ZF$Umi.count)))), pch =2)
    legend(4, 58/100, legend=c("TCR alpha", "TCR beta"),
           pch=c(15, 2), lty=1:2, cex=1)
  }
  dev.off()  
  
}

