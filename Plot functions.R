library(ggplot2)
library(ggseqlogo)
library(gridExtra)
library(ggrepel)
library(msa)

V.col = 'royalblue'
J.col = 'red'
N.col = 'grey'
#VJ.col = 'darkorchid4'
VJ.col = 'black'
RSS.col = 'yellow1'
P.nt.col = "chartreuse1"
J1.col = 'darkgoldenrod2'
J2.col = 'darkorange1'
J3.col = 'brown2'
J4.col = 'deeppink2'
J5.col = 'darkorchid3'
D.col = 'dark green'

J1.col = 'yellow1'
J2.col = 'orange1'
J3.col = 'firebrick2'
J4.col = 'deeppink1'
J5.col = 'darkorchid3'

in.frame.col = 'dark green'
out.of.frame.col = 'grey'
J.RSS.col = c(J1.col,J2.col,J3.col,J4.col,J5.col)
rules.col = c('white','sienna4','magenta','dark green','cyan')

cs.V = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(V.col,4))
cs.J = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(J.col,4))
cs.VJ = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=as.character(rep(VJ.col,4)))
cs.RSS = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(RSS.col,4))
cs.black = make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep("black",4))

bf.col = 'cadetblue2'
mm.col = "red"
sh.col = 'darkseagreen 2'
tl.col = "royal blue"
tr.col = 'dark orchid'
spp.col.xt = c('darkseagreen 2','cadetblue2',"royal blue",'dark orchid','chocolate1','coral1','brown3','red')
spp.col = c(sh.col,bf.col,tl.col,tr.col,mm.col)
cs.spp = lapply(spp.col.xt,function(x)make_col_scheme(chars=c('A', 'T', 'C', 'G'), cols=rep(x,4)))


VDJ.plot = function(df,estetica = 1,estetica.vj = F,Umi.lim =1,lim = 60,alfa = 0.55 ,titulo ='',V.col = 'royalblue',J.col ='red',N.col ='grey',D.col = "dark green",VJ.col = "purple",plot.D = T,every.other.x.label = T){
  {
    df.sel = df[,with(df,c("Umi.count","V.end","J.start","L","D5.end","D3.end"))]
    temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
    
    temp.agg2 = temp.agg[temp.agg$Umi.count>=Umi.lim,]
    if(estetica ==1) temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,J.start,V.end,J.start-V.end))),]
    if(estetica ==2)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,J.start+V.end))),]
    if(estetica ==3)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,(J.start>V.end),J.start+V.end))),]
    if(estetica ==4)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,-D3.end+J.start,J.start,D3.end,(J.start>V.end),J.start+V.end))),]
    
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
    d.x1 = rep(-1,length(temp.agg3$D5.end))
    d.x2 = rep(-1,length(temp.agg3$D5.end))
    if(plot.D){
      d.x1 = temp.agg3$D5.end-1
      d.x2 = temp.agg3$D3.end}
    d.filter = d.x1 >0
    y1=c(0,cumsum(100*prop.table(temp.agg3$Umi.count))[-nrow(temp.agg3)])
    y2=cumsum(100*prop.table(temp.agg3$Umi.count))
    
    
    d1=data.frame(x1=c(v.x1,j.x1),x2=c(v.x2,j.x2),y1=c(y1,y1),y2=c(y2,y2), t=factor(c(rep(c('V','J'),each=length(v.x1))),levels = c("V","J")))
    d2=rbind(d1,data.frame(x1=nn.x1,x2=nn.x2,y1=y1[nn.filter],y2=y2[nn.filter], t=factor(rep("N",sum(nn.filter)),levels = c("V","J","N"))))
    
    d3 = rbind(d2,data.frame(x1 = d.x1[d.filter],x2 = d.x2[d.filter],y1=y1[d.filter],y2=y2[d.filter],t=factor(rep("D",sum(d.filter)),levels = c("V","J","N","D"))))
    colores = c(V.col,J.col,N.col,D.col)
    if(estetica.vj){d3 = rbind(d2,data.frame(x1=vj.x1,x2=vj.x2,y1=y1[vj.filter],y2=y2[vj.filter], t=factor(rep("VJ",sum(vj.filter)),levels = c("V","J","N","VJ"))))
    colores = c(V.col,J.col,N.col,VJ.col,D.col)}
    d4 = with(temp.agg3,data.frame(L = unique(L), Umi.count = rev(tapply(Umi.count,L,sum))))
    d4$y1 = c(0,cumsum(100*prop.table(d4$Umi.count))[-nrow(d4)])
    d4$y2 = c(cumsum(100*prop.table(d4$Umi.count)))
    
    
    labels.x = seq(0, lim, by = 3)
    if(every.other.x.label)labels.x[labels.x%%6 != 0 ] = ""
    
    
    p1 = ggplot() + 
      scale_x_continuous(name="CDR3 nucleotide length",breaks = seq(0, lim, by = 3),limits = c(0,lim +5),labels = labels.x) + 
      scale_y_continuous(name="Fraction of the repertoire",labels = function(x)round(x/100,2),breaks = seq(0, 100, by = 25)) +
      geom_rect(data=d3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=NA,size = 0.1, alpha=alfa) +
      geom_rect(data=d4, mapping=aes(xmin=0, xmax=L, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.3) +
      scale_fill_manual(values=colores) +
      theme_minimal() +
      ggtitle(label = titulo) +
      theme(panel.grid.major.x = element_line(color = "grey",size = 0.3,linetype = 2),
            panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
            text=element_text(family="Helvetica",size = 15),
            axis.ticks.y = element_line(),
            axis.text = element_text(size = 12),
      ) 
    
  }
  return(p1)
}


VDJ.plot.2 = function(df,estetica = 1,estetica.vj = F,Umi.lim =1,lim = 60,alfa = 0.55 ,titulo ='',V.col = 'royalblue',J.col ='red',N.col ='grey',D.col = "dark orchid",VJ.col = "purple",plot.D = T,every.other.x.label = T){
  #df = TCR.a.spp.VJ$ZF
  df$D5.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0] + 1
  df$D3.end[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0] - 1
  df$V.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0]
  df$J.start[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0]
  
  {
    {
      df.sel = df[,with(df,c("Umi.count","V.end","J.start","L","D5.end","D3.end"))]
      temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
      temp.agg2 = temp.agg[temp.agg$Umi.count>=Umi.lim,]
      if(estetica ==1) temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,J.start,V.end,J.start-V.end))),]
      if(estetica ==2)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,J.start+V.end))),]
      if(estetica ==3)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,(J.start>V.end),J.start+V.end))),]
      if(estetica ==4)temp.agg3 = temp.agg2[rev(with(temp.agg2,order(L,-D3.end+J.start,J.start,D3.end,(J.start>V.end),J.start+V.end))),]
      
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
      d.x1 = rep(-1,length(temp.agg3$D5.end))
      d.x2 = rep(-1,length(temp.agg3$D5.end))
      d.x1 = temp.agg3$D5.end-1
      d.x2 = temp.agg3$D3.end
      d.filter = d.x1 >0
      y1=c(0,cumsum(100*prop.table(temp.agg3$Umi.count))[-nrow(temp.agg3)])
      y2=cumsum(100*prop.table(temp.agg3$Umi.count))
      
      
      d1=data.frame(x1=c(v.x1,j.x1),x2=c(v.x2,j.x2),y1=c(y1,y1),y2=c(y2,y2), t=factor(c(rep(c('V','J'),each=length(v.x1))),levels = c("V","J")))
      d2=rbind(d1,data.frame(x1=nn.x1,x2=nn.x2,y1=y1[nn.filter],y2=y2[nn.filter], t=factor(rep("N",sum(nn.filter)),levels = c("V","J","N"))))
      
      d3 = rbind(d2,data.frame(x1 = d.x1[d.filter],x2 = d.x2[d.filter],y1=y1[d.filter],y2=y2[d.filter],t=factor(rep("D",sum(d.filter)),levels = c("V","J","N","D"))))
      colores = c(V.col,J.col,N.col,VJ.col)
      
      #d3 = rbind(d2,data.frame(x1=vj.x1,x2=vj.x2,y1=y1[vj.filter],y2=y2[vj.filter], t=factor(rep("VJ",sum(vj.filter)),levels = c("V","J","N","VJ"))))
      #colores = c(V.col,J.col,N.col,VJ.col,D.col)
      
      d4 = with(temp.agg3,data.frame(L = unique(L), Umi.count = rev(tapply(Umi.count,L,sum))))
      d4$y1 = c(0,cumsum(100*prop.table(d4$Umi.count))[-nrow(d4)])
      d4$y2 = c(cumsum(100*prop.table(d4$Umi.count)))
    }
    
    labels.x = seq(0, lim, by = 3)
    if(every.other.x.label)labels.x[labels.x%%6 != 0 ] = ""
    
    
    p1 = ggplot() + 
      scale_x_continuous(name="CDR3 nucleotide length",breaks = seq(0, lim, by = 3),limits = c(0,lim +5),labels = labels.x) + 
      scale_y_continuous(name="Fraction of the repertoire",labels = function(x)round(x/100,2),breaks = seq(0, 100, by = 25)) +
      geom_rect(data=d3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=NA,size = 0.1, alpha=alfa) +
      geom_rect(data=d4, mapping=aes(xmin=0, xmax=L, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.3) +
      scale_fill_manual(values=colores) +
      theme_minimal() +
      ggtitle(label = titulo) +
      theme(panel.grid.major.x = element_line(color = "grey",size = 0.3,linetype = 2),
            panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
            text=element_text(family="Helvetica",size = 15),
            axis.ticks.y = element_line(),
            axis.text = element_text(size = 12),
      ) 
    
  }
  return(p1)
}





