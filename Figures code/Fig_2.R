# Fig 2 A 
{

  p.gg.ZF.V = ggseqlogo(consensusMatrix((get.V.end.Seq("ZF","a")[spp.TCR.dict[["ZF"]]$a$V.GR.dict$heptamer.start ==15]))[,1:14],col_scheme = cs.black) + 
    theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())
  
  ggsave(filename = paste0(spp.fig.folder,"/Fig2.A_V.pdf"),p.gg.ZF.V,width = 8,height = 2)  
  
  p.gg.ZF.J =   ggseqlogo(consensusMatrix(get.J.Seq("ZF","a")[spp.TCR.dict[["ZF"]]$a$J.GR.dict$heptamer.end == -66])[DNA_BASES,136:149],col_scheme = cs.black) +
    theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())
  
  ggsave(filename = paste0(spp.fig.folder,"/Fig2.A_J.pdf"),p.gg.ZF.J,width = 8,height = 2) 
  
}

# Fig 2 B
{
  {
    ex4 = data.frame(x1 = c(0,14,0,13,0,11,0,7,0,13,0,10,0,10,0,7),x2 = c(15,45,14,42,12,42,9,39,14,45,14,42,14,42,14,39) , y1 = c(c(0,0,7,7,14,14,21,21), c(0-4,0,7-4,7,14-4,14,21-4,21)+2), y2 = c(c(0,0,7,7,14,14,21,21)+1,c(0-4,0,7-4,7,14-4,14,21-4,21)+3), t = factor(c("V","J"),levels = c("V","J")))
    
    ex4.1 = data.frame(x1=c(14,15,15,14), y1 = c(0,0,1,1) - 2, t = factor(c("V"),levels = c("V","J")))
    ex4.2 = data.frame(x1=c(14,16,14), y1 = c(0,0.5,1) - 2 +7, t = factor(c("V"),levels = c("V","J")))
    ex4.3 = data.frame(x1=c(14,16,14), y1 = c(0,0.5,1) - 2 + 14, t = factor(c("V"),levels = c("V","J")))
    ex4.4 = data.frame(x1=c(14,16,14), y1 = c(0,0.5,1) - 2 + 21, t = factor(c("V"),levels = c("V","J")))
    ex4.5 = data.frame(x1=c(13,11,13), y1 = c(0,0.5,1) - 2 + 4, t = factor(c("J"),levels = c("V","J")))
    ex4.6 = data.frame(x1=c(13,11,13) - 3, y1 = c(0,0.5,1) - 2 +7 +4, t = factor(c("J"),levels = c("V","J")))
    ex4.7 = data.frame(x1=c(13,11,13) - 3, y1 = c(0,0.5,1) - 2 + 14 +4, t = factor(c("J"),levels = c("V","J")))
    ex4.8 = data.frame(x1=c(13,11,13) - 6, y1 = c(0,0.5,1) - 2 + 21 +4, t = factor(c("J"),levels = c("V","J")))
    
  }  
  
  p = ggplot() + 
    scale_x_continuous(name="x",breaks = seq(0, 50, by = 3)) + 
    scale_y_continuous(name="y") +
    scale_fill_manual(values=c(V.col,J.col)) +
    
    geom_rect(data=ex4, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=NA,size = 0.1, alpha=0.55) +
    geom_rect(data=ex4, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),fill = NA, color="black",size = 0.3, alpha=0.55) +
    
    geom_polygon(data=ex4.1, mapping=aes(x=x1,y=y1), fill=P.nt.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.2, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.3, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.4, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.5, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.6, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.7, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    geom_polygon(data=ex4.8, mapping=aes(x=x1,y=y1), fill=RSS.col, color='black',size = 0.1, alpha=0.55) +
    
    geom_segment(aes(x = 0, y = -2, xend = 0, yend = 1),linetype =3) +
    geom_segment(aes(x = 15, y = -2, xend = 15, yend = 1),linetype =3) +
    geom_segment(aes(x = 14, y = 0, xend = 14, yend = 3),linetype =3) +
    geom_segment(aes(x = 45, y = 0, xend = 45, yend = 3),linetype =3) +
    
    geom_segment(aes(x = 0, y = 5, xend = 0, yend = 5+3),linetype =3) +
    geom_segment(aes(x = 14, y = 5, xend = 14, yend = 5+3),linetype =3) +
    geom_segment(aes(x = 13, y = 10, xend = 13, yend = 10-3),linetype =3) +
    geom_segment(aes(x = 42, y = 10, xend = 42, yend = 10-3),linetype =3) +
    
    geom_segment(aes(x = 0, y = 12, xend = 0, yend = 12+3),linetype =3) +
    geom_segment(aes(x = 12, y = 12, xend = 12, yend = 12+3),linetype =3) +
    geom_segment(aes(x = 11, y = 17, xend = 11, yend = 17-3),linetype =3) +
    geom_segment(aes(x = 42, y = 17, xend = 42, yend = 17-3),linetype =3) +
    
    geom_segment(aes(x = 0, y = 19, xend = 0, yend = 12+7+3),linetype =3) +
    geom_segment(aes(x = 9, y = 19, xend = 9, yend = 12+7+3),linetype =3) +
    geom_segment(aes(x = 7, y = 24, xend = 7, yend = 17+7-3),linetype =3) +
    geom_segment(aes(x = 39, y = 24, xend = 39, yend = 17+7-3),linetype =3) +
    
    theme(panel.grid.major.x = element_line(color = "grey",size = 0.3,linetype = 2),
          panel.background = element_blank(), axis.line =  element_blank())
  
  
  p.Fig2b = p + theme(panel.grid.major.x  = element_blank(), 
                      panel.grid.minor.x  = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      axis.text.x = element_blank(),
                      legend.position = "none")
}
ggsave(paste0(spp.fig.folder,"Fig2.B.pdf"),p.Fig2b,width = 10,height = 10)
# Fig 2 C



{
  sp = "ZF"
  df = TCR.a.spp.VJ[[sp]]
  df$D5.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0] + 1
  df$D3.end[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0] - 1
  df$V.end[df$Total.insertions < 0] = df$V.end[df$Total.insertions < 0] + df$Total.insertions[df$Total.insertions < 0]
  df$J.start[df$Total.insertions < 0] = df$J.start[df$Total.insertions < 0] - df$Total.insertions[df$Total.insertions < 0]
  rules.col = c('white','sienna4','magenta','dark green','darkorange 1')
  J.RSS.col = c(J1.col,J2.col,J3.col,J4.col,J5.col)
  lim = 64
  Umi.lim = 1
  plot.D = F
  titulo = ''
  estetica = 5
  
  df = df[df$L %% 3 ==0,]
  df$rules[!df$rules %in% 1:4] = "r0"
  df$rules[df$rules %in% 1:4] = paste0("r",  df$rules[df$rules %in% 1:4])
  df$rules = factor(df$rules,levels = c("r1","r2","r3","r4","r0"))
  df.sel = df[,with(df,c("Umi.count","V.end","J.start","L","D5.end","D3.end",'rules','RSS.J'))]
  estetica.vj = T
  temp.agg = aggregate(cbind(Umi.count) ~ ., data = df.sel, FUN = sum, na.rm = TRUE)
  temp.agg2 = temp.agg[temp.agg$Umi.count>=Umi.lim,]
  temp.agg3 = temp.agg2[rev(with(temp.agg2,order(rules,L))),]
  
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
  
  rules.x1 = rep(60,nrow(temp.agg3))
  rules.x2 = rep(63,nrow(temp.agg3))
  rule = temp.agg3$rules
  d1=data.frame(x1=c(v.x1,j.x1),x2=c(v.x2,j.x2),y1=c(y1,y1),y2=c(y2,y2), t=factor(c(rep(c('V','J'),each=length(v.x1))),levels = c("V","J")))
  d2=rbind(d1,data.frame(x1=nn.x1,x2=nn.x2,y1=y1[nn.filter],y2=y2[nn.filter], t=factor(rep("N",sum(nn.filter)),levels = c("V","J","N"))))
  d3 = rbind(d2,data.frame(x1 = d.x1[d.filter],x2 = d.x2[d.filter],y1=y1[d.filter],y2=y2[d.filter],t=factor(rep("D",sum(d.filter)),levels = c("V","J","N","D"))))
  colores = c(V.col,J.col,N.col,VJ.col,rules.col,J.RSS.col)
  d4 = with(temp.agg3,data.frame(L = unique(rule), Umi.count = rev(tapply(Umi.count,rule,sum))))
  d4$y1 = c(0,cumsum(100*prop.table(d4$Umi.count))[-nrow(d4)])
  d4$y2 = c(cumsum(100*prop.table(d4$Umi.count)))
  d5 = rbind(d3,data.frame(x1 = rules.x1,x2 = rules.x2,y1 = y1,y2 = y2,t = factor(rule,levels = c("V","J","N","r0","r1","r2","r3","r4"))))
  
  labels.x = seq(0, lim, by = 3)
  labels.x[labels.x%%6 != 0 ] = ""
  {
    p.ZF.rules= ggplot() + 
      scale_x_continuous(name="CDR3 nucleotide length",breaks = seq(0, lim, by = 3), limits = c(0,lim + 10),labels= labels.x) + 
      scale_y_continuous(name="Fraction of the repertoire",labels = function(x)round(x/100,2)) +
      geom_rect(data=d5, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=NA,size = 0.1, alpha=0.66) +
      geom_rect(data=d4, mapping=aes(xmin=0, xmax=63, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.4) +
      geom_rect(data=d4, mapping=aes(xmin=60, xmax=63, ymin=y1, ymax=y2), fill=NA, color="black",size = 0.4, alpha = 1) +
      scale_fill_manual(values=colores) +
      theme_minimal() +
      theme(panel.grid.major.x = element_line(color = "grey",size = 0.3,linetype = 2),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks.y = element_line()) +
      ggtitle(label = titulo)
  }
  
}
ggsave(paste0(spp.fig.folder,"Fig2.C.pdf"),p.ZF.rules,width = 10,height = 10)



p.teleost.rules = ggplot(rules.df[rules.df$sp %in% c("MF","ZF","RT","MK"),]) + 
  geom_col(mapping = aes(x = factor(sp,levels =c("MF","ZF","RT","MK")), y = value , fill = factor(rules,levels = 1:4)), position = "stack", alpha = 0.7) +
  scale_color_manual(values = rules.col[2:5]) +
  scale_fill_manual(values = rules.col[2:5]) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(labels = c("P. progenetica","D. rerio","O. mykiss","O. latipes")) +
  
  ylab("Fraction of the repertoire explained") +
  xlab(element_blank()) +
  labs(fill = "pattern") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 15, face = 'italic'),
        axis.text.y = element_text(size = 12),
        axis.title =  element_text(size = 15)
  )               # Rotate axis labels
p.teleost.rules
ggsave(paste0(spp.fig.folder,"Fig2.D.pdf"),p.teleost.rules,width = 4,height = 6)
