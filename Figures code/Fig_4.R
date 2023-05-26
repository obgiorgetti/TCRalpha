sel.spp =c("Danio rerio","Mus musculus", "Protopterus annectens","Acipenser ruthenus","Chiloscyllium punctatum","Gallus gallus",'Oncorhynchus mykiss',"Polypterus senegalus","Oryzias latipes")
spp.df.2 = spp.df[!is.na(spp.df$J_frame),]
spp.df.2 = spp.df.2[!duplicated(spp.df.2$spp),]

spp.df.2$orden =  order(order(as.numeric(spp.df.2$J_frame),decreasing = T))


Fig4.a = ggplot(spp.df.2) + 
  scale_y_continuous(name="J RSS frame conservation",limits = c(0,1.07),breaks = seq(0,1,by = 0.25)) + 
  scale_x_continuous(name="",limits = c(0,312)) + 
  geom_point(mapping = aes(x = orden,y=J_frame,color = group), size = 1) +
  
  geom_segment(mapping = aes(x = orden, xend = orden, y=J_frame,yend=J_frame+0.05-0.16*grepl("O",spp)-0.11*grepl("Or",spp)+0.05*grepl("Pr|Da|Ch",spp)), data = ~ subset(.,spp %in%sel.spp)) +
  geom_label(mapping = aes(x = orden,y=J_frame+0.05-0.16*grepl("O",spp)-0.11*grepl("Or",spp)+0.05*grepl("Pr|Da|Ch",spp),label = spp),hjust="left", size = 4.5, data = ~ subset(.,spp %in%sel.spp)) +
  
  scale_color_manual(values=spp.col.xt) +
  scale_fill_manual(values=spp.col.xt) +
  
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size= 15))
ggsave(paste0(spp.fig.folder,"Fig4.a.pdf"),Fig4.a,width = 16,height = 8,bg = "white")

spp.df.fig4 = spp.df

levels(spp.df$group)[levels(spp.df$group) == "Transition"] = "Lobe-finned fishes"
{
  J.n.label = with(spp.df.fig4[!is.na(spp.df.fig4$J_frame),],tapply(J_nr,group,sum))
  J.n.label = paste(J.n.label,"Js")
  J.n.labels = J.n.label[spp.df.fig4$group[!is.na(spp.df$J_frame)]]
  J.n.labels[duplicated(J.n.labels)] = ""
  
  sp.n.label = with(spp.df.fig4[!is.na(spp.df.fig4$J_frame),],tapply(J_nr,group,length))
  sp.n.label = paste("n=",sp.n.label,"species")
  sp.n.labels = sp.n.label[spp.df.fig4$group[!is.na(spp.df$J_frame)]]
  sp.n.labels[duplicated(J.n.labels)] = ""
  
  
  p.box.f4b = ggplot(spp.df.fig4[!is.na(spp.df.fig4$J_frame),]) + geom_boxplot(mapping = aes(x = group, y = J_frame,fill = group)) +
    geom_rect(mapping = aes(xmin = 0, xmax=9, ymin= 1,ymax = 1.25),fill = "white") +
    geom_text(mapping = aes(x = group, y  = 1.20,label = J.n.labels),size = 5) +
    geom_text(mapping = aes(x = group, y  = 1.05,label = sp.n.labels)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,1.25)) +
    scale_fill_manual(values=c('cadetblue2',"red",'darkseagreen 2',"royal blue",'dark orchid','chocolate1','coral1','brown3')[c(3,1,4,5,6,7,8,2)]) +
    theme_minimal() +
    labs(y = "RSS J frame conservation") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_line(),
          #plot.title = element_text(size = 12, face = "bold"),
          panel.border = element_rect(color = "black",fill = NA,size = 1),
          axis.text.y = element_text(size = 16),
          axis.title=element_text(size=18),
          legend.text = element_text(size = 12),
          legend.title = element_blank())
  

  ggsave(paste0(spp.fig.folder,"Fig4.b.pdf"),p.box.f4b,width = 16,height = 8,bg = "white")
}


{
  selected.Js.group.sp = spp.df$group[match( gsub("_"," ",gsub(".*/|_alpha.*|_[0-9]|_HDR|_HNI|_HSOK","",names(selected.Js))),spp.df$common.name)]
  selected.Js.by.group = lapply(levels(spp.df$group),function(x)spp.RSS.xt[selected.Js.group.sp %in% x])
  names(selected.Js.by.group) = levels(spp.df$group)
  selected.Js.by.group.cM = lapply(selected.Js.by.group,function(x)consensusMatrix(unlist(lapply(x,function(y)y$ss.RSS)))[c("A","C","G","T"),])
  
  
  pdf(paste0(spp.fig.folder,"Fig4_B_consensus.pdf"),width = 20,height = 2)  
  
  grid.arrange(grobs = lapply(1:8,function(i)ggseqlogo(selected.Js.by.group.cM[[i]][,29:33],col_scheme = cs.spp[[i]]) +
                                scale_y_continuous(limits = c(0,2)) +
                                #ggtitle(paste(colSums(selected.Js.by.group.cM[[i]])[[1]],"Js from",length(selected.Js.by.group[[i]]),"species") ) +
                                theme(axis.ticks.x = element_blank(),
                                      axis.text.x=element_blank(), 
                                      axis.title.y = element_blank(),
                                      axis.text.y=element_blank(),
                                      #plot.title = element_text(size = 12, face = "bold"),
                                      panel.border = element_rect(color = "black",fill = NA,size = 0.3))),nrow = 1)
  
  
  dev.off()
}

