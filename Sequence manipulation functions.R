

tr.DNA = function(st)translate(DNAStringSet(subseq(as.character(st),1,width=floor(nchar(as.character(st))/3)*3-3)),if.fuzzy.codon = "solve", no.init.codon = TRUE)

tr.DNA3 = function(st)sapply(1:3,function(x)translate(DNAStringSet(subseq(st,x,width=floor(nchar(st)/3)*3-3)),if.fuzzy.codon = "solve", no.init.codon = TRUE))

tr.DNA6 = function(st)return(unlist(list(forward=tr.DNA3(st),reverse=tr.DNA3(reverseComplement(DNAStringSet(st))))))

exclude.alleles = function(sD,distancia = 1,remove.since.nth = 2){
  gb = graph.builder(which(matrix.builder(sD)<=distancia,arr.ind = T))
  return(unlist(sapply(gb,function(x)if(length(x)>=remove.since.nth)x[remove.since.nth:length(x)])))}

graph.split = function(sD,distancia = 1)lapply(graph.builder(which(matrix.builder(sD)<=distancia,arr.ind = T)),as.numeric)

graph.builder = function(nodes){
  conflicting.UMIs.list = list()
  if (length(nodes!=0)){
    c = as.vector(t(nodes))
    g <- graph(c)
    cl <- clusters(g)
    conflicting.UMIs.list <- lapply(seq_along(cl$csize)[cl$csize > 1], function(x) V(g)[cl$membership %in% x])
  } else{conflicting.UMIs.list=NULL}
  return (conflicting.UMIs.list)
}

matrix.builder = function(distancia){
  matriz = as.matrix(distancia)
  matriz[upper.tri(matriz)] <- NA
  diag(matriz) = NA
  return(matriz)  
}

char.cumsum = function(char.vector,p)as.numeric(names(which(cumsum(prop.table(rev(table(nchar(as.character(char.vector))))))>=p)[1]))

nodes.finder.UMI = function(secuencia){
  secuencia=as.character(secuencia)
  oR = lapply(1:14,function(i)paste0(subseq(secuencia,1,width = (i-1)),subseq(secuencia,i+1,width = (14-i))))
  nodos= rowSums(sapply(oR,function(x)duplicated(x)|duplicated(x,fromLast = TRUE)))>0
  return(nodos)
}

pc.generator.nt = function(df)paste(df$CDR3.nucleotide.sequence,df$V.gene,df$J.gene,sep = ";")

pc.chart = function(pc.list,i,j,k=0)sapply(pc.list[ sapply(pc.list,length)>=(i+k)],function(x)sapply(pc.list[ sapply(pc.list,length)>=(i+k)],function(y)length(intersect(x[sample((1+k):(i+k),j)],y[sample((1+k):(i+k),j)]))))

get.V.Seq = function(sp,gene,largo = 100,data.list = spp.genome.data.list,TCR.dict = spp.TCR.dict,xt = 0)getSeq(data.list[[sp]][[gene]],TCR.dict[[sp]][[gene]]$V.GR.dict)

get.V.end.Seq = function(sp,gene,largo = 100,data.list = spp.genome.data.list,TCR.dict = spp.TCR.dict,xt = 0)getSeq(data.list[[sp]][[gene]],flank(TCR.dict[[sp]][[gene]]$V.GR.dict-3-xt,start = F,width = largo))

get.V.RSS.Seq  = function(sp,gene,largo = 100)subseq(getSeq(spp.genome.data.list[[sp]][[gene]],flank(spp.TCR.dict[[sp]][[gene]]$V.GR.dict-3,start = F,width = largo)),spp.TCR.dict[[sp]][[gene]]$V.GR.dict$heptamer.start)

get.J.Seq = function(sp,gene,largo = 200)getSeq(spp.genome.data.list[[sp]][[gene]],resize(spp.TCR.dict[[sp]][[gene]]$J.GR.dict,width = largo,fix = "end"))

get.V.RSS.matrix = function(sp,gene,largo = 100)consensusMatrix(subseq(get.V.RSS.Seq(sp,gene,largo),1,38))[1:4,]

get.J.RSS.Seq = function(sp,gene,extend = 0)subseq(get.J.Seq(sp,gene),end = spp.TCR.dict[[sp]][[gene]]$J.GR.dict$heptamer.end+extend, width = 28+extend)

get.J.RSS.matrix = function(sp,gene,extend = 0)consensusMatrix(subseq(get.J.Seq(sp,gene),end = spp.TCR.dict[[sp]][[gene]]$J.GR.dict$heptamer.end+extend, width = 28+extend))[1:4,]

get.D.Seq = function(sp,gene,full = F)if(full)getSeq(spp.genome.data.list[[sp]][[gene]],full.D.converter( spp.TCR.dict[[sp]][[gene]]$D.GR.dict))else getSeq(spp.genome.data.list[[sp]][[gene]], spp.TCR.dict[[sp]][[gene]]$D.GR.dict)

full.D.converter = function(D.GR,left.RSS = 28, right.RSS = 39)resize(flank(D.GR,left.RSS),left.RSS+width(D.GR)+right.RSS)

kmer.generator = function(secuencias,largo)lapply(as.character(secuencias),function(x)sapply(1:(nchar(x)-(largo-1)),function(y)subseq(x,y,width = largo)))

dictionary.extractor = function(genome.data,V.GR){
  genomic.sequences = getSeq(genome.data,V.GR)
  genomic.dictionary = DNAStringSet(sapply(1:length(genomic.sequences),function(i)genomic.sequences[[i]][-start(V.GR[i]$introns):-end(V.GR[i]$introns)]))
  genomic.dictionary[width(genomic.dictionary)==0] = genomic.sequences[width(genomic.dictionary)==0]
  return(list(genomic.sequences=genomic.sequences,genomic.dictionary=genomic.dictionary))
}

nonstop = function(sec)sec[!grepl("\\*",sec)]

ns.filter = function(sec)!grepl("\\*",sec)
