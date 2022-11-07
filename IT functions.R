

TCR.entropy.VorJ = function(database,L.min = 18,L.max=72, threshold = 100){
  {
    #generate new variables for making it readable
    CDR3.nt = database$CDR3.nucleotide.sequence
    CDR3.aa = database$CDR3.amino.acid.sequence
    
    L=as.numeric(nchar(CDR3.nt))
    tL = table(L)
    Ls = sort(unique(L))[tL>=threshold & as.numeric(names(tL))>=L.min & as.numeric(names(tL)) <= L.max]
    Ls.aa = Ls[Ls%%3==0]
    V=database$V.gene
    Ve=database$V.end
    
    J=database$J.gene
    Js=database$J.start
  }
  if (length(Ls)==0) return(list())
  nt.H= list()
  nt.H.V= list()
  nt.H.J= list()
  nt.H.VJ = list()
  FXYs = list()
  consensus.nt = list()
  aa.H= list()
  aa.H.V= list()
  aa.H.J= list()
  aa.H.VJ = list()
  FXYs.aa = list()
  consensus.aa = list()
  
  #please make sure L is numeric!
  
  for (l in Ls){
    nt.H[[l]] = apply(consensusMatrix(CDR3.nt[L==l]),2,function(x)entropia(x[x!=0]))
    nt.H.V.table = tapply(CDR3.nt[L==l],V[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
    nt.H.V[[l]] = colSums(do.call(rbind,nt.H.V.table)/sum(L==l))
    nt.H.J.table = tapply(CDR3.nt[L==l],J[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
    nt.H.J[[l]] = colSums(do.call(rbind,nt.H.J.table)/sum(L==l))  
    FXYs[[l]] = t(apply(data.frame(nt.H[[l]]-nt.H.V[[l]],nt.H[[l]]-nt.H.J[[l]]),1,function(x)c(0,max(x))[c(which.min(x),which.max(x))])) # this line means: take 0 and the maximum mutual information, if the order in the dataframe is the same do nothing, if the maximum was on the left, flip
    consensus.nt[[l]] = strsplit(consensusString(CDR3.nt[L==l]),"")[[1]]
    
    
    if (l%%3 ==0){
      aa.H[[l]] = apply(consensusMatrix(CDR3.aa[L==l]),2,function(x)entropia(x[x!=0]))
      aa.H.V.table = tapply(CDR3.aa[L==l],V[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
      aa.H.V[[l]] = colSums(do.call(rbind,aa.H.V.table)/sum(L==l))
      aa.H.J.table = tapply(CDR3.aa[L==l],J[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
      aa.H.J[[l]] = colSums(do.call(rbind,aa.H.J.table)/sum(L==l))
      FXYs.aa[[l]] = t(apply(data.frame(aa.H[[l]]-aa.H.V[[l]],aa.H[[l]]-aa.H.J[[l]]),1,function(x)c(0,max(x))[c(which.min(x),which.max(x))])) # this line means: take 0 and the maximum mutual information, if the order in the dataframe is the same do nothing, if the maximum was on the left, flip
      consensus.aa[[l]] = strsplit(consensusString(CDR3.aa[L==l]),"")[[1]]
      
    }
  }
  entropy.nt.VorJ = sum(Hx(L[L%in%Ls]),sum(sapply(Ls,function(l)(sum(nt.H[[l]])-sum( FXYs[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls]==l))))
  #entropy.nt.VorJ.c = c(Hx(L[L%in%Ls]),c(sapply(Ls,function(l)(sum(nt.H[[l]])-sum( FXYs[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls]==l))))
  entropy.nt.VorJ.c = c(Hx(L[L%in%Ls]),c(sapply(Ls,function(l)(sum(nt.H[[l]])-sum( FXYs[[l]])+Hxy(V[L==l],J[L==l])))))
  names(entropy.nt.VorJ.c) = c("L",Ls)
  if (length (Ls.aa)>0){
    entropy.aa.VorJ = sum(Hx(L[L%in%Ls.aa]),sum(sapply(Ls[Ls%%3 ==0],function(l)(sum(aa.H[[l]])-sum( FXYs.aa[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls[Ls%%3 ==0]]==l))))
    #entropy.aa.VorJ.c = c(Hx(L[L%in%Ls & L%%3 ==0]),c(sapply(Ls[Ls%%3 ==0],function(l)(sum(aa.H[[l]])-sum( FXYs.aa[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls[Ls%%3 ==0]]==l))))
    entropy.aa.VorJ.c = c(Hx(L[L%in%Ls.aa]),c(sapply(Ls[Ls%%3 ==0],function(l)(sum(aa.H[[l]])-sum( FXYs.aa[[l]])+Hxy(V[L==l],J[L==l])))))
    names(entropy.aa.VorJ.c) = c("L",Ls.aa)
    
  } else {entropy.aa.VorJ = c()
  entropy.aa.VorJ.c = c()}
  
  res.nt = list(H.CDR3 = nt.H, H.V = nt.H.V, H.J = nt.H.J,FXYs = FXYs, consensus = consensus.nt)
  res.aa = list(H.CDR3 = aa.H, H.V = aa.H.V, H.J = aa.H.J, FXYs.aa = FXYs.aa, consensus = consensus.aa)
  
  return(list(nt=apply(do.call(cbind,res.nt),1,as.list),aa=apply(do.call(cbind,res.aa),1,as.list),statistics = list(entropy.nt.VorJ = entropy.nt.VorJ, entropy.aa.VorJ = entropy.aa.VorJ, entropy.nt.VorJ.c = entropy.nt.VorJ.c, entropy.aa.VorJ.c = entropy.aa.VorJ.c)))
}

TCR.entropy.VorJ.umi = function(database){
  {
    #generate new variables for making it readable
    clone.nr = database$Umi.count
    CDR3.nt = rep(database$CDR3.nucleotide.sequence,clone.nr)
    CDR3.aa = rep(database$CDR3.amino.acid.sequence,clone.nr)
    
    
    L=nchar(CDR3.nt)
    Ls = sort(unique(L))[table(L)>100]
    Ls.aa = Ls[Ls%%3==0]
    V=rep(database$V.gene,clone.nr)
    Ve=rep(database$V.end,clone.nr)
    
    J=rep(database$J.gene,clone.nr)
    Js=rep(database$J.start,clone.nr)
  }
  nt.H= list()
  nt.H.V= list()
  nt.H.J= list()
  nt.H.VJ = list()
  FXYs = list()
  consensus.nt = list()
  aa.H= list()
  aa.H.V= list()
  aa.H.J= list()
  aa.H.VJ = list()
  FXYs.aa = list()
  consensus.aa = list()
  
  #please make sure L is numeric!
  
  for (l in Ls){
    nt.H[[l]] = apply(consensusMatrix(CDR3.nt[L==l]),2,function(x)entropia(x[x!=0]))
    nt.H.V.table = tapply(CDR3.nt[L==l],V[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
    nt.H.V[[l]] = colSums(do.call(rbind,nt.H.V.table)/sum(L==l))
    nt.H.J.table = tapply(CDR3.nt[L==l],J[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
    nt.H.J[[l]] = colSums(do.call(rbind,nt.H.J.table)/sum(L==l))  
    FXYs[[l]] = t(apply(data.frame(nt.H[[l]]-nt.H.V[[l]],nt.H[[l]]-nt.H.J[[l]]),1,function(x)c(0,max(x))[c(which.min(x),which.max(x))])) # this line means: take 0 and the maximum mutual information, if the order in the dataframe is the same do nothing, if the maximum was on the left, flip
    consensus.nt[[l]] = strsplit(consensusString(CDR3.nt[L==l]),"")[[1]]
    
    
    if (l%%3 ==0){
      aa.H[[l]] = apply(consensusMatrix(CDR3.aa[L==l]),2,function(x)entropia(x[x!=0]))
      aa.H.V.table = tapply(CDR3.aa[L==l],V[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
      aa.H.V[[l]] = colSums(do.call(rbind,aa.H.V.table)/sum(L==l))
      aa.H.J.table = tapply(CDR3.aa[L==l],J[L==l],function(x)as.numeric(length(x)*apply(consensusMatrix(x),2,function(y)entropia(y[y!=0]))))
      aa.H.J[[l]] = colSums(do.call(rbind,aa.H.J.table)/sum(L==l))
      FXYs.aa[[l]] = t(apply(data.frame(aa.H[[l]]-aa.H.V[[l]],aa.H[[l]]-aa.H.J[[l]]),1,function(x)c(0,max(x))[c(which.min(x),which.max(x))])) # this line means: take 0 and the maximum mutual information, if the order in the dataframe is the same do nothing, if the maximum was on the left, flip
      consensus.aa[[l]] = strsplit(consensusString(CDR3.aa[L==l]),"")[[1]]
      
    }
  }
  entropy.nt.VorJ = sum(Hx(L[L%in%Ls]),sum(sapply(Ls,function(l)(sum(nt.H[[l]])-sum( FXYs[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls]==l))))
  entropy.nt.VorJ.c = c(Hx(L[L%in%Ls]),c(sapply(Ls,function(l)(sum(nt.H[[l]])-sum( FXYs[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls]==l))))
  if (length (Ls.aa)>0){
    entropy.aa.VorJ = sum(Hx(L[L%in%Ls & L%%3 ==0]),sum(sapply(Ls[Ls%%3 ==0],function(l)(sum(aa.H[[l]])-sum( FXYs.aa[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls[Ls%%3 ==0]]==l))))
    entropy.aa.VorJ.c = c(Hx(L[L%in%Ls & L%%3 ==0]),c(sapply(Ls[Ls%%3 ==0],function(l)(sum(aa.H[[l]])-sum( FXYs.aa[[l]])+Hxy(V[L==l],J[L==l]))*mean(L[L%in%Ls[Ls%%3 ==0]]==l))))
  } else {entropy.aa.VorJ = c()
  entropy.aa.VorJ.c = c()}
  
  res.nt = list(H.CDR3 = nt.H, H.V = nt.H.V, H.J = nt.H.J,FXYs = FXYs, consensus = consensus.nt)
  res.aa = list(H.CDR3 = aa.H, H.V = aa.H.V, H.J = aa.H.J, FXYs.aa = FXYs.aa, consensus = consensus.aa)
  
  return(list(nt=apply(do.call(cbind,res.nt),1,as.list),aa=apply(do.call(cbind,res.aa),1,as.list),statistics = list(entropy.nt.VorJ = entropy.nt.VorJ, entropy.aa.VorJ = entropy.aa.VorJ, entropy.nt.VorJ.c = entropy.nt.VorJ.c, entropy.aa.VorJ.c = entropy.aa.VorJ.c)))
}

entropia = function(tabla1)sum(-log2(tabla1/sum(tabla1))*tabla1/sum(tabla1))

Hx = function(x)entropia(table(x))

Hxy = function(x,y)entropia(table(paste(x,y)))

general.entropy.report = function(db){
  el = entropia(table(nchar(as.character(db$CDR3.nucleotide.sequence))))
  ev = entropia(table(db$V.gene))
  ej = entropia(table(db$J.gene))
  evj = entropia(table(paste(db$V.gene,db$J.gene)))
  elv = entropia(table(paste(nchar(as.character(db$CDR3.nucleotide.sequence)),db$V.gene)))
  elj = entropia(table(paste(nchar(as.character(db$CDR3.nucleotide.sequence)),db$J.gene)))
  return (round(data.frame (H.length = el, H.V = ev, H.J = ej, MI.V_J = ev+ej-evj, MI.l_V = el+ev-elv,MI.l_J = el+ej-elj),2))
}

entropy.output.VorJ = function(TCR.entropy.VorJ.x,plot = 1,entropy.limit = 2,titulo = NULL, leyenda = F,letras = NULL,cex.letras = 1){#letras is a consensus string from sp
  H.CDR3 = TCR.entropy.VorJ.x[[1]]
  H.V = TCR.entropy.VorJ.x[[2]]
  H.J = TCR.entropy.VorJ.x[[3]]
  FXY = TCR.entropy.VorJ.x[[4]]
  if (leyenda)letras = TCR.entropy.VorJ.x[[5]]
  X = H.CDR3 - H.V #(I.xv)
  Y = H.CDR3 - H.J #(I.xj)
  XdY = H.CDR3
  
  
  if (plot==1){barplot(t(as.matrix(data.frame(XdY))),col=c("grey"),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  if (plot==2){barplot(t(as.matrix(data.frame(X,XdY-X))),col=c("darkblue","grey"),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  if (plot==3){barplot(t(as.matrix(data.frame(Y,XdY-Y))),col=c("red","grey"),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  #      if (plot==4){barplot(t(as.matrix(data.frame(FXY[,1],FXY[,2],XdY-rowSums(FXY)))),col=(c("darkblue","red","gray")),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  if (plot==4){barplot(t(as.matrix(data.frame(FXY[,1],FXY[,2],XdY-rowSums(FXY)))),col=(alpha(c("darkblue","red","gray"),alpha = 1)),names.arg = letras,ylim = c(0,entropy.limit),main = titulo,las = 1,cex.names = cex.letras, cex.axis = cex.letras)}
  
  #if (plot==5){barplot(t(as.matrix(data.frame(XY,(XdY-XY)))),col=c("purple","grey"),names.arg = letras,ylim = c(0,entropy.limit))}
  if (plot==7){barplot(t(as.matrix(data.frame(FXY[,1]*(entropy.limit/XdY),FXY[,2]*(entropy.limit/XdY),(XdY-rowSums(FXY))*(entropy.limit/XdY)))),col=(c("darkblue","red","gray")),names.arg = letras,ylim = c(0,entropy.limit))}
  if (plot==5){barplot(t(as.matrix(data.frame(XdY-rowSums(FXY)))),col=c("grey"),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  if (plot==8){barplot(t(as.matrix(data.frame(XdY-XY))),col=c("grey"),names.arg = letras,ylim = c(0,entropy.limit))}
  if (plot==6){barplot(t(as.matrix(data.frame(XdY-rowSums(FXY),FXY[,1],FXY[,2]))),col=(c("gray","darkblue","red")),names.arg = letras,ylim = c(0,entropy.limit),main = titulo)}
  
  
}

entropy.content = function(csmat){
  proportions = sapply(1:ncol(csmat),function(col)csmat[,col]/colSums(csmat)[col])
  entropy = (-proportions*log2(proportions))
  entropy[is.na(entropy)]=0
  return (colSums(entropy))
}

for(x in list.files(paste0(spp.data.folder,"IT"))) load(file = paste0(paste0(spp.data.folder,"IT/"),x))

