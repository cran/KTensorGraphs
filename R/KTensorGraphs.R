# funciones auxiliares:

read<-function(X,Y=NULL)
{
  filas<-dim(X)[1]
  columnas<-dim(X)[2]
  if((length(dim(X)))==3){repeticiones<-dim(X)[3]}else{repeticiones<-NULL}
  if(!is.null(dimnames(X)[[1]])){namesf<-dimnames(X)[[1]]}else{namesf<-paste(1:filas,sep="")}
  if(!is.null(dimnames(X)[[2]])){namesc<-dimnames(X)[[2]]}else{namesc<-paste("V",1:columnas,sep="")}
  if((length(dim(X)))==3)
  {
    if(!is.null(dimnames(X)[[3]])){namesr<-dimnames(X)[[3]]}else{namesr<-paste("R",1:repeticiones,sep="")}
  }else{namesr<-NULL}
  nas<-apply(is.na(X),1,function(v){return(all(!v))})
  columnas2<-NULL
  namesc2<-NULL
  conf<-FALSE
  if(!is.null(Y))
  {
    columnas2<-dim(Y)[2]
    if(!is.null(dimnames(Y)[[2]])){namesc2<-dimnames(Y)[[2]]}else{namesc2<-paste("V",1:columnas2,sep="")}
    if(filas!=dim(Y)[1])
    {
      conf<-TRUE
      print("te has confundido, el numero de filas de las dos matrices es distinto")
    }
    if((length(dim(X)))==3)
    {
      if(repeticiones!=dim(Y)[3])
      {
        conf<-TRUE
        print("te has confundido, el numero de repeticiones de los dos cubos es distinto")
      }
    }
    nasY<-apply(is.na(Y),1,function(v){return(all(!v))})
    nas<-nas&nasY
  }
  if(!all(nas))
  {
    if((length(dim(X)))==2){X<-X[nas,]}else{X<-X[nas,,]}
    if((length(dim(Y)))==2){Y<-Y[nas,]}else{Y<-Y[nas,,]}
    namesf<-namesf[nas]
    print(paste("se han suprimido algunas filas (",filas-sum(nas),") con datos faltantes",sep=""))
    filas<-sum(nas)
  }
  return(list(X,filas,columnas,repeticiones,namesf,namesc,namesr,Y,columnas2,namesc2,conf))
}

colores<-function(filas,columnas,coloresf,coloresc,conf,columnas2=NULL,coloresc2=NULL)
{
  if(is.null(coloresf)){coloresf<-rep("black",filas)}
  if(is.null(coloresc)){coloresc<-rep("black",columnas)}
  if((length(coloresf))!=filas)
  {
    conf<-TRUE
    print("te has confundido, el numero de etiquetas de colores de las filas es distinto del numero de filas")
  }
  if((length(coloresc))!=columnas)
  {
    conf<-TRUE
    print("te has confundido, el numero de etiquetas de colores de las columnas es distinto del numero de columnas")
  }
  if(is.null(columnas2))
  {
    return(list(coloresf,coloresc,conf))
  } else {
    if(is.null(coloresc2)){coloresc2<-rep("black",columnas2)}
    if((length(coloresc2))!=columnas2)
    {
      conf<-TRUE
      print("te has confundido, el numero de etiquetas de colores de las columnas de la segunda matriz
            es distinto del numero de columnas de la segunda matriz")
    }
    return(list(coloresf,coloresc,coloresc2,conf))
    }
}

preproc<-function(X,norm,cubo=FALSE,capas=FALSE)
{
  if(!cubo)
  {
    X<-apply(X,2,function(v){return(v-mean(v))})
    if(norm){X<-apply(X,2,function(v){return(v/sqrt(mean(v^2)))})}
  } else {
    filas<-dim(X)[1]
    columnas<-dim(X)[2]
    repeticiones<-dim(X)[3]
    X<-apply(X,2:3,function(v){return(v-mean(v))})
    if(norm)
    {
      if(!capas){X<-apply(X,2:3,function(v){return(v/sqrt(mean(v^2)))})
      } else {
        X<-aperm(array(apply(X,2,function(m){return(m/sqrt(mean(m^2)))}),
                       dim=c(filas,repeticiones,columnas)),c(1,3,2))
      }
    }
  }
  X[is.nan(X)]<-0
  return(X)
}

pesos<-function(t)
{
  Dn<-diag((1/t)/(sum(1/t)))
  Dn2<-diag(sqrt(diag(Dn)))
  return(list(Dn,Dn2))
}

contributions<-function(M,names,title,first=FALSE)
{
  A<-t(round(1000*apply(M,1,function(v){return((v^2)/(sum(v^2)))})))
  if(dim(M)[2]==1){A<-t(A)}
  rownames(A)<-names
  if(!first){write(paste("\n",paste(rep("-",100),collapse=""),"\n",sep=""),file="results.txt",append=TRUE)}
  write(paste("\n",paste(rep(" ",20),collapse=""),"Contribuciones de ",title,"\n",sep=""),file="results.txt",
        append=!first)
  write(paste("Axis",paste(1:(dim(M)[2]),collapse="\t"),sep="\t"),file="results.txt",append=TRUE)
  write.table(A,file="results.txt",append=TRUE,sep="\t",col.names=FALSE)
}

contributions2<-function(M,dimension)
{
  A<-t(round(1000*apply(M,1,function(v){return((v^2)/(sum(v^2)))})))
  if(dim(M)[2]==1){A<-t(A)}
  A<-matrix(A[,dimension],ncol=length(dimension))
  return(apply(A,1,sum))
}

screeplot<-function(d)
{
  e<-100*d^2/sum(d^2)
  barplot(e,space=0.2,names.arg=1:length(d),cex.names=0.7)
  text((1:length(d))+0.2*(1:length(d))-0.5,0.5,paste(round(e,3),"%",sep=""),adj=0,srt=90)
}

compr<-function(dimX,dimY,d,c=NULL,tucker=FALSE)
{
  if(!tucker)
  {
    if((dimX!=(floor(dimX)))||(dimY!=(floor(dimY)))||(dimX>(length(d)))||(dimY>(length(d)))||(dimX>=dimY))
    {
      print("te has confundido, ejes incorrectos")
      return(FALSE)
    } else {return(TRUE)}
  } else {
    if((dimX!=(floor(dimX)))||(dimY!=(floor(dimY)))||(dimX>d)||(dimY>d)||((c!=1)&&(dimX>=dimY))||
       ((c==1)&&((dimX!=1)||(dimY!=1))))
    {
      print("te has confundido, ejes incorrectos")
      return(FALSE)
    } else {return(TRUE)}
  }
}

# funcion general para representar todos los graficos de los distintos analisis, excepto para el Tucker3 y el CoTucker3

plotm<-function(dim,d,M1,M2,M3=NULL,M4=NULL,M5=NULL,M6=NULL,lim1,lim2,names1,names2,colores1,colores2,contf,contc,
                cotucker=FALSE)
{
  contf[is.nan(contf)]<-0
  contc[is.nan(contc)]<-0
  mcontf<-max(contf)/2
  mcontc<-max(contc)/2
  s<-function(i,M,colores){segments(M[,1,i],M[,2,i],M[,1,i+1],M[,2,i+1],col=colores)}
  dev.new(width=14,height=7,noRStudioGD=TRUE)
  layout(matrix(1:2,nrow=1))
  layout.show(2)
  e<-round((100*d^2/sum(d^2))[dim],3)
  plot(M1[,1],M1[,2],type="n",xlim=c(min(0,lim1[,1]),max(0,lim1[,1])),ylim=c(min(0,lim1[,2]),max(0,lim1[,2])),
       xlab=paste("Axis ",dim[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dim[2]," (",e[2],"%)",sep=""))
  text(M1[,1],M1[,2],names1,cex=0.8,col=colores1)
  box(lwd=2)
  abline(h=0,lwd=2)
  abline(v=0,lwd=2)
  if(!cotucker)
  {
    if(!is.null(M3)){arrows(M1[,1],M1[,2],M3[,1],M3[,2],col=colores1,angle=10,length=0.1)}
  } else {
    if(!is.null(M3))
    {
      cual1<-cbind(M1[,1],M3[,1])
      cual2<-cbind(M1[,2],M3[,2])
      cual<-(apply(cual1,1,function(v){return(abs(v[1]-v[2])>1e-15)}))|
        (apply(cual2,1,function(v){return(abs(v[1]-v[2])>1e-15)}))
      if(sum(cual)!=0){arrows(M1[cual,1],M1[cual,2],M3[cual,1],M3[cual,2],col=colores1[cual],angle=10,length=0.1)}
    }
  }
  if(!is.null(M4)){lapply(1:(dim(M4)[3]-1),s,M=M4,colores=colores1)}
  if(!is.null(M6))
  {
    lcolores1<-colores1[!apply(M6==0,1,all)]
    lM6<-M6[!apply(M6==0,1,all),]
    arrows(0,0,lM6[,1],lM6[,2],col=lcolores1,angle=10,length=0.1)
  }
  plot(M2[,1],M2[,2],type="n",xlim=c(min(0,lim2[,1]),max(0,lim2[,1])),ylim=c(min(0,lim2[,2]),max(0,lim2[,2])),
       xlab=paste("Axis ",dim[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dim[2]," (",e[2],"%)",sep=""))
  text(M2[,1],M2[,2],names2,cex=0.8,col=colores2)
  box(lwd=2)
  abline(h=0,lwd=2)
  abline(v=0,lwd=2)
  lcolores2<-colores2[!apply(M2==0,1,all)]
  if(!cotucker){lM2<-M2[!apply(M2==0,1,all),]} else {lM2<-M2[!apply(abs(M2)<1e-15,1,all),]}
  arrows(0,0,lM2[,1],lM2[,2],col=lcolores2,angle=10,length=0.1)
  if(!is.null(M5)){lapply(1:(dim(M5)[3]-1),s,M=M5,colores=colores2)}
  h<-function(m){return(abs(c(min(0,m[,1]),max(0,m[,1]),min(0,m[,2]),max(0,m[,2]))))}
  k<-h(lim1)/h(lim2)
  k<-min(k[k!=0])
  dev.new()
  plot(M1[,1],M1[,2],type="n",xlim=c(min(0,lim1[,1],k*lim2[,1]),max(0,lim1[,1],k*lim2[,1])),
       ylim=c(min(0,lim1[,2],k*lim2[,2]),max(0,lim1[,2],k*lim2[,2])),
       xlab=paste("Axis ",dim[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dim[2]," (",e[2],"%)",sep=""))
  lcolores1<-colores1[contf>=mcontf]
  lnames1<-names1[contf>=mcontf]
  lM1<-matrix(M1[contf>=mcontf,],ncol=2)
  box(lwd=2)
  abline(h=0,lwd=2)
  abline(v=0,lwd=2)
  title(main=paste("Elementos con mayor contribucion.\nVariables multiplicadas por una constante (",round(k,3),
                   ")\npara que sean representables conjuntamente.",sep=""),cex.main=0.75,font.main=1)
  if(!cotucker)
  {
    if(!is.null(M3))
    {
      lM3<-matrix(M3[contf>=mcontf,],ncol=2)
      arrows(lM1[,1],lM1[,2],lM3[,1],lM3[,2],col=lcolores1,angle=10,length=0.1)
    }
  } else {
    if(!is.null(M3))
    {
      lM3<-matrix(M3[contf>=mcontf,],ncol=2)
      cual1<-cbind(lM1[,1],lM3[,1])
      cual2<-cbind(lM1[,2],lM3[,2])
      cual<-(apply(cual1,1,function(v){return(abs(v[1]-v[2])>1e-15)}))|
        (apply(cual2,1,function(v){return(abs(v[1]-v[2])>1e-15)}))
      if(sum(cual)!=0)
      {
        arrows(lM1[cual,1],lM1[cual,2],lM3[cual,1],lM3[cual,2],col=lcolores1[cual],angle=10,length=0.1)
      }
    }
  }
  if(!is.null(M4))
  {
    lM4<-array(M4[contf>=mcontf,,],dim=c(length(lcolores1),2,dim(M4)[3]))
    lapply(1:(dim(lM4)[3]-1),s,M=lM4,colores=lcolores1)
  }
  if(!is.null(M6))
  {
    lM6<-matrix(lM6[contf>=mcontf,],ncol=2)
    lcolores1<-lcolores1[!apply(lM6==0,1,all)]
    lM6<-matrix(lM6[!apply(lM6==0,1,all),],ncol=2)
    arrows(0,0,lM6[,1],lM6[,2],col=lcolores1,angle=10,length=0.1)
  }
  lcolores2<-colores2[contc>=mcontc]
  lnames2<-names2[contc>=mcontc]
  lM2<-matrix(M2[contc>=mcontc,],ncol=2)
  lcolores2<-lcolores2[!apply(lM2==0,1,all)]
  lnames2<-lnames2[!apply(lM2==0,1,all)]
  lM2<-matrix(lM2[!apply(lM2==0,1,all),],ncol=2)
  text(k*lM2[,1],k*lM2[,2],lnames2,cex=0.8,col=lcolores2)
  arrows(0,0,k*lM2[,1],k*lM2[,2],col=lcolores2,angle=10,length=0.1)
  if(!is.null(M5))
  {
    lM5<-array(M5[contc>=mcontc,,],dim=c(length(lcolores2),2,dim(M5)[3]))
    lapply(1:(dim(lM5)[3]-1),s,M=k*lM5,colores=lcolores2)
  }
  if(is.null(M6))
  {
    points(lM1[,1],lM1[,2],pch=19,col=lcolores1,cex=0.75)
    identify(lM1[,1],lM1[,2],lnames1,cex=0.8,tolerance=3)
  }else{text(lM1[,1],lM1[,2],lnames1,cex=0.8,col=lcolores1)}
}

tensorial<-function(X,dimension,U)
{
  m<-matrix(c(1,2,3,2,1,3,3,2,1),nrow=3,byrow=TRUE)
  X<-aperm(X,m[dimension,])
  return(aperm(array(apply(X,2:3,function(v,b){return(t(b)%*%v)},b=U),dim=c(dim(U)[2],dim(X)[2:3])),
               m[dimension,]))
}

desplegar<-function(X,dimension)
{
  m<-matrix(c(1,2,3,2,1,3,3,1,2),nrow=3)
  return(matrix(aperm(X,m[,dimension]),nrow=dim(X)[dimension]))
}

results<-function(X,title,names=NULL,axis=FALSE,first=FALSE)
{
  if(!is.null(names))
  {
    X<-t(X)
    colnames(X)<-names
    X<-t(X)
  }
  if(!first){write(paste("\n",paste(rep("-",100),collapse=""),"\n",sep=""),file="results.txt",append=TRUE)}
  write(paste("\n",paste(rep(" ",20),collapse=""),title,"\n",sep=""),file="results.txt",append=!first)
  if(axis){write(paste("Axis",paste(1:(dim(X)[2]),collapse="\t"),sep="\t"),file="results.txt",append=TRUE)}
  write.table(X,file="results.txt",append=TRUE,sep="\t",col.names=FALSE)
}

# funcion auxiliar para realizar el TUCKER3 fijado el numero de componentes para cada dimension

every<-function(p,q,r,T,iter,tol,mas=FALSE)
{

  # calcula los vectores singulares por la izquierda de los desplegamientos del tensor a lo largo de cada una
  # de las dimensiones

  A1<-(svd(desplegar(T,1),nu=p,nv=0))$u
  A2<-A1
  B1<-(svd(desplegar(T,2),nu=q,nv=0))$u
  B2<-B1
  C1<-(svd(desplegar(T,3),nu=r,nv=0))$u
  C2<-C1

  # calcula el tensor core

  G1<-tensorial(tensorial(tensorial(T,1,A1),2,B1),3,C1)
  G2<-G1

  # paso de iteracion:

  i<-function(l)
  {
    seguir<-l[[14]]
    rep<-l[[15]]
    iter<-l[[16]]
    if((seguir)&&(rep<iter))
    {
      T<-l[[1]]
      p<-l[[2]]
      q<-l[[3]]
      r<-l[[4]]
      A1<-l[[5]]
      B1<-l[[6]]
      C1<-l[[7]]
      G1<-l[[8]]
      A2<-l[[9]]
      B2<-l[[10]]
      C2<-l[[11]]
      G2<-l[[12]]
      tol<-l[[13]]

      # calcula la matriz optima fijando los vectores para las dimensiones dos y tres

      A<-(svd(desplegar(tensorial(tensorial(T,2,B1),3,C1),1),nu=p,nv=0))$u

      # calcula la matriz optima fijando los vectores para las dimensiones uno y tres

      B<-(svd(desplegar(tensorial(tensorial(T,1,A),3,C1),2),nu=q,nv=0))$u

      # calcula la matriz optima fijando los vectores para las dimensiones uno y dos

      C<-(svd(desplegar(tensorial(tensorial(T,1,A),2,B),3),nu=r,nv=0))$u

      # calcula el tensor core

      G<-tensorial(tensorial(tensorial(T,1,A),2,B),3,C)

      # calcula la norma de los vectores diferencias entre una iteracion y las dos anteriores

      normA1<-sum((A-A1)^2)
      normA2<-sum((A-A2)^2)
      normB1<-sum((B-B1)^2)
      normB2<-sum((B-B2)^2)
      normC1<-sum((C-C1)^2)
      normC2<-sum((C-C2)^2)
      normG1<-sum((G-G1)^2)
      normG2<-sum((G-G2)^2)

      # si los pares de normas anteriores son menores que la tolerancia definida desde un principio se detienen las iteraciones

      if(((normA1<=tol)||(normA2<=tol))&&((normB1<=tol)||(normB2<=tol))&&((normC1<=tol)||(normC2<=tol))&&
         ((normG1<=tol)||(normG2<=tol)))
      {
        seguir<-FALSE
      }
      rep<-rep+1
      return(list(T,p,q,r,A,B,C,G,A1,B1,C1,G1,tol,seguir,rep,iter))
    } else {return(l)}
  }
  l<-Reduce(function(i,...){return(i(...))},rep.int(list(i),iter),
            list(T,p,q,r,A1,B1,C1,G1,A2,B2,C2,G2,tol,TRUE,0,iter),right=TRUE)
  p<-l[[2]]
  q<-l[[3]]
  r<-l[[4]]
  A<-l[[5]]
  B<-l[[6]]
  C<-l[[7]]
  G<-l[[8]]
  rep<-l[[15]]

  # calcula el error: la suma de los cuadrados de las diferencias entre los tensores original y aproximacion

  E<-sum((T-tensorial(tensorial(tensorial(G,1,t(A)),2,t(B)),3,t(C)))^2)

  # finalmente, la funcion devuelve la suma de las componentes, el error, el ajuste en forma de porcentaje y
  # el numero de iteraciones llevadas a cabo
  # cuando se haya elegido la combinacion de componentes devuelve las matrices A,B,C y el tensor core

  if(!mas){return(c(p+q+r,E,100-((100*E)/(sum(T^2))),rep))}else{return(list(A,B,C,G))}
}

# funcion para crear los graficos para el Tucker3, divididos en tres para representar cada una de las tres dimensiones,
# incluso si en alguna se ha retenido solo una componente

plotmt<-function(dim1,dim2,dim3,sos,M1,M2,M3,lim1,lim2,lim3,names1,names2,names3,colores1,colores2,
                 titles=FALSE)
{
  sos<-round(sos,3)
  dev.new(width=14,height=4.666,noRStudioGD=TRUE)
  layout(matrix(1:3,nrow=1))
  layout.show(3)
  if(dim1[2]!=1)
  {
    plot(M1[,1],M1[,2],type="n",xlim=c(min(0,lim1[,1]),max(0,lim1[,1])),
         ylim=c(min(0,lim1[,2]),max(0,lim1[,2])),
         xlab=paste("Component ",dim1[1]," (",sos[dim1[1],1],"%)",sep=""),
         ylab=paste("Component ",dim1[2]," (",sos[dim1[2],1],"%)",sep=""))
    text(M1[,1],M1[,2],names1,cex=0.8,col=colores1)
    abline(v=0,lwd=2)
  } else {
    plot(1:(dim(M1)[1]),M1[,1],type="n",xlim=c(0,dim(M1)[1]+1),
         ylim=c(min(0,M1[,1]),max(0,M1[,1])),xlab="",ylab=paste("Component 1 (",sos[1,1],"%)",sep=""))
    text(1:(dim(M1)[1]),M1[,1],names1,cex=0.8,col=colores1)
  }
  box(lwd=2)
  abline(h=0,lwd=2)
  if(titles){title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)}
  if(dim2[2]!=1)
  {
    plot(M2[,1],M2[,2],type="n",xlim=c(min(0,lim2[,1]),max(0,lim2[,1])),
         ylim=c(min(0,lim2[,2]),max(0,lim2[,2])),
         xlab=paste("Component ",dim2[1]," (",sos[dim2[1],2],"%)",sep=""),
         ylab=paste("Component ",dim2[2]," (",sos[dim2[2],2],"%)",sep=""))
    text(M2[,1],M2[,2],names2,cex=0.8,col=colores2)
    abline(v=0,lwd=2)
    arrows(0,0,M2[,1],M2[,2],col=colores2,angle=10,length=0.1)
  } else {
    plot(1:(dim(M2)[1]),M2[,1],type="n",xlim=c(0,dim(M2)[1]+1),
         ylim=c(min(0,M2[,1]),max(0,M2[,1])),xlab="",ylab=paste("Component 1 (",sos[1,2],"%)",sep=""))
    text(1:(dim(M2)[1]),M2[,1],names2,cex=0.8,col=colores2)
  }
  box(lwd=2)
  abline(h=0,lwd=2)
  if(titles){title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)}
  if(dim3[2]!=1)
  {
    plot(M3[,1],M3[,2],type="n",xlim=c(min(0,lim3[,1]),max(0,lim3[,1])),
         ylim=c(min(0,lim3[,2]),max(0,lim3[,2])),
         xlab=paste("Component ",dim3[1]," (",sos[dim3[1],3],"%)",sep=""),
         ylab=paste("Component ",dim3[2]," (",sos[dim3[2],3],"%)",sep=""))
    text(M3[,1],M3[,2],names3,cex=0.8)
    abline(v=0,lwd=2)
    arrows(0,0,M3[,1],M3[,2],angle=10,length=0.1)
  } else {
    plot(1:(dim(M3)[1]),M3[,1],type="n",xlim=c(0,dim(M3)[1]+1),ylim=c(min(0,M3[,1]),max(0,M3[,1])),
         xlab="",ylab=paste("Component 1 (",sos[1,3],"%)",sep=""))
    text(1:(dim(M3)[1]),M3[,1],names3,cex=0.8)
  }
  box(lwd=2)
  abline(h=0,lwd=2)
  if(titles){title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)}
}

# inicio del programa para realizar un PCA

PCA<-function(X,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc=NULL,norm=FALSE,contr=FALSE)
{

  # lee la matriz de datos con las etiquetas de las filas y las columnas (si no estan incluidas, se nombran por
  # defecto) y suprime las filas que tengan datos faltantes

  l<-read(X)
  X<-l[[1]]
  filas<-l[[2]]
  columnas<-l[[3]]
  namesf<-l[[5]]
  namesc<-l[[6]]

  # lee las etiquetas de los colores para las filas y columnas y comprueba que hay tantas como filas y columnas
  # (si no se dan, se asignan por defecto de color negro)

  l<-colores(filas,columnas,coloresf,coloresc,FALSE)
  coloresf<-l[[1]]
  coloresc<-l[[2]]
  conf<-l[[3]]

  # centra la matriz de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm)

    # crea la matriz con los pesos uniformes para las filas y su raiz cuadrada

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    Dn2<-l[[2]]

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares segun el metodo explicado

    c<-svd(Dn2%*%X)
    d<-c$d
    u<-c$u
    v<-c$v

    # calcula las coordenadas para las filas y las columnas mediante el diagrama de dualidad

    F<-X%*%v%*%diag(1/sqrt(diag(t(v)%*%v)))
    C<-t(X)%*%Dn%*%(solve(Dn2)%*%u%*%diag(1/sqrt(diag(t(u)%*%u))))

    # si se ha elegido, calcula las contribuciones para las filas y las columnas

    if(contr)
    {
      contributions(F,namesf,"las filas",TRUE)
      contributions(C,namesc,"las columnas")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(4.375,2.625),heights=c(3.5,3.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma de la matriz de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,1.25),ylim=c(-1.25,1.5),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      text(0,0,"X")
      text(-1.25,0,filas)
      text(0,1.25,columnas)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        Fd<-F[,dimension]
        Cd<-C[,dimension]

        # representa el grafico para las filas en la izquierda con las etiquetas de colores
        # representa las columnas en la derecha con las etiquetas y vectores desde el origen de colores

        plotm(dim=dimension,d=d,M1=Fd,M2=Cd,lim1=Fd,lim2=Cd,names1=namesf,names2=namesc,colores1=coloresf,
              colores2=coloresc,contf=contributions2(F,dimension),contc=contributions2(C,dimension))
      }
    }
  }

  # fin del programa para el PCA

}

# inicio del programa para realizar un BGA

BGA<-function(X,gruposf,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc=NULL,norm=FALSE,contr=FALSE)
{

  # lee la matriz de datos con las etiquetas de las filas y las columnas (si no estan incluidas, se nombran por
  # defecto) y suprime las filas que tengan datos faltantes

  l<-read(X)
  X<-l[[1]]
  filas<-l[[2]]
  columnas<-l[[3]]
  namesf<-l[[5]]
  namesc<-l[[6]]

  # lee los grupos para las filas y las etiquetas (mas de uno)

  namesg<-unique(gruposf)
  conf<-FALSE
  if(length(namesg)<2)
  {
    conf<-TRUE
    print("te has confundido, el numero de grupos para las filas tiene que ser mayor que uno")
  }

  # lee las etiquetas de los colores para las filas y columnas y comprueba que hay tantas como filas y columnas
  # (si no se dan, se asignan por defecto de color negro)

  if(!conf)
  {
    l<-colores(filas,columnas,coloresf,coloresc,conf)
    coloresf<-l[[1]]
    coloresc<-l[[2]]
    conf<-l[[3]]
  }

  # centra la matriz de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm)

    # crea la matriz con los pesos para las filas segun los tamanos de los grupos y su raiz cuadrada

    l<-pesos(unlist(lapply(1:(length(namesg)),function(x,g,n){return(sum(g==n[x]))},g=gruposf,n=namesg)))
    Dg<-l[[1]]
    Dg2<-l[[2]]

    # calcula la matriz con las medias para los grupos

    XB<-matrix(unlist(lapply(1:length(namesg),function(x,X1,g1,n){return(apply(X1[g1==n[x],],2,mean))},X1=X,
                             g1=gruposf,n=namesg)),ncol=columnas,byrow=TRUE)

    # centra la matriz de las medias por columnas y si se introdujo TRUE en normalizacion normaliza la matriz de medias
    # por columnas

    XB<-preproc(XB,norm)

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares segun el metodo explicado

    c<-svd(Dg2%*%XB)
    d<-c$d
    u<-c$u
    v<-c$v
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para los grupos, las filas y las columnas mediante el diagrama de dualidad

    FB<-XB%*%Vr
    F<-X%*%Vr
    C<-t(XB)%*%Dg%*%solve(Dg2)%*%u%*%diag(1/sqrt(diag(t(u)%*%u)))

    # si se ha elegido, calcula las contribuciones para los grupos, las filas y las columnas

    if(contr)
    {
      contributions(FB,namesg,"los grupos",TRUE)
      contributions(F,namesf,"las filas")
      contributions(C,namesc,"las columnas")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(4.375,2.625),heights=c(3.5,3.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma de las matrices de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,1.25),ylim=c(-2.25,1.5),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-2,1,-1)
      text(0,0,"X")
      text(0,-1.5,"XB")
      text(-1.25,0,filas)
      text(0,1.25,columnas)
      text(-1.25,-1.5,length(namesg))

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        FBd<-FB[,dimension]
        Fd<-F[,dimension]
        Cd<-C[,dimension]
        F2<-rbind(FBd,Fd)

        # representa el grafico para los grupos en la izquierda con las etiquetas
        # representa las columnas en la derecha con las etiquetas y vectores desde el origen segun los colores de los grupos
        # a los que pertenecen

        plotm(dim=dimension,d=d,M1=FBd,M2=Cd,lim1=F2,lim2=Cd,names1=namesg,names2=namesc,
              colores1=rep("black",length(namesg)),colores2=coloresc,contf=contributions2(FB,dimension),
              contc=contributions2(C,dimension))

        # representa el grafico para las filas en la izquierda con las etiquetas segun los colores de los grupos a los que
        # pertenecen
        # representa las columnas en la derecha con las etiquetas y vectores desde el origen segun los colores de los grupos
        # a los que pertenecen

        plotm(dim=dimension,d=d,M1=Fd,M2=Cd,lim1=F2,lim2=Cd,names1=namesf,names2=namesc,colores1=coloresf,
              colores2=coloresc,contf=contributions2(F,dimension),contc=contributions2(C,dimension))
      }
    }
  }

  # fin del programa para el BGA

}

# inicio del programa para realizar un COIA

COIA<-function(X,Y,dimXx=NULL,dimYx=NULL,dimXy=NULL,dimYy=NULL,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc1=NULL,
               coloresc2=NULL,norm=FALSE,contr=FALSE,cotucker=FALSE,tcotucker=FALSE)
{

  # lee las dos matrices de datos con las etiquetas de las filas y de las columnas de ambas matrices (si no estan
  # incluidas, se nombran por defecto), suprime las filas que tengan datos faltantes y comprueba que las dos matrices
  # tengan las mismas filas

  l<-read(X,Y)
  X<-l[[1]]
  filas<-l[[2]]
  columnas1<-l[[3]]
  namesf<-l[[5]]
  namesc1<-l[[6]]
  Y<-l[[8]]
  columnas2<-l[[9]]
  namesc2<-l[[10]]
  conf<-l[[11]]

  # lee los colores de las filas y de las columnas de las dos matrices y comprueba que hay tantos como filas y como
  # columnas (si no se dan, se asignan por defecto en negro)

  if(!conf)
  {
    l<-colores(filas,columnas1,coloresf,coloresc1,conf,columnas2,coloresc2)
    coloresf<-l[[1]]
    coloresc1<-l[[2]]
    coloresc2<-l[[3]]
    conf<-l[[4]]
  }

  # centra las dos matrices de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm)
    Y<-preproc(Y,norm)

    # crea la matriz con los pesos uniformes para las filas y su raiz cuadrada

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    Dn2<-l[[2]]

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares para las dos matrices y para
    # la coinercia segun el metodo explicado

    cX<-svd(Dn2%*%X)
    dx<-cX$d
    ux<-cX$u
    vx<-cX$v
    cY<-svd(Dn2%*%Y)
    dy<-cY$d
    uy<-cY$u
    vy<-cY$v
    c<-svd(t(Y)%*%Dn%*%X)
    d<-c$d
    u<-c$u
    v<-c$v
    Ur<-u%*%diag(1/sqrt(diag(t(u)%*%u)))
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas y las columnas por separado y para la coinercia mediante los diagramas de
    # dualidad

    FX<-X%*%vx%*%diag(1/sqrt(diag(t(vx)%*%vx)))
    CX<-t(X)%*%Dn%*%solve(Dn2)%*%(ux)%*%diag(1/sqrt(diag(t(ux)%*%ux)))
    FY<-Y%*%vy%*%diag(1/sqrt(diag(t(vy)%*%vy)))
    CY<-t(Y)%*%Dn%*%solve(Dn2)%*%uy%*%diag(1/sqrt(diag(t(uy)%*%uy)))
    FXc<-X%*%Vr
    CXc<-t(X)%*%Dn%*%Y%*%Ur
    FYc<-Y%*%Ur
    CYc<-t(Y)%*%Dn%*%X%*%Vr

    # si se ha elegido, calcula las contribuciones para las filas y las columnas de ambas matrices por separado y para la
    # coinercia

    if(contr)
    {
      if(!cotucker)
      {
        contributions(FX,namesf,"las filas segun la primera matriz",TRUE)
        contributions(CX,namesc1,"las columnas de la primera matriz")
        contributions(FY,namesf,"las filas segun la segunda matriz")
        contributions(CY,namesc2,"las columnas de la segunda matriz")
      }
      if(!cotucker||!tcotucker){contributions(FXc,namesf,"las filas segun la primera matriz en la co-inercia",cotucker)}
      if(!cotucker||tcotucker){contributions(CXc,namesc1,"las columnas de la primera matriz en la co-inercia",tcotucker)}
      if(!cotucker||!tcotucker){contributions(FYc,namesf,"las filas segun la segunda matriz en la co-inercia")}
      if(!cotucker||tcotucker){contributions(CYc,namesc2,"las columnas de la segunda matriz en la co-inercia")}
    }

    # si no se ha elegido representar ninguno de los ejes, representa los valores singulares en tres diagramas de
    # sedimentacion

    if(((is.null(dimXx))||(is.null(dimYx)))&&((is.null(dimXy))||(is.null(dimYy)))&&((is.null(dimX))||
                                                                                    (is.null(dimY))))
    {
      if(!cotucker)
      {
        dev.new()
        layout(matrix(c(1,1,2,1),ncol=2),widths=c(4.375,2.625),heights=c(3.5,3.5))
        layout.show(2)
        screeplot(dx)

        # en cada uno de los diagramas representa un esquema con la forma de las matrices de datos correspondientes

        plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,1.25),ylim=c(-1.25,1.5),bty="n",xaxt="n",yaxt="n")
        rect(-1,-1,1,1)
        text(0,0,"X")
        text(-1.25,0,filas)
        text(0,1.25,columnas1)
        dev.new()
        layout(matrix(c(1,1,2,1),ncol=2),widths=c(4.375,2.625),heights=c(3.5,3.5))
        layout.show(2)
        screeplot(dy)
        plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,1.25),ylim=c(-1.25,1.5),bty="n",xaxt="n",yaxt="n")
        rect(-1,-1,1,1)
        text(0,0,"Y")
        text(-1.25,0,filas)
        text(0,1.25,columnas2)
      }
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(35/9,28/9),heights=c(3.5,3.5))
      layout.show(2)
      screeplot(d)
      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,3.25),ylim=c(-4.25,1.5),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(1,-1,3,1)
      rect(0,-4,2,-2)
      text(0,0,"X")
      text(2,0,"Y")
      text(1,-3,"Y' Dn X")
      text(-1.25,0,filas)
      text(0,1.25,columnas1)
      text(2,1.25,columnas2)
      text(-0.25,-3,columnas2)
      text(1,-1.75,columnas1)
    }

    # si se ha elegido representar los ejes para la primera matriz, comprueba que los ejes a representar en los graficos
    # sean correctos

    if(!((is.null(dimXx))||(is.null(dimYx))))
    {
      if(compr(dimXx,dimYx,dx))
      {
        dimensionx<-c(dimXx,dimYx)
        FXd<-FX[,dimensionx]
        CXd<-CX[,dimensionx]

        # representa el grafico para las filas de la primera matriz en la izquierda con las etiquetas segun los colores de
        # los grupos a los que pertenecen
        # representa las columnas de la primera matriz en la derecha con las etiquetas y vectores desde el origen con su
        # color

        plotm(dim=dimensionx,d=dx,M1=FXd,M2=CXd,lim1=FXd,lim2=CXd,names1=namesf,names2=namesc1,colores1=coloresf,
              colores2=coloresc1,contf=contributions2(FX,dimensionx),contc=contributions2(CX,dimensionx))
      }
    }
  }

  # si se ha elegido representar los ejes para la segunda matriz, comprueba que los ejes a representar en los graficos
  # sean correctos

  if(!((is.null(dimXy))||(is.null(dimYy))))
  {
    if(compr(dimXy,dimYy,dy))
    {
      dimensiony<-c(dimXy,dimYy)
      FYd<-FY[,dimensiony]
      CYd<-CY[,dimensiony]

      # representa el grafico de las filas segun la segunda matriz en la izquierda con las etiquetas segun los colores de
      # los grupos a los que pertenecen
      # representa las columnas de la segunda matriz en la derecha con las etiquetas y vectores desde el origen con su
      # color

      plotm(dim=dimensiony,d=dy,M1=FYd,M2=CYd,lim1=FYd,lim2=CYd,names1=namesf,names2=namesc2,colores1=coloresf,
            colores2=coloresc2,contf=contributions2(FY,dimensiony),contc=contributions2(CY,dimensiony))
    }
  }

  # si se ha elegido representar los ejes para la coinercia, comprueba que los ejes a representar en los graficos sean
  # correctos

  if(!((is.null(dimX))||(is.null(dimY))))
  {
    if(compr(dimX,dimY,d))
    {
      dimension<-c(dimX,dimY)
      FXcd<-FXc[,dimension]
      CXcd<-CXc[,dimension]
      FYcd<-FYc[,dimension]
      CYcd<-CYc[,dimension]
      F<-rbind(FXcd,FYcd)

      # representa el grafico para la coinercia de las filas en la izquierda con las etiquetas y vectores de la primera a
      # la segunda matriz segun los colores de los grupos a los que pertenecen
      # representa las columnas de la primera matriz en la derecha con las etiquetas y vectores desde el origen con su
      # color

      plotm(dim=dimension,d=d,M1=FXcd,M2=CXcd,M3=FYcd,lim1=F,lim2=CXcd,names1=namesf,names2=namesc1,
            colores1=coloresf,colores2=coloresc1,contf=contributions2(FXc,dimension),
            contc=contributions2(CXc,dimension),cotucker=cotucker)

      # representa el grafico para la coinercia de las filas en la izquierda con las etiquetas y vectores de la segunda a
      # la primera matriz segun los colores de los grupos a los que pertenecen
      # representa las columnas de la segunda matriz en la derecha con las etiquetas y vectores desde el origen con su
      # color

      plotm(dim=dimension,d=d,M1=FYcd,M2=CYcd,M3=FXcd,lim1=F,lim2=CYcd,names1=namesf,names2=namesc2,
            colores1=coloresf,colores2=coloresc2,contf=contributions2(FYc,dimension),
            contc=contributions2(CYc,dimension),cotucker=cotucker)
    }
  }

  # fin del programa para el COIA

}

# inicio del programa para realizar un PTA

PTA<-function(X,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc=NULL,norm=FALSE,contr=FALSE)
{

  # lee el cubo de datos con las etiquetas de las filas, las columnas y las repeticiones (si no estan incluidas, se
  # nombran por defecto) y suprime las filas que tengan datos faltantes

  l<-read(X)
  X<-l[[1]]
  filas<-l[[2]]
  columnas<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc<-l[[6]]
  namesr<-l[[7]]
  conf<-l[[11]]

  # lee las etiquetas de los colores para las filas y columnas y comprueba que hay tantas como filas y columnas
  # (si no se dan, se asignan por defecto de color negro)

  if(!conf)
  {
    l<-colores(filas,columnas,coloresf,coloresc,conf)
    coloresf<-l[[1]]
    coloresc<-l[[2]]
    conf<-l[[3]]
  }

  # centra el cubo de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm,TRUE)

    # crea la matriz con los pesos uniformes para las filas y las repeticiones y sus raices cuadradas

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    Dn2<-l[[2]]
    l<-pesos(rep(1,repeticiones))
    Dk<-l[[1]]
    Dk2<-l[[2]]

    # calcula la matriz de varianzas-covarianzas vectoriales

    l<-as.matrix(expand.grid(1:repeticiones,1:repeticiones))
    Covv<-matrix(apply(l,1,function(v,X1,D){return(sum(diag(t(X1[,,v[1]])%*%D%*%X1[,,v[2]])))},X1=X,D=Dn),
                 nrow=repeticiones,ncol=repeticiones)

    # extrae los vectores propios segun el metodo explicado y calcula las coordenadas para la interestructura

    a<-(eigen(Covv%*%Dk,symmetric=TRUE))$vectors
    VI<-solve(Dk2)%*%(a[,1:2])%*%diag(1/sqrt(diag(t(a[,1:2])%*%(a[,1:2]))))
    VI[,1]<-abs(VI[,1])
    I<-Covv%*%Dk%*%VI

    # calcula la matriz compromiso

    Xc<-apply(X,1:2,function(v,v2){return(sum(v*v2))},v2=VI[,1]/sum(VI[,1]))

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares segun el metodo explicado

    c<-svd(Dn2%*%Xc)
    d<-c$d
    u<-c$u
    v<-c$v
    Ur<-(solve(Dn2))%*%u%*%diag(1/sqrt(diag(t(u)%*%u)))
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas y las columnas de la matriz compromiso y de las trayectorias mediante el
    # diagrama de dualidad

    Fc<-Xc%*%Vr
    Cc<-t(Xc)%*%Dn%*%Ur
    Ft<-array(apply(X,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(filas,columnas,repeticiones))
    Ct<-array(apply(X,3,function(m,m2){return(t(m)%*%m2)},m2=Dn%*%Ur),dim=c(columnas,columnas,repeticiones))
    F<-matrix(aperm(Ft,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas,byrow=TRUE)
    C<-matrix(aperm(Ct,c(2,1,3)),nrow=columnas*repeticiones,ncol=columnas,byrow=TRUE)

    # si se ha elegido, calcula las contribuciones para las filas y las columnas de la matriz compromiso y de las
    # trayectorias

    if(contr)
    {
      contributions(Fc,namesf,"las filas en el compromiso",TRUE)
      contributions(Cc,namesc,"las columnas en el compromiso")
      contributions(F,rep(namesf,repeticiones),"las filas en todas las repeticiones")
      contributions(C,rep(namesc,repeticiones),"las columnas en todas las repeticiones")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(31.5/9,24.5/9),heights=c(4.5,2.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma del cubo de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-4.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-4,1,-2)
      segments(-1,1,0,2)
      segments(0,2,2,2)
      segments(2,2,2,0)
      segments(2,0,1,-1)
      segments(1,1,2,2)
      text(0,0,"X")
      text(0,-3,"Xc")
      text(-1.25,0,filas)
      text(0,1.25,columnas)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-3,filas)
      text(0,-1.75,columnas)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        Fcd<-Fc[,dimension]
        Ccd<-Cc[,dimension]
        Ftd<-Ft[,dimension,]
        Ctd<-Ct[,dimension,]
        Fd<-F[,dimension]
        Cd<-C[,dimension]
        F2<-rbind(Fcd,Fd)
        C2<-rbind(Ccd,Cd)

        # representa el grafico de la interestructura con las etiquetas para las repeticiones y vectores desde el origen

        dev.new()
        plot(I[,1],I[,2],type="n",xlim=c(min(0,I[,1]),max(0,I[,1])),ylim=c(min(0,I[,2]),max(0,I[,2])),
             xlab="",ylab="")
        text(I[,1],I[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,I[,1],I[,2],angle=10,length=0.1)

        # representa el grafico para las filas del compromiso en la izquierda con las etiquetas segun los colores de los
        # grupos a los que pertenecen
        # representa las columnas del compromiso en la derecha con las etiquetas y vectores desde el origen segun los colores
        # de los grupos a los que pertenecen

        plotm(dim=dimension,d=d,M1=Fcd,M2=Ccd,lim1=F2,lim2=C2,names1=namesf,names2=namesc,colores1=coloresf,
              colores2=coloresc,contf=contributions2(Fc,dimension),contc=contributions2(Cc,dimension))

        # representa el grafico para las trayectorias de las filas en la izquierda con las etiquetas y segmentos segun los
        # colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas en la derecha con las etiquetas, segmentos y vectores desde el origen
        # segun los colores de los grupos a los que pertenecen

        plotm(dim=dimension,d=d,M1=Ftd[,,1],M2=Ctd[,,1],M4=Ftd,M5=Ctd,lim1=F2,lim2=C2,names1=namesf,
              names2=namesc,colores1=coloresf,colores2=coloresc,contf=contributions2(Ft[,,1],dimension),
              contc=contributions2(Ct[,,1],dimension))
      }
    }
  }

  # fin del programa para el PTA

}

# inicio del programa para realizar un STATIS

STATIS<-function(X,dimX=NULL,dimY=NULL,coloresf=NULL,norm=FALSE,contr=FALSE)
{

  # comprueba que todas las matrices del cubo tengan el mismo numero de filas

  conf<-FALSE
  if(length(unique(unlist(lapply(X,function(m){return(dim(m)[1])}))))!=1)
  {
    conf<-TRUE
    print("te has confundido, el numero de filas no es el mismo para todas las matrices")
  }

  # lee el cubo de datos con las etiquetas de las filas y las repeticiones (si no estan incluidas, se
  # nombran por defecto) y suprime las filas que tengan datos faltantes

  if(!conf)
  {
    filas<-dim(X[[1]])[1]
    repeticiones<-length(X)
    if(!is.null(dimnames(X[[1]])[[1]])){namesf<-dimnames(X[[1]])[[1]]}else{namesf<-paste(1:filas,sep="")}
    if(!is.null(names(X))){namesr<-names(X)}else{namesr<-paste("R",1:repeticiones,sep="")}
    nas<-Reduce("&",lapply(lapply(X,is.na),function(m){return(apply(m,1,function(v){return(all(!v))}))}))
    if(!all(nas))
    {
      X<-lapply(X,function(m){return(m[nas,])})
      namesf<-namesf[nas]
      print(paste("se han suprimido algunas filas (",filas-sum(nas),") con datos faltantes",sep=""))
      filas<-sum(nas)
    }

    # lee las etiquetas de los colores para las filas y comprueba que hay tantas como filas
    # (si no se dan, se asignan por defecto de color negro)

    if(is.null(coloresf)){coloresf<-rep("black",filas)}
    if((length(coloresf))!=filas)
    {
      conf<-TRUE
      print("te has confundido, el numero de etiquetas de colores de las filas es distinto del numero de filas")
    }
  }

  # centra cada matriz de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-lapply(X,function(m){return(preproc(m,norm))})
    lapply(X,function(m){m[is.nan(m)]<-0})

    # crea la matriz con los pesos uniformes para las filas y las repeticiones y sus raices cuadradas

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    Dn2<-l[[2]]
    l<-pesos(rep(1,repeticiones))
    Dk<-l[[1]]
    Dk2<-l[[2]]

    # calcula el cubo con las matrices de productos cruzados

    X<-array(unlist(lapply(X,function(m){return(m%*%(t(m)))})),dim=c(filas,filas,repeticiones))

    # calcula la matriz de varianzas-covarianzas vectoriales

    l<-as.matrix(expand.grid(1:repeticiones,1:repeticiones))
    Covv<-matrix(apply(l,1,function(v,X1,D){return(sum(diag(t(X1[,,v[1]])%*%D%*%X1[,,v[2]])))},X1=X,D=Dn),
                 nrow=repeticiones,ncol=repeticiones)

    # extrae los vectores propios segun el metodo explicado y calcula las coordenadas para la interestructura

    a<-(eigen(Covv%*%Dk,symmetric=TRUE))$vectors
    VI<-solve(Dk2)%*%(a[,1:2])%*%diag(1/sqrt(diag(t(a[,1:2])%*%(a[,1:2]))))
    VI[,1]<-abs(VI[,1])
    I<-Covv%*%Dk%*%VI

    # calcula la matriz compromiso

    Xc<-apply(X,1:2,function(v,v2){return(sum(v*v2))},v2=VI[,1]/sum(VI[,1]))

    # extrae los vectores singulares por la derecha y los valores singulares segun el metodo explicado

    c<-svd(Dn2%*%Xc)
    d<-c$d
    v<-c$v
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas de la matriz compromiso y de las trayectorias mediante el
    # diagrama de dualidad

    Fc<-Xc%*%Vr
    Ft<-array(apply(X,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(filas,filas,repeticiones))
    F<-matrix(aperm(Ft,c(2,1,3)),nrow=filas*repeticiones,ncol=filas,byrow=TRUE)

    # si se ha elegido, calcula las contribuciones para las filas de la matriz compromiso y de las
    # trayectorias

    if(contr)
    {
      contributions(Fc,namesf,"las filas en el compromiso",TRUE)
      contributions(F,rep(namesf,repeticiones),"las filas en todas las repeticiones")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(31.5/9,24.5/9),heights=c(4.5,2.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma del cubo de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-4.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-4,1,-2)
      segments(-1,1,0,2)
      segments(0,2,2,2)
      segments(2,2,2,0)
      segments(2,0,1,-1)
      segments(1,1,2,2)
      text(0,0,"X")
      text(0,-3,"Xc")
      text(-1.25,0,filas)
      text(0,1.25,filas)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-3,filas)
      text(0,-1.75,filas)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        Fcd<-Fc[,dimension]
        Ftd<-Ft[,dimension,]
        Fd<-F[,dimension]
        F2<-rbind(Fcd,Fd)

        # representa el grafico de la interestructura con las etiquetas para las repeticiones y vectores desde el origen

        dev.new()
        plot(I[,1],I[,2],type="n",xlim=c(min(0,I[,1]),max(0,I[,1])),ylim=c(min(0,I[,2]),max(0,I[,2])),
             xlab="",ylab="")
        text(I[,1],I[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,I[,1],I[,2],angle=10,length=0.1)

        # representa el grafico para las filas del compromiso con las etiquetas segun los colores de los
        # grupos a los que pertenecen

        contf<-contributions2(Fc,dimension)
        contf[is.nan(contf)]<-0
        mcontf<-max(contf)/2
        dev.new()
        e<-round((100*d^2/sum(d^2))[dimension],3)
        plot(Fcd[,1],Fcd[,2],type="n",xlim=c(min(0,F2[,1]),max(0,F2[,1])),ylim=c(min(0,F2[,2]),max(0,F2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        text(Fcd[,1],Fcd[,2],namesf,cex=0.8,col=coloresf)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        dev.new()
        plot(Fcd[,1],Fcd[,2],type="n",xlim=c(min(0,F2[,1]),max(0,F2[,1])),
             ylim=c(min(0,F2[,2]),max(0,F2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        lcoloresf<-coloresf[contf>=mcontf]
        lnamesf<-namesf[contf>=mcontf]
        lFcd<-matrix(Fcd[contf>=mcontf,],ncol=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)
        points(lFcd[,1],lFcd[,2],pch=19,col=lcoloresf,cex=0.75)
        identify(lFcd[,1],lFcd[,2],lnamesf,cex=0.8,tolerance=3)

        # representa el grafico para las trayectorias de las filas con las etiquetas y segmentos segun los
        # colores de los grupos a los que pertenecen

        contf<-contributions2(Ft[,,1],dimension)
        contf[is.nan(contf)]<-0
        mcontf<-max(contf)/2
        s<-function(i,M,colores){segments(M[,1,i],M[,2,i],M[,1,i+1],M[,2,i+1],col=colores)}
        dev.new()
        plot(Ftd[,1,1],Ftd[,2,1],type="n",xlim=c(min(0,F2[,1]),max(0,F2[,1])),ylim=c(min(0,F2[,2]),max(0,F2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        text(Ftd[,1,1],Ftd[,2,1],namesf,cex=0.8,col=coloresf)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        lapply(1:(dim(Ftd)[3]-1),s,M=Ftd,colores=coloresf)
        dev.new()
        plot(Ftd[,1,1],Ftd[,2,1],type="n",xlim=c(min(0,F2[,1]),max(0,F2[,1])),
             ylim=c(min(0,F2[,2]),max(0,F2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        lcoloresf<-coloresf[contf>=mcontf]
        lnamesf<-namesf[contf>=mcontf]
        lM1<-matrix(Ftd[contf>=mcontf,,1],ncol=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)
        lM4<-array(Ftd[contf>=mcontf,,],dim=c(length(lcoloresf),2,dim(Ftd)[3]))
        lapply(1:(dim(Ftd)[3]-1),s,M=lM4,colores=lcoloresf)
        points(lM1[,1],lM1[,2],pch=19,col=lcoloresf,cex=0.75)
        identify(lM1[,1],lM1[,2],lnamesf,cex=0.8,tolerance=3)
      }
    }
  }

  # fin del programa para el STATIS

}

# inicio del programa para realizar un STATISDUAL

STATISDUAL<-function(X,dimX=NULL,dimY=NULL,coloresc=NULL,norm=FALSE,contr=FALSE)
{

  # comprueba que todas las matrices del cubo tengan el mismo numero de columnas

  conf<-FALSE
  if(length(unique(unlist(lapply(X,function(m){return(dim(m)[2])}))))!=1)
  {
    conf<-TRUE
    print("te has confundido, el numero de columnas no es el mismo para todas las matrices")
  }

  # lee el cubo de datos con las etiquetas de las columnas y las repeticiones (si no estan incluidas, se
  # nombran por defecto) y suprime las filas que tengan datos faltantes

  if(!conf)
  {
    columnas<-dim(X[[1]])[2]
    repeticiones<-length(X)
    if(!is.null(dimnames(X[[1]])[[2]])){namesc<-dimnames(X[[1]])[[2]]}else{namesc<-paste(1:columnas,sep="")}
    if(!is.null(names(X))){namesr<-names(X)}else{namesr<-paste("R",1:repeticiones,sep="")}
    l<-lapply(X,function(m){
      nas<-apply(is.na(m),1,function(v){return(all(!v))})
      return(list(m[nas,],dim(m)[1]-sum(nas)))
    })
    X<-lapply(l,function(l2){return(l2[[1]])})
    s<-Reduce(sum,lapply(l,function(l2){return(l2[[2]])}))
    if(s>0){print(paste("se han suprimido algunas filas (",s,") con datos faltantes",sep=""))}

    # lee las etiquetas de los colores para las columnas y comprueba que hay tantas como columnas
    # (si no se dan, se asignan por defecto de color negro)

    if(is.null(coloresc)){coloresc<-rep("black",columnas)}
    conf<-FALSE
    if((length(coloresc))!=columnas)
    {
      conf<-TRUE
      print("te has confundido, el numero de etiquetas de colores de las columnas es distinto del numero de columnas")
    }
  }

  # centra cada matriz de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-lapply(X,function(m){return(preproc(m,norm))})
    lapply(X,function(m){m[is.nan(m)]<-0})

    # crea la matriz con los pesos uniformes para las filas de cada matriz y las repeticiones y sus raices cuadradas

    l<-lapply(X,function(m){return(pesos(rep(1,dim(m)[1])))})
    Dn<-lapply(l,function(l2){return(l2[[1]])})
    l<-pesos(rep(1,repeticiones))
    Dk<-l[[1]]
    Dk2<-l[[2]]

    # calcula el cubo con las matrices de productos cruzados

    X<-array(unlist(lapply(1:(length(X)),function(i,X1,D){return((t(X1[[i]]))%*%D[[i]]%*%X1[[i]])},X1=X,D=Dn)),dim=c(columnas,columnas,repeticiones))

    # calcula la matriz de varianzas-covarianzas vectoriales

    l<-as.matrix(expand.grid(1:repeticiones,1:repeticiones))
    Covv<-matrix(apply(l,1,function(v,X1,D){return(sum(diag(t(X1[,,v[1]])%*%X1[,,v[2]])))},X1=X),
                 nrow=repeticiones,ncol=repeticiones)

    # extrae los vectores propios segun el metodo explicado y calcula las coordenadas para la interestructura

    a<-(eigen(Covv%*%Dk,symmetric=TRUE))$vectors
    VI<-solve(Dk2)%*%(a[,1:2])%*%diag(1/sqrt(diag(t(a[,1:2])%*%(a[,1:2]))))
    VI[,1]<-abs(VI[,1])
    I<-Covv%*%Dk%*%VI

    # calcula la matriz compromiso

    Xc<-apply(X,1:2,function(v,v2){return(sum(v*v2))},v2=VI[,1]/sum(VI[,1]))

    # extrae los vectores singulares por la izquierda y los valores singulares segun el metodo explicado

    c<-svd(Xc)
    d<-c$d
    u<-c$u
    Ur<-u%*%diag(1/sqrt(diag(t(u)%*%u)))

    # calcula las coordenadas para las columnas de la matriz compromiso y de las trayectorias mediante el
    # diagrama de dualidad

    Cc<-t(Xc)%*%Ur
    Ct<-array(apply(X,3,function(m,m2){return(t(m)%*%m2)},m2=Ur),dim=c(columnas,columnas,repeticiones))
    C<-matrix(aperm(Ct,c(2,1,3)),nrow=columnas*repeticiones,ncol=columnas,byrow=TRUE)

    # si se ha elegido, calcula las contribuciones para las columnas de la matriz compromiso y de las
    # trayectorias

    if(contr)
    {
      contributions(Cc,namesc,"las columnas en el compromiso",TRUE)
      contributions(C,rep(namesc,repeticiones),"las columnas en todas las repeticiones")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(31.5/9,24.5/9),heights=c(4.5,2.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma del cubo de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-4.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-4,1,-2)
      segments(-1,1,0,2)
      segments(0,2,2,2)
      segments(2,2,2,0)
      segments(2,0,1,-1)
      segments(1,1,2,2)
      text(0,0,"X")
      text(0,-3,"Xc")
      text(-1.25,0,columnas)
      text(0,1.25,columnas)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-3,columnas)
      text(0,-1.75,columnas)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        Ccd<-Cc[,dimension]
        Ctd<-Ct[,dimension,]
        Cd<-C[,dimension]
        C2<-rbind(Ccd,Cd)

        # representa el grafico de la interestructura con las etiquetas para las repeticiones y vectores desde el origen

        dev.new()
        plot(I[,1],I[,2],type="n",xlim=c(min(0,I[,1]),max(0,I[,1])),ylim=c(min(0,I[,2]),max(0,I[,2])),
             xlab="",ylab="")
        text(I[,1],I[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,I[,1],I[,2],angle=10,length=0.1)

        # representa las columnas del compromiso con las etiquetas y vectores desde el origen segun los colores
        # de los grupos a los que pertenecen

        contc<-contributions2(Cc,dimension)
        contc[is.nan(contc)]<-0
        mcontc<-max(contc)/2
        dev.new()
        e<-round((100*d^2/sum(d^2))[dimension],3)
        plot(Ccd[,1],Ccd[,2],type="n",xlim=c(min(0,C2[,1]),max(0,C2[,1])),ylim=c(min(0,C2[,2]),max(0,C2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        text(Ccd[,1],Ccd[,2],namesc,cex=0.8,col=coloresc)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        lcoloresc<-coloresc[!apply(Ccd==0,1,all)]
        lCcd<-Ccd[!apply(Ccd==0,1,all),]
        arrows(0,0,lCcd[,1],lCcd[,2],col=lcoloresc,angle=10,length=0.1)
        dev.new()
        plot(Ccd[,1],Ccd[,2],type="n",xlim=c(min(0,C2[,1]),max(0,C2[,1])),
             ylim=c(min(0,C2[,2]),max(0,C2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)
        lcoloresc<-coloresc[contc>=mcontc]
        lnamesc<-namesc[contc>=mcontc]
        lCcd<-matrix(Ccd[contc>=mcontc,],ncol=2)
        lcoloresc<-lcoloresc[!apply(lCcd==0,1,all)]
        lnamesc<-lnamesc[!apply(lCcd==0,1,all)]
        lCcd<-matrix(lCcd[!apply(lCcd==0,1,all),],ncol=2)
        text(lCcd[,1],lCcd[,2],lnamesc,cex=0.8,col=lcoloresc)
        arrows(0,0,lCcd[,1],lCcd[,2],col=lcoloresc,angle=10,length=0.1)

        # representa las trayectorias de las columnas con las etiquetas, segmentos y vectores desde el origen
        # segun los colores de los grupos a los que pertenecen

        contc<-contributions2(Ct[,,1],dimension)
        contc[is.nan(contc)]<-0
        mcontc<-max(contc)/2
        s<-function(i,M,colores){segments(M[,1,i],M[,2,i],M[,1,i+1],M[,2,i+1],col=colores)}
        dev.new()
        plot(Ctd[,1,1],Ctd[,2,1],type="n",xlim=c(min(0,C2[,1]),max(0,C2[,1])),ylim=c(min(0,C2[,2]),max(0,C2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        text(Ctd[,1,1],Ctd[,2,1],namesc,cex=0.8,col=coloresc)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        lcoloresc<-coloresc[!apply(Ctd[,,1]==0,1,all)]
        lM2<-Ctd[!apply(Ctd[,,1]==0,1,all),,1]
        arrows(0,0,lM2[,1],lM2[,2],col=lcoloresc,angle=10,length=0.1)
        lapply(1:(dim(Ctd)[3]-1),s,M=Ctd,colores=coloresc)
        dev.new()
        plot(Ctd[,1,1],Ctd[,2,1],type="n",xlim=c(min(0,C2[,1]),max(0,C2[,1])),
             ylim=c(min(0,C2[,2]),max(0,C2[,2])),
             xlab=paste("Axis ",dimension[1]," (",e[1],"%)",sep=""),ylab=paste("Axis ",dimension[2]," (",e[2],"%)",sep=""))
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        title(main="Elementos con mayor contribucion.",cex.main=0.75,font.main=1)
        lcoloresc<-coloresc[contc>=mcontc]
        lnamesc<-namesc[contc>=mcontc]
        lM2<-matrix(Ctd[contc>=mcontc,,1],ncol=2)
        lcoloresc<-lcoloresc[!apply(lM2==0,1,all)]
        lnamesc<-lnamesc[!apply(lM2==0,1,all)]
        lM2<-matrix(lM2[!apply(lM2==0,1,all),],ncol=2)
        text(lM2[,1],lM2[,2],lnamesc,cex=0.8,col=lcoloresc)
        arrows(0,0,lM2[,1],lM2[,2],col=lcoloresc,angle=10,length=0.1)
        lCtd<-array(Ctd[contc>=mcontc,,],dim=c(length(lcoloresc),2,dim(Ctd)[3]))
        lapply(1:(dim(lCtd)[3]-1),s,M=lCtd,colores=lcoloresc)
        return(NULL)
      }
    }
  }

  # fin del programa para el STATIS DUAL

}

# inicio del programa para realizar un BGCOIA

BGCOIA<-function(X,Y,dimXx=NULL,dimYx=NULL,dimXy=NULL,dimYy=NULL,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc1=NULL,
                 coloresc2=NULL,norm=FALSE,contr=FALSE)
{

  # lee los dos cubos de datos con las etiquetas de las filas, las columnas y las repeticiones (si no estan incluidas,
  # se nombran por defecto), comprueba que los dos cubos tengan las mismas filas y repeticiones y suprime las filas que
  # tengan datos faltantes

  l<-read(X,Y)
  X<-l[[1]]
  filas<-l[[2]]
  columnas1<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc1<-l[[6]]
  namesr<-l[[7]]
  Y<-l[[8]]
  columnas2<-l[[9]]
  namesc2<-l[[10]]
  conf<-l[[11]]

  # lee los colores de las filas y de las columnas de los dos cubos y comprueba que hay tantos como filas y como
  # columnas (si no se dan, se asignan por defecto en negro)

  if(!conf)
  {
    l<-colores(filas,columnas1,coloresf,coloresc1,conf,columnas2,coloresc2)
    coloresf<-l[[1]]
    coloresc1<-l[[2]]
    coloresc2<-l[[3]]
    conf<-l[[4]]
  }

  if(!conf)
  {

    # crea la matriz con los pesos uniformes para las repeticiones y su raiz cuadrada

    l<-pesos(rep(1,repeticiones))
    Dg<-l[[1]]
    Dg2<-l[[2]]

    # calcula las matrices con las medias por repeticiones de ambos cubos

    XB<-matrix(unlist(lapply(1:repeticiones,function(x,X1){return(apply(X1[,,x],2,mean))},X1=X)),ncol=columnas1,
               byrow=TRUE)
    YB<-matrix(unlist(lapply(1:repeticiones,function(x,X1){return(apply(X1[,,x],2,mean))},X1=Y)),ncol=columnas2,
               byrow=TRUE)

    # centra las matrices de las medias por columnas y si se introdujo TRUE en normalizacion normaliza las matrices de
    # medias por columnas

    XB<-preproc(XB,norm)
    YB<-preproc(YB,norm)

    # centra los cubo de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

    X<-preproc(X,norm,TRUE)
    Y<-preproc(Y,norm,TRUE)

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares para las dos matrices y para
    # la coinercia segun el metodo explicado

    cX<-svd(Dg2%*%XB)
    dx<-cX$d
    ux<-cX$u
    vx<-cX$v
    cY<-svd(Dg2%*%YB)
    dy<-cY$d
    uy<-cY$u
    vy<-cY$v
    c<-svd(t(YB)%*%Dg%*%XB)
    d<-c$d
    u<-c$u
    v<-c$v
    Ur<-u%*%diag(1/sqrt(diag(t(u)%*%u)))
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas, las columnas y las repeticiones para las matrices de las medias, para la
    # coinercia y para las trayectorias mediante el diagrama de dualidad

    FXB<-XB%*%vx%*%diag(1/sqrt(diag(t(vx)%*%vx)))
    CXB<-t(XB)%*%Dg%*%solve(Dg2)%*%ux%*%diag(1/sqrt(diag(t(ux)%*%ux)))
    FYB<-YB%*%vy%*%diag(1/sqrt(diag(t(vy)%*%vy)))
    CYB<-t(YB)%*%Dg%*%solve(Dg2)%*%uy%*%diag(1/sqrt(diag(t(uy)%*%uy)))
    FXBc<-XB%*%Vr
    CXBc<-t(XB)%*%Dg%*%YB%*%Ur
    FYBc<-YB%*%Ur
    CYBc<-t(YB)%*%Dg%*%XB%*%Vr
    FX<-array(apply(X,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(filas,columnas2,repeticiones))
    FY<-array(apply(Y,3,function(m,m2){return(m%*%m2)},m2=Ur),dim=c(filas,columnas2,repeticiones))
    PFX<-matrix(aperm(FX,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)
    PFY<-matrix(aperm(FY,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)

    # si se ha elegido, calcula las contribuciones para las filas, las columnas y las repeticiones para las matrices
    # de las medias, para la coinercia y para las trayectorias

    if(contr)
    {
      contributions(FXB,namesr,"las repeticiones segun el primer cubo",TRUE)
      contributions(CXB,namesc1,"las columnas del primer cubo")
      contributions(FYB,namesr,"las repeticiones segun el segundo cubo")
      contributions(CYB,namesc2,"las columnas del segundo cubo")
      contributions(FXBc,namesr,"las repeticiones segun primer cubo en la co-inercia")
      contributions(CXBc,namesc1,"las columnas del primer cubo en la co-inercia")
      contributions(FYBc,namesr,"las repeticiones segun el segundo cubo en la co-inercia")
      contributions(CYBc,namesc2,"las columnas del segundo cubo en la co-inercia")
      contributions(PFX,rep(namesf,repeticiones),"las filas segun el primer cubo en todas las repeticiones")
      contributions(PFY,rep(namesf,repeticiones),"las filas segun el segundo cubo en todas las repeticiones")
    }

    # si no se ha elegido representar ninguno de los ejes, representa los valores singulares en tres diagramas de
    # sedimentacion

    if(((is.null(dimXx))||(is.null(dimYx)))&&((is.null(dimXy))||(is.null(dimYy)))&&((is.null(dimX))||
                                                                                    (is.null(dimY))))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(3.9375,3.0625),heights=c(3.5,3.5))
      layout.show(2)
      screeplot(dx)

      # en cada uno de los diagramas representa un esquema con la forma de los cubos de datos correspondientes

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-2.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-2,1,-1)
      segments(-1,1,0,2)
      segments(0,2,2,2)
      segments(2,2,2,0)
      segments(2,0,1,-1)
      segments(1,1,2,2)
      text(0,0,"X")
      text(0,-1.5,"XB")
      text(-1.25,0,filas)
      text(0,1.25,columnas1)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-1.5,repeticiones)
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(3.9375,3.0625),heights=c(3.5,3.5))
      layout.show(2)
      screeplot(dy)
      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-2.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(-1,-2,1,-1)
      segments(-1,1,0,2)
      segments(0,2,2,2)
      segments(2,2,2,0)
      segments(2,0,1,-1)
      segments(1,1,2,2)
      text(0,0,"Y")
      text(0,-1.5,"YB")
      text(-1.25,0,filas)
      text(0,1.25,columnas2)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-1.5,repeticiones)
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(2.85,4.15),heights=c(4.5,2.5))
      layout.show(2)
      screeplot(d)
      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,4.25),ylim=c(-5.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(1,-1,3,1)
      rect(-1,-2,1,-1)
      rect(1,-2,3,1)
      rect(0,-5,2,-3)
      segments(-1,1,0,2)
      segments(0,2,4,2)
      segments(4,2,4,0)
      segments(4,0,3,-1)
      segments(1,1,2,2)
      segments(3,1,4,2)
      text(0,0,"X")
      text(2,0,"Y")
      text(0,-1.5,"XB")
      text(2,-1.5,"YB")
      text(1,-4,"YB' Dg XB")
      text(-1.25,0,filas)
      text(0,1.25,columnas1)
      text(2,1.25,columnas2)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-1.5,repeticiones)
      text(-0.25,-4,columnas2)
      text(1,-2.75,columnas1)
    }

    # si se ha elegido representar los ejes para la primera matriz, comprueba que los ejes a representar en los
    # graficos sean correctos

    if(!((is.null(dimXx))||(is.null(dimYx))))
    {
      if(compr(dimXx,dimYx,dx))
      {
        dimensionx<-c(dimXx,dimYx)
        FXBd<-FXB[,dimensionx]
        CXBd<-CXB[,dimensionx]

        # representa el grafico para las repeticiones del primer cubo en la izquierda con las etiquetas
        # representa las columnas de la primera matriz de medias en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimensionx,d=dx,M1=FXBd,M2=CXBd,lim1=FXBd,lim2=CXBd,names1=namesr,names2=namesc1,
              colores1=rep("black",repeticiones),colores2=coloresc1,contf=contributions2(FXB,dimensionx),
              contc=contributions2(CXB,dimensionx))
      }
    }

    # si se ha elegido representar los ejes para la segunda matriz, comprueba que los ejes a representar en los
    # graficos sean correctos

    if(!((is.null(dimXy))||(is.null(dimYy))))
    {
      if(compr(dimXy,dimYy,dy))
      {
        dimensiony<-c(dimXy,dimYy)
        FYBd<-FYB[,dimensiony]
        CYBd<-CYB[,dimensiony]

        # representa el grafico para las repeticiones del segundo cubo en la izquierda con las etiquetas
        # representa las columnas de la segunda matriz de medias en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimensiony,d=dy,M1=FYBd,M2=CYBd,lim1=FYBd,lim2=CYBd,names1=namesr,names2=namesc2,
              colores1=rep("black",repeticiones),colores2=coloresc2,contf=contributions2(FYB,dimensiony),
              contc=contributions2(CYB,dimensiony))
      }
    }

    # si se ha elegido representar los ejes para la coinercia y las trayectorias, comprueba que los ejes a representar
    # en los graficos sean correctos

    if(!((is.null(dimX))||(is.null(dimY))))
    {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        FXBcd<-FXBc[,dimension]
        CXBcd<-CXBc[,dimension]
        FYBcd<-FYBc[,dimension]
        CYBcd<-CYBc[,dimension]
        FXd<-FX[,dimension,]
        FYd<-FY[,dimension,]
        PFXd<-PFX[,dimension]
        PFYd<-PFY[,dimension]
        Fc<-rbind(FXBcd,FYBcd)

        # representa el grafico para la coinercia de las repeticiones en la izquierda con las etiquetas y vectores de la
        # primera a la segunda matriz de medias segun los colores de los grupos a los que pertenecen
        # representa las columnas de la primera matriz en la derecha con las etiquetas y vectores desde el origen con su
        # color

        plotm(dim=dimension,d=d,M1=FXBcd,M2=CXBcd,M3=FYBcd,lim1=Fc,lim2=CXBcd,names1=namesr,names2=namesc1,
              colores1=rep("black",repeticiones),colores2=coloresc1,contf=contributions2(FXBc,dimension),
              contc=contributions2(CXBc,dimension))

        # representa el grafico para la coinercia de las repeticiones en la izquierda con las etiquetas y vectores de la
        # segunda a la primera matriz de medias segun los colores de los grupos a los que pertenecen
        # representa las columnas de la segunda matriz en la derecha con las etiquetas y vectores desde el origen con su
        # color

        plotm(dim=dimension,d=d,M1=FYBcd,M2=CYBcd,M3=FXBcd,lim1=Fc,lim2=CYBcd,names1=namesr,names2=namesc2,
              colores1=rep("black",repeticiones),colores2=coloresc2,contf=contributions2(FYBc,dimension),
              contc=contributions2(CYBc,dimension))

        # representa el grafico para las trayectorias de las filas segun el primer cubo en la izquierda con las etiquetas
        # y segmentos segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del primer cubo en la derecha con las etiquetas, segmentos y vectores
        # desde el origen segun su color

        plotm(dim=dimension,d=d,M1=FXd[,,1],M2=CXBcd,M4=FXd,lim1=PFXd,lim2=CXBcd,names1=namesf,names2=namesc1,
              colores1=coloresf,colores2=coloresc1,contf=contributions2(FX[,,1],dimension),
              contc=contributions2(CXBc,dimension))

        # representa el grafico para las trayectorias de las filas segun el segundo cubo en la izquierda con las etiquetas
        # y segmentos segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del segundo cubo en la derecha con las etiquetas, segmentos y vectores
        # desde el origen segun su color

        plotm(dim=dimension,d=d,M1=FYd[,,1],M2=CYBcd,M4=FYd,lim1=PFYd,lim2=CYBcd,names1=namesf,names2=namesc2,
              colores1=coloresf,colores2=coloresc2,contf=contributions2(FY[,,1],dimension),
              contc=contributions2(CYBc,dimension))
      }
    }
  }

  # fin del programa para el BGCOIA

}

# inicio del programa para realizar un STATICO

STATICO<-function(X,Y,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc1=NULL,coloresc2=NULL,norm=FALSE,contr=FALSE)
{

  # lee los dos cubo de datos con las etiquetas de las filas, las columnas y las repeticiones (si no estan incluidas,
  # se nombran por defecto), comprueba que los dos cubos tengan las mismas filas y repeticiones y suprime las filas
  # que tengan datos faltantes

  l<-read(X,Y)
  X<-l[[1]]
  filas<-l[[2]]
  columnas1<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc1<-l[[6]]
  namesr<-l[[7]]
  Y<-l[[8]]
  columnas2<-l[[9]]
  namesc2<-l[[10]]
  conf<-l[[11]]

  # lee los colores de las filas y de las columnas de los dos cubos y comprueba que hay tantos como filas y como
  # columnas (si no se dan, se asignan por defecto en negro)

  if(!conf)
  {
    l<-colores(filas,columnas1,coloresf,coloresc1,conf,columnas2,coloresc2)
    coloresf<-l[[1]]
    coloresc1<-l[[2]]
    coloresc2<-l[[3]]
    conf<-l[[4]]
  }

  # centra los cubos de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm,TRUE)
    Y<-preproc(Y,norm,TRUE)

    # crea la matriz con los pesos uniformes para las filas, las repeticiones y la raiz cuadrada de esta

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    l<-pesos(rep(1,repeticiones))
    Dk<-l[[1]]
    Dk2<-l[[2]]

    # crea el cubo con las covarianzas entre los dos originales y calcula la matriz de varianzas-covarianzas vectoriales

    Z<-array(unlist(lapply(1:repeticiones,function(x,A,B,D){return(t(A[,,x])%*%D%*%B[,,x])},A=Y,B=X,D=Dn)),
             dim=c(columnas2,columnas1,repeticiones))
    l<-as.matrix(expand.grid(1:repeticiones,1:repeticiones))
    Covv<-matrix(apply(l,1,function(v,X1){return(sum(diag(t(X1[,,v[1]])%*%X1[,,v[2]])))},X1=Z),nrow=repeticiones,
                 ncol=repeticiones)

    # extrae los vectores propios segun el metodo explicado y calcula las coordenadas para la interestructura

    a<-(eigen(Covv%*%Dk,symmetric=TRUE))$vectors
    VI<-solve(Dk2)%*%(a[,1:2])%*%diag(1/sqrt(diag(t(a[,1:2])%*%(a[,1:2]))))
    VI[,1]<-abs(VI[,1])
    I<-Covv%*%Dk%*%VI

    # calcula la matriz compromiso

    Zc<-apply(Z,1:2,function(v,v2){return(sum(v*v2))},v2=VI[,1]/sum(VI[,1]))

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares segun el metodo explicado

    c<-svd(Zc)
    d<-c$d
    u<-c$u
    v<-c$v
    Ur<-u%*%diag(1/sqrt(diag(t(u)%*%u)))
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas y las columnas del compromiso y de las trayectorias mediante el diagrama
    # de dualidad

    Fc<-Zc%*%Vr
    Cc<-t(Zc)%*%Ur
    FZ<-array(apply(Z,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(columnas2,columnas2,repeticiones))
    CZ<-array(apply(Z,3,function(m,m2){return(t(m)%*%m2)},m2=Ur),dim=c(columnas1,columnas2,repeticiones))
    PFZ<-matrix(aperm(FZ,c(2,1,3)),nrow=columnas2*repeticiones,ncol=columnas2,byrow=TRUE)
    PCZ<-matrix(aperm(CZ,c(2,1,3)),nrow=columnas1*repeticiones,ncol=columnas2,byrow=TRUE)
    FX<-array(apply(X,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(filas,columnas2,repeticiones))
    FY<-array(apply(Y,3,function(m,m2){return(m%*%m2)},m2=Ur),dim=c(filas,columnas2,repeticiones))
    PFX<-matrix(aperm(FX,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)
    PFY<-matrix(aperm(FY,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)
    CZ2<-rbind(Cc,PCZ)
    FZ2<-rbind(Fc,PFZ)

    # si se ha elegido, calcula las contribuciones para las filas y las columnas de la matriz compromiso y de las
    # trayectorias de las filas y las columnas

    if(contr)
    {
      contributions(Fc,namesc2,"las columnas del segundo cubo en el compromiso",TRUE)
      contributions(Cc,namesc1,"las columnas del primer cubo en el compromiso")
      contributions(PFZ,rep(namesc2,repeticiones),"las columnas del segundo cubo en todas las repeticiones")
      contributions(PCZ,rep(namesc1,repeticiones),"las columnas del primer cubo en todas las repeticiones")
      contributions(PFX,rep(namesf,repeticiones),"las filas segun el primer cubo en todas las repeticiones")
      contributions(PFY,rep(namesf,repeticiones),"las filas segun el segundo cubo en todas las repeticiones")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(3.375,3.625),heights=c(5.5,1.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma de los cubos de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,4.25),ylim=c(-8.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(1,-1,3,1)
      rect(0,-5,2,-3)
      rect(0,-8,2,-6)
      segments(-1,1,0,2)
      segments(0,2,4,2)
      segments(4,2,4,0)
      segments(4,0,3,-1)
      segments(1,1,2,2)
      segments(3,1,4,2)
      segments(0,-3,1,-2)
      segments(1,-2,3,-2)
      segments(3,-2,3,-4)
      segments(3,-4,2,-5)
      segments(2,-3,3,-2)
      text(0,0,"X")
      text(2,0,"Y")
      text(1,-4,"Z")
      text(1,-7,"Zc")
      text(-1.25,0,filas)
      text(0,1.25,columnas1)
      text(2,1.25,columnas2)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-0.25,-4,columnas2)
      text(1,-2.75,columnas1)
      text(0.25,-2.25,repeticiones,srt=45)
      text(-0.25,-7,columnas2)
      text(1,-5.75,columnas1)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        Fcd<-Fc[,dimension]
        Ccd<-Cc[,dimension]
        FZd<-FZ[,dimension,]
        CZd<-CZ[,dimension,]
        FXd<-FX[,dimension,]
        FYd<-FY[,dimension,]
        PFXd<-PFX[,dimension]
        PFYd<-PFY[,dimension]
        FZ2d<-FZ2[,dimension]
        CZ2d<-CZ2[,dimension]

        # representa el grafico de la interestructura con las etiquetas para las repeticiones y vectores desde el origen

        dev.new()
        plot(I[,1],I[,2],type="n",xlim=c(min(0,I[,1]),max(0,I[,1])),ylim=c(min(0,I[,2]),max(0,I[,2])),xlab="",
             ylab="")
        text(I[,1],I[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,I[,1],I[,2],angle=10,length=0.1)

        # representa el grafico para las columnas del primer cubo en la izquierda con las etiquetas y vectores desde el
        # origen segun su color
        # representa las columnas del segundo cubo en la derecha con las etiquetas y vectores desde el origen segun su color

        plotm(dim=dimension,d=d,M1=Ccd,M2=Fcd,M6=Ccd,lim1=CZ2d,lim2=FZ2d,names1=namesc1,names2=namesc2,
              colores1=coloresc1,colores2=coloresc2,contf=contributions2(Cc,dimension),
              contc=contributions2(Fc,dimension))

        # representa el grafico para las trayectorias de las columnas del primer cubo en la izquierda con las etiquetas y
        # vectores segun su color
        # representa las trayectorias de las columnas del segundo cubo en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimension,d=d,M1=CZd[,,1],M2=FZd[,,1],M4=CZd,M5=FZd,M6=CZd[,,1],lim1=CZ2d,lim2=FZ2d,
              names1=namesc1,names2=namesc2,colores1=coloresc1,colores2=coloresc2,
              contf=contributions2(CZ[,,1],dimension),contc=contributions2(FZ[,,1],dimension))

        # representa el grafico para las trayectorias de las filas segun el primer cubo en la izquierda con las etiquetas
        # segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del primer cubo en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimension,d=d,M1=FXd[,,1],M2=CZd[,,1],M4=FXd,M5=CZd,lim1=PFXd,lim2=CZ2d,names1=namesf,
              names2=namesc1,colores1=coloresf,colores2=coloresc1,contf=contributions2(FX[,,1],dimension),
              contc=contributions2(CZ[,,1],dimension))

        # representa el grafico para las trayectorias de las filas segun el segundo cubo en la izquierda con las etiquetas
        # segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del segundo cubo en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimension,d=d,M1=FYd[,,1],M2=FZd[,,1],M4=FYd,M5=FZd,lim1=PFYd,lim2=FZ2d,names1=namesf,
              names2=namesc2,colores1=coloresf,colores2=coloresc2,contf=contributions2(FY[,,1],dimension),
              contc=contributions2(FZ[,,1],dimension))
      }
    }
  }

  # fin del programa para el STATICO

}

# inicio del programa para realizar un COSTATIS

COSTATIS<-function(X,Y,dimX=NULL,dimY=NULL,coloresf=NULL,coloresc1=NULL,coloresc2=NULL,norm=FALSE,contr=FALSE)
{

  # lee los dos cubo de datos con las etiquetas de las filas, las columnas y las repeticiones (si no estan incluidas,
  # se nombran por defecto), comprueba que los dos cubos tengan las mismas filas y repeticiones y suprime las filas
  # que tengan datos faltantes

  l<-read(X,Y)
  X<-l[[1]]
  filas<-l[[2]]
  columnas1<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc1<-l[[6]]
  namesr<-l[[7]]
  Y<-l[[8]]
  columnas2<-l[[9]]
  namesc2<-l[[10]]
  conf<-l[[11]]

  # lee los colores de las filas y de las columnas de los dos cubos y comprueba que hay tantos como filas y como
  # columnas (si no se dan, se asignan por defecto en negro)

  if(!conf)
  {
    l<-colores(filas,columnas1,coloresf,coloresc1,conf,columnas2,coloresc2)
    coloresf<-l[[1]]
    coloresc1<-l[[2]]
    coloresc2<-l[[3]]
    conf<-l[[4]]
  }

  # centra los cubos de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por columnas

  if(!conf)
  {
    X<-preproc(X,norm,TRUE)
    Y<-preproc(Y,norm,TRUE)

    # crea la matriz con los pesos uniformes para las filas, las repeticiones y la raiz cuadrada de esta

    l<-pesos(rep(1,filas))
    Dn<-l[[1]]
    l<-pesos(rep(1,repeticiones))
    Dk<-l[[1]]
    Dk2<-l[[2]]

    # calcula las matrices de varianzas-covarianzas vectoriales

    l<-as.matrix(expand.grid(1:repeticiones,1:repeticiones))
    CovvX<-matrix(apply(l,1,function(v,X1,D){return(sum(diag(t(X1[,,v[1]])%*%D%*%X1[,,v[2]])))},X1=X,D=Dn),
                  nrow=repeticiones,ncol=repeticiones)
    CovvY<-matrix(apply(l,1,function(v,X1,D){return(sum(diag(t(X1[,,v[1]])%*%D%*%X1[,,v[2]])))},X1=Y,D=Dn),
                  nrow=repeticiones,ncol=repeticiones)

    # extrae los vectores propios segun el metodo explicado y calcula las coordenadas para las interestructuras

    aX<-(eigen(CovvX%*%Dk,symmetric=TRUE))$vectors
    VIX<-solve(Dk2)%*%(aX[,1:2])%*%diag(1/sqrt(diag(t(aX[,1:2])%*%(aX[,1:2]))))
    VIX[,1]<-abs(VIX[,1])
    IX<-CovvX%*%Dk%*%VIX
    aY<-(eigen(CovvY%*%Dk,symmetric=TRUE))$vectors
    VIY<-solve(Dk2)%*%(aY[,1:2])%*%diag(1/sqrt(diag(t(aY[,1:2])%*%(aY[,1:2]))))
    VIY[,1]<-abs(VIY[,1])
    IY<-CovvY%*%Dk%*%VIY

    # calcula las matrices compromiso

    Xc<-apply(X,1:2,function(v,v2){return(sum(v*v2))},v2=VIX[,1]/sum(VIX[,1]))
    Yc<-apply(Y,1:2,function(v,v2){return(sum(v*v2))},v2=VIY[,1]/sum(VIY[,1]))

    # si se introdujo TRUE en normalizacion, centra y normaliza los compromisos por columnas

    Xc<-preproc(Xc,norm)
    Yc<-preproc(Yc,norm)

    # extrae los vectores singulares por la izquierda y la derecha y los valores singulares segun el metodo explicado

    c<-svd(t(Yc)%*%Dn%*%Xc)
    d<-c$d
    u<-c$u
    v<-c$v
    Ur<-u%*%diag(1/sqrt(diag(t(u)%*%u)))
    Vr<-v%*%diag(1/sqrt(diag(t(v)%*%v)))

    # calcula las coordenadas para las filas y las columnas de los compromisos y de las trayectorias mediante el
    # diagrama de dualidad

    FXc<-Xc%*%Vr
    CXc<-t(Xc)%*%Dn%*%Yc%*%Ur
    FYc<-Yc%*%Ur
    CYc<-t(Yc)%*%Dn%*%Xc%*%Vr
    FXt<-array(apply(X,3,function(m,m2){return(m%*%m2)},m2=Vr),dim=c(filas,columnas2,repeticiones))
    CXt<-array(apply(X,3,function(m,m2){return(t(m)%*%m2)},m2=Dn%*%Yc%*%Ur),dim=c(columnas1,columnas2,repeticiones))
    FYt<-array(apply(Y,3,function(m,m2){return(m%*%m2)},m2=Ur),dim=c(filas,columnas2,repeticiones))
    CYt<-array(apply(Y,3,function(m,m2){return(t(m)%*%m2)},m2=Dn%*%Xc%*%Vr),dim=c(columnas2,columnas2,repeticiones))
    FX<-matrix(aperm(FXt,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)
    CX<-matrix(aperm(CXt,c(2,1,3)),nrow=columnas1*repeticiones,ncol=columnas2,byrow=TRUE)
    FY<-matrix(aperm(FYt,c(2,1,3)),nrow=filas*repeticiones,ncol=columnas2,byrow=TRUE)
    CY<-matrix(aperm(CYt,c(2,1,3)),nrow=columnas2*repeticiones,ncol=columnas2,byrow=TRUE)
    CX2<-rbind(CXc,CX)
    CY2<-rbind(CYc,CY)

    # si se ha elegido, calcula las contribuciones para las filas y las columnas de las matrices compromiso y de las
    # trayectorias de las filas y las columnas

    if(contr)
    {
      contributions(FXc,namesf,"las filas segun el primer cubo en el compromiso",TRUE)
      contributions(CXc,namesc1,"las columnas del primer cubo en el compromiso")
      contributions(FYc,namesf,"las filas segun el segundo cubo en el compromiso")
      contributions(CYc,namesc2,"las columnas del segundo cubo en el compromiso")
      contributions(FX,rep(namesf,repeticiones),"las filas segun el primer cubo en todas las repeticiones")
      contributions(CX,rep(namesc1,repeticiones),"las columnas del primer cubo en todas las repeticiones")
      contributions(FY,rep(namesf,repeticiones),"las filas segun el segundo cubo en todas las repeticiones")
      contributions(CY,rep(namesc2,repeticiones),"las columnas del segundo cubo en todas las repeticiones")
    }

    # si no se ha elegido representar los ejes, representa los valores singulares en un diagrama de sedimentacion

    if((is.null(dimX))||(is.null(dimY)))
    {
      dev.new()
      layout(matrix(c(1,1,2,1),ncol=2),widths=c(310.5/99,382.5/99),heights=c(5.5,1.5))
      layout.show(2)
      screeplot(d)

      # representa un esquema con la forma de los cubos de datos

      plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,4.25),ylim=c(-7.25,2.25),bty="n",xaxt="n",yaxt="n")
      rect(-1,-1,1,1)
      rect(1,-1,3,1)
      rect(-1,-4,1,-2)
      rect(1,-4,3,-2)
      rect(0,-7,2,-5)
      segments(-1,1,0,2)
      segments(0,2,4,2)
      segments(4,2,4,0)
      segments(4,0,3,-1)
      segments(1,1,2,2)
      segments(3,1,4,2)
      text(0,0,"X")
      text(2,0,"Y")
      text(0,-3,"Xc")
      text(2,-3,"Yc")
      text(1,-6,"Yc' Dn Xc")
      text(-1.25,0,filas)
      text(0,1.25,columnas1)
      text(2,1.25,columnas2)
      text(-0.75,1.75,repeticiones,srt=45)
      text(-1.25,-3,filas)
      text(0,-1.75,columnas1)
      text(2,-1.75,columnas2)
      text(-0.25,-6,columnas2)
      text(1,-4.75,columnas1)

      # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

    } else {
      if(compr(dimX,dimY,d))
      {
        dimension<-c(dimX,dimY)
        FXcd<-FXc[,dimension]
        CXcd<-CXc[,dimension]
        FYcd<-FYc[,dimension]
        CYcd<-CYc[,dimension]
        FXtd<-FXt[,dimension,]
        CXtd<-CXt[,dimension,]
        FYtd<-FYt[,dimension,]
        CYtd<-CYt[,dimension,]
        FXd<-FX[,dimension]
        FYd<-FY[,dimension]
        CX2d<-CX2[,dimension]
        CY2d<-CY2[,dimension]
        Fc<-rbind(FXcd,FYcd)

        # crea el entorno del grafico de las interestructuras, y se divide en dos

        dev.new(width=14,height=7,noRStudioGD=TRUE)
        layout(matrix(1:2,nrow=1))
        layout.show(2)

        # representa el grafico para las repeticiones segun el primer cubo en la izquierda con las etiquetas y vectores
        # desde el origen

        plot(IX[,1],IX[,2],type="n",xlim=c(min(0,IX[,1]),max(0,IX[,1])),ylim=c(min(0,IX[,2]),max(0,IX[,2])),
             xlab="",ylab="")
        text(IX[,1],IX[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,IX[,1],IX[,2],angle=10,length=0.1)

        # representa las repeticiones segun el segundo cubo en la derecha con las etiquetas y vectores desde el origen

        plot(IY[,1],IY[,2],type="n",xlim=c(min(0,IY[,1]),max(0,IY[,1])),ylim=c(min(0,IY[,2]),max(0,IY[,2])),
             xlab="",ylab="")
        text(IY[,1],IY[,2],namesr,cex=0.8,pos=2)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        arrows(0,0,IY[,1],IY[,2],angle=10,length=0.1)

        # representa el grafico para la coinercia de las filas de los compromisos en la izquierda con las etiquetas y
        # vectores de la segunda a la primera matriz segun los colores de los grupos
        # representa las columnas de la primera matriz compromiso en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimension,d=d,M1=FXcd,M2=CXcd,M3=FYcd,lim1=Fc,lim2=CX2d,names1=namesf,names2=namesc1,
              colores1=coloresf,colores2=coloresc1,contf=contributions2(FXc,dimension),
              contc=contributions2(CXc,dimension))

        # representa el grafico para la coinercia de las filas de los compromisos en la izquierda con las etiquetas y
        # vectores de la primera a la segunda matriz segun los colores de los grupos
        # representa las columnas de la segunda matriz compromiso en la derecha con las etiquetas y vectores desde el
        # origen segun su color

        plotm(dim=dimension,d=d,M1=FYcd,M2=CYcd,M3=FXcd,lim1=Fc,lim2=CY2d,names1=namesf,names2=namesc2,
              colores1=coloresf,colores2=coloresc2,contf=contributions2(FYc,dimension),
              contc=contributions2(CYc,dimension))

        # representa el grafico para las trayectorias de las filas del primer cubo en la izquierda con las etiquetas y
        # segmentos segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del primer cubo en la derecha con las etiquetas, segmentos y
        # vectores desde el origen segun su color

        plotm(dim=dimension,d=d,M1=FXtd[,,1],M2=CXtd[,,1],M4=FXtd,M5=CXtd,lim1=FXd,lim2=CX2d,names1=namesf,
              names2=namesc1,colores1=coloresf,colores2=coloresc1,contf=contributions2(FXt[,,1],dimension),
              contc=contributions2(CXt[,,1],dimension))

        # representa el grafico para las trayectorias de las filas del segundo cubo en la izquierda con las etiquetas y
        # segmentos segun los colores de los grupos a los que pertenecen
        # representa las trayectorias de las columnas del segundo cubo en la derecha con las etiquetas, segmentos y
        # vectores desde el origen segun su color

        plotm(dim=dimension,d=d,M1=FYtd[,,1],M2=CYtd[,,1],M4=FYtd,M5=CYtd,lim1=FYd,lim2=CY2d,names1=namesf,
              names2=namesc2,colores1=coloresf,colores2=coloresc2,contf=contributions2(FYt[,,1],dimension),
              contc=contributions2(CYt[,,1],dimension))
      }
    }
  }

  # fin del programa para el COSTATIS

}

# inicio del programa para realizar un TUCKER3

TUCKER3<-function(X,p=NULL,q=NULL,r=NULL,P1=NULL,P2=NULL,Q1=NULL,Q2=NULL,R1=NULL,R2=NULL,coloresf=NULL,
                  coloresc=NULL,maximo=5,iter=100,tol=10^-8,norm=FALSE,contr=FALSE)
{

  # lee el cubo de datos con las etiquetas de las filas, las columnas y las repeticiones (si no estan incluidas,
  # se nombran por defecto) y suprime las filas que tengan datos faltantes

  l<-read(X)
  X<-l[[1]]
  filas<-l[[2]]
  columnas<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc<-l[[6]]
  namesr<-l[[7]]
  conf<-l[[11]]

  # lee las etiquetas de los colores para las filas y columnas y comprueba que hay tantas como filas y columnas
  # (si no se dan, se asignan por defecto de color negro)

  if(!conf)
  {
    l<-colores(filas,columnas,coloresf,coloresc,conf)
    coloresf<-l[[1]]
    coloresc<-l[[2]]
    conf<-l[[3]]
  }

  # centra el cubo de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por capas laterales

  if(!conf)
  {
    X<-preproc(X,norm,TRUE,TRUE)

    # comprueba los valores de maximo y de iteraciones

    if((maximo<=1)||(maximo!=(floor(maximo))))
    {
      conf<-TRUE
      print("te has confundido, numero maximo de componentes incorrecto")
    }
    if(!conf){if((iter<1)||(iter!=(floor(iter))))
    {
      conf<-TRUE
      print("te has confundido, numero maximo de iteraciones incorrecto")
    }
    }
    if(!conf)
    {
      c1<-min(maximo,filas)
      c2<-min(maximo,columnas)
      c3<-min(maximo,repeticiones)

      # si no se ha elegido una combinacion de componentes, realiza el analisis para todas las combinaciones posibles

      if((is.null(p))||(is.null(q))||(is.null(r)))
      {

        # crea una tabla con tantas combinaciones de componentes como sean posibles

        l<-as.matrix(expand.grid(1:c3,1:c2,1:c1))
        l<-l[,3:1]
        names<-apply(l,1,function(v){return(paste(v[1],"x",v[2],"x",v[3],sep=""))})

        # para cada combinacion que cumpla la regla del maximo producto realiza la funcion auxiliar y construye la tabla
        # con los resultados

        a<-function(cont,l,T,iter,tol,c1,c2,c3)
        {
          if(((l[cont,1])<=((l[cont,2])*(l[cont,3])))&&((l[cont,2])<=((l[cont,1])*(l[cont,3])))&&
             ((l[cont,3])<=((l[cont,1])*(l[cont,2]))))
          {
            return(c(cont,every(l[cont,1],l[cont,2],l[cont,3],T,iter,tol)))
          } else {return(rep(0,5))}
        }
        table1<-matrix(unlist(lapply(1:(dim(l)[1]),a,l=l,T=X,iter=iter,tol=tol,c1=c1,c2=c2,c3=c3)),ncol=5,
                       byrow=TRUE)

        # elimina de la tabla las combinaciones que no cumplen la regla del maximo producto

        cual<-apply(table1,1,function(v){return(any(v!=0))})
        table1<-table1[cual,]
        names<-names[cual]
        filas1<-dim(table1)[1]
        ini<-rbind(c("Number","Model_Size","Sum","Best_given_Sum","SS(Res)","Prop._SS(Fit)","Number_of_iterations"),
                   cbind(table1[,1],names,table1[,2],rep("",filas1),table1[,3:5]))

        # crea otra tabla con las combinaciones con un mejor ajuste para cada valor de la suma de componentes

        cual<-unlist(lapply(3:(c1+c2+c3),function(i,m){if(any(m[,2]==i)){min(which(m[,4]==max(m[m[,2]==i,4])))}},
                            m=table1))
        table2<-table1[cual,]
        names2<-names[cual]
        ini[cual+1,4]<-"*"
        filas2<-dim(table2)[1]

        # guarda las tablas como resultados

        results(ini,"todas las combinaciones",first=TRUE)
        ini2<-rbind(c("Number","Model_Size","S","SS(Res)","DifFit","Prop._SS(Fit)","Number_of_iterations"),
                    cbind(table2[,1],names2,table2[,2:3],table2[,4]-c(0,table2[-filas2,4]),table2[,4:5]))
        results(ini2,"combinaciones con mejor ajuste")

        # crea una tabla con las combinaciones que pertenecen a la envolvente convexa de entre todas

        table3<-table2[table2[,3]<=(table2[,2]-3)*(table2[filas2,3]-table2[1,3])/(c1+c2+c3-3)+table2[1,3],]
        cual<-chull(table3[,2:3])
        table3<-table3[cual[order(cual)],]
        filas3<-dim(table3)[1]

        # elimina todas las combinaciones mas estables menos la mas simple de entre estas

        cual<-ifelse(any((table3[,3]-c(table3[-1,3],0))<table3[,3]/100),
                     min(which((table3[,3]-c(table3[-1,3],0))<table3[,3]/100)),filas3)
        table4<-table3[1:cual,]
        filas4<-dim(table4)[1]
        if(filas4>=3)
        {
          st<-c(0,((table4[-c(filas4-1,filas4),3]-table4[-c(1,filas4),3])/
                     (table4[-c(1,filas4),2]-table4[-c(filas4-1,filas4),2]))/
                  ((table4[-c(1,filas4),3]-table4[-c(1,2),3])/(table4[-c(1,2),2]-table4[-c(1,filas4),2])),0)
        } else {st<-0}
        if(any(st!=0)){vertical<-table4[min(which(st==min(st[st!=0])))-1,2]+0.5}

        # representa el grafico scree-plot con todas las combinaciones originales, con vectores para las que pertenecen a
        # la envolvente convexa

        dev.new(width=14,height=7,noRStudioGD=TRUE)
        layout(matrix(c(1,1,2,1),ncol=2),widths=c(9.9375,4.0625),heights=c(3.5,3.5))
        layout.show(2)
        par(mar=c(5,8,4,2)+0.1)
        plot(table1[,2],table1[,3],type="n",xlim=c(min(table1[,2])-1,max(table1[,2])+1),
             ylim=c(min(table1[,3])-1,max(table1[,3])+1),xlab="Suma del numero de componentes",
             ylab="Suma de cuadrados residual")
        text(table1[,2],table1[,3],names,cex=0.8,pos=4)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        points(table1[,2],table1[,3],pch=3,col=rgb(255,0,0,maxColorValue=255),cex=0.8,lwd=2)
        lapply(1:(filas3-1),function(i,M)
        {
          arrows(M[i,2],M[i,3],M[i+1,2],M[i+1,3],col=rgb(0,0,255,maxColorValue = 255),angle=10,length=0.1,lwd=1.5)
        },M=table3)

        # representa una recta vertical separando las combinaciones estables de las no estables, asi la primera a su
        # derecha podria ser la elegida para el resto del analisis

        if(any(st!=0)){abline(v=vertical,col=rgb(255,0,127,maxColorValue=255),lwd=2)}

        # representa un esquema con la forma del cubo de datos

        plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,2.25),ylim=c(-1.25,2.25),bty="n",xaxt="n",yaxt="n")
        rect(-1,-1,1,1)
        segments(-1,1,0,2)
        segments(0,2,2,2)
        segments(2,2,2,0)
        segments(2,0,1,-1)
        segments(1,1,2,2)
        text(0,0,"X")
        text(-1.25,0,filas)
        text(0,1.25,columnas)
        text(-0.75,1.75,repeticiones,srt=45)

        # si se ha elegido una combinacion de componentes, comprueba que sean correctos

      } else {
        if((p<1)||(p!=(floor(p)))||(p>c1))
        {
          conf<-TRUE
          print("te has confundido, numero de componentes incorrecto")
        }
        if(!conf)
        {
          if((q<1)||(q!=(floor(q)))||(q>c2))
          {
            conf<-TRUE
            print("te has confundido, numero de componentes incorrecto")
          }
        }
        if(!conf)
        {
          if((r<1)||(r!=(floor(r)))||(r>c3))
          {
            conf<-TRUE
            print("te has confundido, numero de componentes incorrecto")
          }
        }
        if(!conf)
        {
          comp<-c(p,q,r)

          # realiza la funcion auxiliar para la combinacion de componentes establecida

          l<-every(p,q,r,X,iter,tol,mas=TRUE)
          A<-l[[1]]
          B<-l[[2]]
          C<-l[[3]]
          G<-l[[4]]
          Aprox<-tensorial(tensorial(tensorial(G,1,t(A)),2,t(B)),3,t(C))

          # crea una tabla con los porcentajes de ajuste para cada componente de cada dimension, y el ajuste total para
          # cada dimension, y la guarda como resultado

          sos<-100*(matrix(unlist(lapply(1:3,function(i,T,comp)
          {
            return(c(apply(T^2,i,sum),rep(0,max(comp)-comp[i])))
          },T=G,comp=comp)),ncol=3))/(sum(X^2))
          sos<-rbind(sos,apply(sos,2,sum))
          ini3<-rbind(c("Component","Dimension_1","Dimension_2","Dimension_3"),
                      cbind(c(1:(max(comp)),"Total_explained_variation"),sos))
          results(ini3,"porcentajes de ajuste",first=TRUE)

          # crea otra tabla con el tensor core y el porcentaje de ajuste para cada combinacion de componentes de las tres
          # dimensiones, y la guarda como resultado

          core<-cbind(matrix(aperm(G,c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE),
                      matrix(aperm(100*(G^2)/(sum(X^2)),c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE))
          ini4<-rbind(c(rep("",3),rep(c("Mode_2_components",rep("",q-1)),2)),
                      c(rep("",3),"Residual_Sums_of_Squares",rep("",q-1),"Explained_Variation",rep("",q-1)),
                      c(rep("",3),rep(1:q,2)),cbind(unlist(lapply(1:r,function(i,p)
                      {
                        return(c(paste("Mode_3,_Component_",i,sep=""),rep("",p-1)))
                      },p=p)),rep(c("Mode_1_components",rep("",p-1)),r),rep(1:p,r),core))
          results(ini4,"core array")

          # calcula el ajuste de cada fila, columna y repeticion

          fit1<-apply(Aprox^2,1,sum)
          fit2<-apply(Aprox^2,2,sum)
          fit3<-apply(Aprox^2,3,sum)
          fit1<-cbind(fit1,apply(X^2,1,sum)-fit1)
          fit2<-cbind(fit2,apply(X^2,2,sum)-fit2)
          fit3<-cbind(fit3,apply(X^2,3,sum)-fit3)

          # guarda como resultado las coordenadas de cada fila, columna y repeticion para los graficos

          results(A,"coordenadas de las filas",names=namesf,axis=TRUE)
          results(B,"coordenadas de las columnas",names=namesc,axis=TRUE)
          results(C,"coordenadas de las repeticiones",names=namesr,axis=TRUE)

          # si se ha elegido, calcula las contribuciones para las filas, las columnas y las repeticiones del cubo de datos

          if(contr)
          {
            contributions(A,namesf,"las filas")
            contributions(B,namesc,"las columnas")
            contributions(C,namesr,"las repeticiones")
          }

          # representa el grafico con el ajuste para cada fila, columna y repeticion con las etiquetas segun los colores de
          # los grupos a los que pertenecen
          # ademas pinta una recta separando las filas, columnas y repeticiones con un ajuste superior o inferior a la media

          dev.new(width=14,height=7,noRStudioGD=TRUE)
          layout(matrix(1:3,nrow=1))
          layout.show(3)
          par(mar=c(5,6,4,2)+0.1)
          plot(fit1[,1],fit1[,2],type="n",xlim=c(min(fit1[,1]),max(fit1[,1])),ylim=c(min(fit1[,2]),max(fit1[,2]))
               ,xlab="",ylab="")
          text(fit1[,1],fit1[,2],namesf,col=coloresf)
          box(lwd=2)
          abline(h=0,lwd=2)
          abline(v=0,lwd=2)
          abline(0,sum(fit1[,2])/sum(fit1[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
          plot(fit2[,1],fit2[,2],type="n",xlim=c(min(fit2[,1]),max(fit2[,1])),ylim=c(min(fit2[,2]),max(fit2[,2])),
               xlab="",ylab="")
          text(fit2[,1],fit2[,2],namesc,col=coloresc)
          box(lwd=2)
          abline(h=0,lwd=2)
          abline(v=0,lwd=2)
          abline(0,sum(fit2[,2])/sum(fit2[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
          plot(fit3[,1],fit3[,2],type="n",xlim=c(min(fit3[,1]),max(fit3[,1])),ylim=c(min(fit3[,2]),max(fit3[,2])),
               xlab="",ylab="")
          text(fit3[,1],fit3[,2],namesr)
          box(lwd=2)
          abline(h=0,lwd=2)
          abline(v=0,lwd=2)
          abline(0,sum(fit3[,2])/sum(fit3[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)

          # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

          if((!is.null(P1))&&(!is.null(P2))&&(!is.null(Q1))&&(!is.null(Q2))&&(!is.null(R1))&&(!is.null(R2)))
          {
            if((compr(P1,P2,c1,p,TRUE))&&(compr(Q1,Q2,c2,q,TRUE))&&(compr(R1,R2,c3,r,TRUE)))
            {
              dimensionP<-c(P1,P2)
              dimensionQ<-c(Q1,Q2)
              dimensionR<-c(R1,R2)
              Ad<-A[,dimensionP]
              Bd<-B[,dimensionQ]
              Cd<-C[,dimensionR]

              # representa las filas en la izquierda con las etiquetas segun los colores de los grupos a los que pertenecen
              # representa las columnas en el centro con las etiquetas y vectores desde el origen segun los colores de los
              # grupos a los que pertenecen
              # representa las repeticiones en la derecha con las etiquetas y vectores desde el origen

              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sos,M1=Ad,M2=Bd,M3=Cd,lim1=Ad,lim2=Bd,
                     lim3=Cd,names1=namesf,names2=namesc,names3=namesr,colores1=coloresf,colores2=coloresc)

              contf<-contributions2(A,dimensionP)
              contc<-contributions2(B,dimensionQ)
              contre<-contributions2(C,dimensionR)
              contf[is.nan(contf)]<-0
              contc[is.nan(contc)]<-0
              contre[is.nan(contre)]<-0
              mcontf<-max(contf)/2
              mcontc<-max(contc)/2
              mcontre<-max(contre)/2
              lcoloresf<-coloresf[contf>=mcontf]
              lcoloresc<-coloresc[contc>=mcontc]
              lnamesf<-namesf[contf>=mcontf]
              lnamesc<-namesc[contc>=mcontc]
              lnamesr<-namesr[contre>=mcontre]
              lA<-matrix(Ad[contf>=mcontf,],ncol=2)
              lB<-matrix(Bd[contc>=mcontc,],ncol=2)
              lC<-matrix(Cd[contre>=mcontre,],ncol=2)
              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sos,M1=lA,M2=lB,M3=lC,lim1=Ad,lim2=Bd,
                     lim3=Cd,names1=lnamesf,names2=lnamesc,names3=lnamesr,colores1=lcoloresf,colores2=lcoloresc,
                     titles=TRUE)
            }
          }
        }
      }
    }
  }

  # fin del programa para el TUCKER3

}

# inicio del programa para realizar un COTUCKER3

COTUCKER3<-function(X,Y,p=NULL,q=NULL,r=NULL,P1=NULL,P2=NULL,Q1=NULL,Q2=NULL,R1=NULL,R2=NULL,
                    dimAX=NULL,dimAY=NULL,dimBX=NULL,dimBY=NULL,dimCX=NULL,dimCY=NULL,
                    coloresf=NULL,coloresc1=NULL,coloresc2=NULL,maximo=5,iter=100,tol=10^-8,norm=FALSE,contr=FALSE)
{

  # lee los cubos de datos con las etiquetas de las filas, las columnas y las repeticiones
  # (si no estan incluidas, se nombran por defecto) y suprime las filas que tengan datos faltantes

  l<-read(X,Y)
  X<-l[[1]]
  filas<-l[[2]]
  columnas1<-l[[3]]
  repeticiones<-l[[4]]
  namesf<-l[[5]]
  namesc1<-l[[6]]
  namesr<-l[[7]]
  Y<-l[[8]]
  columnas2<-l[[9]]
  namesc2<-l[[10]]
  conf<-l[[11]]

  # lee las etiquetas de los colores para las filas y columnas y comprueba que hay tantas como filas y columnas
  # (si no se dan, se asignan por defecto de color negro)

  if(!conf)
  {
    l<-colores(filas,columnas1,coloresf,coloresc1,conf,columnas2,coloresc2)
    coloresf<-l[[1]]
    coloresc1<-l[[2]]
    coloresc2<-l[[3]]
    conf<-l[[4]]
  }

  # centra los cubos de datos por columnas y si se introdujo TRUE en normalizacion, normaliza por capas laterales

  if(!conf)
  {
    X<-preproc(X,norm,TRUE,TRUE)
    Y<-preproc(Y,norm,TRUE,TRUE)

    # comprueba los valores de maximo y de iteraciones

    if((maximo<=1)||(maximo!=(floor(maximo))))
    {
      conf<-TRUE
      print("te has confundido, numero maximo de componentes incorrecto")
    }
    if(!conf)
    {
      if((iter<1)||(iter!=(floor(iter))))
      {
        conf<-TRUE
        print("te has confundido, numero maximo de iteraciones incorrecto")
      }
    }
    if(!conf)
    {
      c1<-min(maximo,filas)
      c2<-min(maximo,columnas1,columnas2)
      c3<-min(maximo,repeticiones)

      # si no se ha elegido una combinacion de componentes, realiza el analisis para todas las combinaciones posibles

      if((is.null(p))||(is.null(q))||(is.null(r)))
      {

        # crea una tabla con tantas combinaciones de componentes como sean posibles

        l<-as.matrix(expand.grid(1:c3,1:c2,1:c1))
        l<-l[,3:1]
        names<-apply(l,1,function(v){return(paste(v[1],"x",v[2],"x",v[3],sep=""))})

        # para cada combinacion que cumpla la regla del maximo producto realiza la funcion auxiliar y construye la tabla
        # con los resultados

        a<-function(cont,l,T,iter,tol,c1,c2,c3)
        {
          if(((l[cont,1])<=((l[cont,2])*(l[cont,3])))&&((l[cont,2])<=((l[cont,1])*(l[cont,3])))&&
             ((l[cont,3])<=((l[cont,1])*(l[cont,2]))))
          {
            return(c(cont,every(l[cont,1],l[cont,2],l[cont,3],T,iter,tol)))
          } else {return(rep(0,5))}
        }
        table1X<-matrix(unlist(lapply(1:(dim(l)[1]),a,l=l,T=X,iter=iter,tol=tol,c1=c1,c2=c2,c3=c3)),ncol=5,byrow=TRUE)
        table1Y<-matrix(unlist(lapply(1:(dim(l)[1]),a,l=l,T=Y,iter=iter,tol=tol,c1=c1,c2=c2,c3=c3)),ncol=5,byrow=TRUE)

        # elimina de la tabla las combinaciones que no cumplen la regla del maximo producto

        cual<-apply(table1X,1,function(v){return(any(v!=0))})
        table1X<-table1X[cual,]
        table1Y<-table1Y[cual,]
        names<-names[cual]
        filas1<-dim(table1X)[1]
        table1<-cbind(table1X[,1],table1X[,2],table1X[,3]+table1Y[,3],
                      100-100*(table1X[,3]+table1Y[,3])/(100*(table1X[,3]/(100-table1X[,4])+table1Y[,3]/(100-table1Y[,4]))),
                      apply(cbind(table1X[,5],table1Y[,5]),1,max))
        ini<-rbind(c("Number","Model_Size","Sum","Best_given_Sum","SS(Res)","Prop._SS(Fit)","Number_of_iterations"),
                   cbind(table1[,1],names,table1[,2],rep("",filas1),table1[,3:5]))

        # crea otra tabla con las combinaciones con un mejor ajuste para cada valor de la suma de componentes

        cual<-unlist(lapply(3:(c1+c2+c3),function(i,m){if(any(m[,2]==i)){min(which(m[,4]==max(m[m[,2]==i,4])))}},m=table1))
        table2<-table1[cual,]
        names2<-names[cual]
        ini[cual+1,4]<-"*"
        filas2<-dim(table2)[1]

        # guarda las tablas como resultados

        results(ini,"todas las combinaciones",first=TRUE)
        ini2<-rbind(c("Number","Model_Size","S","SS(Res)","DifFit","Prop._SS(Fit)","Number_of_iterations"),
                    cbind(table2[,1],names2,table2[,2:3],table2[,4]-c(0,table2[-filas2,4]),table2[,4:5]))
        results(ini2,"combinaciones con mejor ajuste")

        # crea una tabla con las combinaciones que pertenecen a la envolvente convexa de entre todas

        cual<-table2[,3]<=((table2[,2]-3)*(table2[filas2,3]-table2[1,3])/(c1+c2+c3-3)+table2[1,3])
        if((table2[filas2,2])==(c1+c2+c3)){cual[filas2]<-TRUE}
        table3<-table2[cual,]
        cual<-chull(table3[,2:3])
        table3<-table3[cual[order(cual)],]
        filas3<-dim(table3)[1]

        # elimina todas las combinaciones mas estables menos la mas simple de entre estas

        cual<-ifelse(any((table3[,3]-c(table3[-1,3],0))<table3[,3]/100),
                     min(which((table3[,3]-c(table3[-1,3],0))<table3[,3]/100)),filas3)
        table4<-table3[1:cual,]
        filas4<-dim(table4)[1]
        if(filas4>=3)
        {
          st<-c(0,((table4[-c(filas4-1,filas4),3]-table4[-c(1,filas4),3])/
                     (table4[-c(1,filas4),2]-table4[-c(filas4-1,filas4),2]))/
                  ((table4[-c(1,filas4),3]-table4[-c(1,2),3])/(table4[-c(1,2),2]-table4[-c(1,filas4),2])),0)
        } else {st<-0}
        if(any(st!=0)){vertical<-table4[min(which(st==min(st[st!=0])))-1,2]+0.5}

        # representa el grafico scree-plot con todas las combinaciones originales, con vectores para las que pertenecen
        # a la envolvente convexa

        dev.new(width=14,height=7,noRStudioGD=TRUE)
        layout(matrix(c(1,1,2,1),ncol=2),widths=c(8.85,5.15),heights=c(3.5,3.5))
        layout.show(2)
        par(mar=c(5,8,4,2)+0.1)
        plot(table1[,2],table1[,3],type="n",xlim=c(min(table1[,2])-1,max(table1[,2])+1),
             ylim=c(min(table1[,3])-1,max(table1[,3])+1),xlab="Suma del numero de componentes",
             ylab="Suma de cuadrados residual")
        text(table1[,2],table1[,3],names,cex=0.8,pos=4)
        box(lwd=2)
        abline(h=0,lwd=2)
        abline(v=0,lwd=2)
        points(table1[,2],table1[,3],pch=3,col=rgb(255,0,0,maxColorValue=255),cex=0.8,lwd=2)
        lapply(1:(filas3-1),function(i,M){arrows(M[i,2],M[i,3],M[i+1,2],M[i+1,3],col=rgb(0,0,255,maxColorValue = 255),
                                                 angle=10,length=0.1,lwd=1.5)},M=table3)

        # representa una recta vertical separando las combinaciones estables de las no estables, asi la primera a su derecha
        # podria ser la elegida para el resto del analisis

        if(any(st!=0)){abline(v=vertical,col=rgb(255,0,127,maxColorValue=255),lwd=2)}

        # representa un esquema con la forma de los cubos de datos

        plot(0,0,type="n",xlab="",ylab="",xlim=c(-1.5,4.25),ylim=c(-1.25,2.25),bty="n",xaxt="n",yaxt="n")
        rect(-1,-1,1,1)
        rect(1,-1,3,1)
        segments(-1,1,0,2)
        segments(0,2,4,2)
        segments(4,2,4,0)
        segments(4,0,3,-1)
        segments(1,1,2,2)
        segments(3,1,4,2)
        text(0,0,"X")
        text(2,0,"Y")
        text(-1.25,0,filas)
        text(0,1.25,columnas1)
        text(2,1.25,columnas2)
        text(-0.75,1.75,repeticiones,srt=45)

        # si se ha elegido una combinacion de componentes, comprueba que sean correctos

      } else {
        if((p<1)||(p!=(floor(p)))||(p>c1))
        {
          conf<-TRUE
          print("te has confundido, numero de componentes incorrecto")
        }
        if(!conf)
        {
          if((q<1)||(q!=(floor(q)))||(q>c2))
          {
            conf<-TRUE
            print("te has confundido, numero de componentes incorrecto")
          }
        }
        if(!conf)
        {
          if((r<1)||(r!=(floor(r)))||(r>c3))
          {
            conf<-TRUE
            print("te has confundido, numero de componentes incorrecto")
          }
        }
        if(!conf)
        {
          comp<-c(p,q,r)

          # realiza la funcion auxiliar para la combinacion de componentes establecida

          l<-every(p,q,r,X,iter,tol,mas=TRUE)
          AX<-l[[1]]
          BX<-l[[2]]
          CX<-l[[3]]
          GX<-l[[4]]
          AproxX<-tensorial(tensorial(tensorial(GX,1,t(AX)),2,t(BX)),3,t(CX))
          l<-every(p,q,r,Y,iter,tol,mas=TRUE)
          AY<-l[[1]]
          BY<-l[[2]]
          CY<-l[[3]]
          GY<-l[[4]]
          AproxY<-tensorial(tensorial(tensorial(GY,1,t(AY)),2,t(BY)),3,t(CY))

          # crea una tabla con los porcentajes de ajuste para cada componente de cada dimension,
          # y el ajuste total para cada dimension, una tabla para cada cubo

          sosX<-100*(matrix(unlist(lapply(1:3,function(i,T,comp)
          {
            return(c(apply(T^2,i,sum),rep(0,max(comp)-comp[i])))
          },T=GX,comp=comp)),ncol=3))/(sum(X^2))
          sosY<-100*(matrix(unlist(lapply(1:3,function(i,T,comp)
          {
            return(c(apply(T^2,i,sum),rep(0,max(comp)-comp[i])))
          },T=GY,comp=comp)),ncol=3))/(sum(Y^2))
          sosX<-rbind(sosX,apply(sosX,2,sum))
          sosY<-rbind(sosY,apply(sosY,2,sum))
          ini3X<-rbind(c("Component","Dimension_1","Dimension_2","Dimension_3"),
                       cbind(c(1:(max(comp)),"Total_explained_variation"),sosX))
          ini3Y<-rbind(c("Component","Dimension_1","Dimension_2","Dimension_3"),
                       cbind(c(1:(max(comp)),"Total_explained_variation"),sosY))

          # crea otra tabla con el tensor core y el porcentaje de ajuste para combinacion de componentes de las tres dimensiones,
          # una tabla para cada cubo

          coreX<-cbind(matrix(aperm(GX,c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE),
                       matrix(aperm(100*(GX^2)/(sum(X^2)),c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE))
          coreY<-cbind(matrix(aperm(GY,c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE),
                       matrix(aperm(100*(GY^2)/(sum(Y^2)),c(2,1,3)),nrow=p*r,ncol=q,byrow=TRUE))
          ini4X<-rbind(c(rep("",3),rep(c("Mode_2_components",rep("",q-1)),2)),
                       c(rep("",3),"Residual_Sums_of_Squares",rep("",q-1),"Explained_Variation",rep("",q-1)),
                       c(rep("",3),rep(1:q,2)),cbind(unlist(lapply(1:r,function(i,p)
                       {
                         return(c(paste("Mode_3,_Component_",i,sep=""),rep("",p-1)))
                       },p=p)),rep(c("Mode_1_components",rep("",p-1)),r),rep(1:p,r),coreX))
          ini4Y<-rbind(c(rep("",3),rep(c("Mode_2_components",rep("",q-1)),2)),
                       c(rep("",3),"Residual_Sums_of_Squares",rep("",q-1),"Explained_Variation",rep("",q-1)),
                       c(rep("",3),rep(1:q,2)),cbind(unlist(lapply(1:r,function(i,p)
                       {
                         return(c(paste("Mode_3,_Component_",i,sep=""),rep("",p-1)))
                       },p=p)),rep(c("Mode_1_components",rep("",p-1)),r),rep(1:p,r),coreY))

          # calcula el ajuste de cada fila, columna de los dos cubos y repeticion

          fit1X<-apply(AproxX^2,1,sum)
          fit2X<-apply(AproxX^2,2,sum)
          fit3X<-apply(AproxX^2,3,sum)
          fit1X<-cbind(fit1X,apply(X^2,1,sum)-fit1X)
          fit2X<-cbind(fit2X,apply(X^2,2,sum)-fit2X)
          fit3X<-cbind(fit3X,apply(X^2,3,sum)-fit3X)
          fit1Y<-apply(AproxY^2,1,sum)
          fit2Y<-apply(AproxY^2,2,sum)
          fit3Y<-apply(AproxY^2,3,sum)
          fit1Y<-cbind(fit1Y,apply(Y^2,1,sum)-fit1Y)
          fit2Y<-cbind(fit2Y,apply(Y^2,2,sum)-fit2Y)
          fit3Y<-cbind(fit3Y,apply(Y^2,3,sum)-fit3Y)

          # si no se ha elegido representar los ejes, guarda como resultados las tablas con los porcentajes de ajuste
          # para cada componente de cada dimension, y el ajuste total para cada dimension

          if(!((!is.null(P1))&&(!is.null(P2))&&(!is.null(Q1))&&(!is.null(Q2))&&(!is.null(R1))&&(!is.null(R2))))
          {
            results(ini3X,"porcentajes de ajuste primer cubo",first=TRUE)
            results(ini3Y,"porcentajes de ajuste segundo cubo")

            # guarda como resultados los tensores core y los porcentajes de ajuste para cada combinacion de componentes
            # de las tres dimensiones

            results(ini4X,"core array primer cubo")
            results(ini4Y,"core array segundo cubo")

            # guarda como resultado las coordenadas de cada fila, columna de los dos cubos y repeticion para los graficos

            results(AX,"coordenadas de las filas primer cubo",names=namesf,axis=TRUE)
            results(BX,"coordenadas de las columnas primer cubo",names=namesc1,axis=TRUE)
            results(CX,"coordenadas de las repeticiones primer cubo",names=namesr,axis=TRUE)
            results(AY,"coordenadas de las filas segundo cubo",names=namesf,axis=TRUE)
            results(BY,"coordenadas de las columnas segundo cubo",names=namesc2,axis=TRUE)
            results(CY,"coordenadas de las repeticiones segundo cubo",names=namesr,axis=TRUE)

            # si se ha elegido, calcula las contribuciones para las filas, las columnas y las repeticiones de los dos cubo de datos

            if(contr)
            {
              contributions(AX,namesf,"las filas primer cubo")
              contributions(BX,namesc1,"las columnas primer cubo")
              contributions(CX,namesr,"las repeticiones primer cubo")
              contributions(AY,namesf,"las filas segundo cubo")
              contributions(BY,namesc2,"las columnas segundo cubo")
              contributions(CY,namesr,"las repeticiones segundo cubo")
            }

            # representa el grafico con el ajuste para cada fila, columna de los dos cubos y repeticion
            # con las etiquetas segun los colores de los grupos a los que pertenecen
            # ademas pinta una recta separando las filas, columnas de los dos cubos y repeticiones con un ajuste superior
            # o inferior a la media

            if((is.null(dimAX)||is.null(dimAY))&&(is.null(dimBX)||is.null(dimBY))&&(is.null(dimCX)||is.null(dimCY)))
            {
              contr=FALSE
              dev.new(width=14,height=7,noRStudioGD=TRUE)
              layout(matrix(1:3,nrow=1))
              layout.show(3)
              par(mar=c(5,6,4,2)+0.1)
              plot(fit1X[,1],fit1X[,2],type="n",xlim=c(min(fit1X[,1]),max(fit1X[,1])),
                   ylim=c(min(fit1X[,2]),max(fit1X[,2])),xlab="",ylab="")
              text(fit1X[,1],fit1X[,2],namesf,col=coloresf)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit1X[,2])/sum(fit1X[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
              plot(fit2X[,1],fit2X[,2],type="n",xlim=c(min(fit2X[,1]),max(fit2X[,1])),
                   ylim=c(min(fit2X[,2]),max(fit2X[,2])),xlab="",ylab="")
              text(fit2X[,1],fit2X[,2],namesc1,col=coloresc1)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit2X[,2])/sum(fit2X[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
              plot(fit3X[,1],fit3X[,2],type="n",xlim=c(min(fit3X[,1]),max(fit3X[,1])),
                   ylim=c(min(fit3X[,2]),max(fit3X[,2])),xlab="",ylab="")
              text(fit3X[,1],fit3X[,2],namesr)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit3X[,2])/sum(fit3X[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
              dev.new(width=14,height=7,noRStudioGD=TRUE)
              layout(matrix(1:3,nrow=1))
              layout.show(3)
              par(mar=c(5,6,4,2)+0.1)
              plot(fit1Y[,1],fit1Y[,2],type="n",xlim=c(min(fit1Y[,1]),max(fit1Y[,1])),
                   ylim=c(min(fit1Y[,2]),max(fit1Y[,2])),xlab="",ylab="")
              text(fit1Y[,1],fit1Y[,2],namesf,col=coloresf)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit1Y[,2])/sum(fit1Y[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
              plot(fit2Y[,1],fit2Y[,2],type="n",xlim=c(min(fit2Y[,1]),max(fit2Y[,1])),
                   ylim=c(min(fit2Y[,2]),max(fit2Y[,2])),xlab="",ylab="")
              text(fit2Y[,1],fit2Y[,2],namesc2,col=coloresc2)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit2Y[,2])/sum(fit2Y[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
              plot(fit3Y[,1],fit3Y[,2],type="n",xlim=c(min(fit3Y[,1]),max(fit3Y[,1])),
                   ylim=c(min(fit3Y[,2]),max(fit3Y[,2])),xlab="",ylab="")
              text(fit3Y[,1],fit3Y[,2],namesr)
              box(lwd=2)
              abline(h=0,lwd=2)
              abline(v=0,lwd=2)
              abline(0,sum(fit3Y[,2])/sum(fit3Y[,1]),col=rgb(255,0,0,maxColorValue = 255),lwd=2)
            }

            # realiza los tres analisis de co-inercia para las tres dimensiones

            dimnames(AX)[[1]]<-namesf
            dimnames(AY)[[1]]<-namesf
            dimnames(BX)[[1]]<-namesc1
            dimnames(BY)[[1]]<-namesc2
            dimnames(CX)[[1]]<-namesr
            dimnames(CY)[[1]]<-namesr
            if((is.null(dimBX)||is.null(dimBY))&&(is.null(dimCX)||is.null(dimCY)))
            {
              COIA(AX,AY,dimX=dimAX,dimY=dimAY,coloresf=coloresf,norm=norm,contr=contr,cotucker=TRUE)
            }
            if((is.null(dimAX)||is.null(dimAY))&&(is.null(dimCX)||is.null(dimCY)))
            {
              COIA(t(BX),t(BY),dimX=dimBX,dimY=dimBY,coloresc1=coloresc1,coloresc2=coloresc2,norm=norm,contr=contr,
                   cotucker=TRUE,tcotucker=TRUE)
            }
            if((is.null(dimAX)||is.null(dimAY))&&(is.null(dimBX)||is.null(dimBY)))
            {
              COIA(CX,CY,dimX=dimCX,dimY=dimCY,norm=norm,contr=contr,cotucker=TRUE)
            }

            # si se ha elegido representar los ejes, comprueba que los ejes a representar en los graficos sean correctos

          } else {
            if((compr(P1,P2,c1,p,TRUE))&&(compr(Q1,Q2,c2,q,TRUE))&&(compr(R1,R2,c3,r,TRUE)))
            {
              dimensionP<-c(P1,P2)
              dimensionQ<-c(Q1,Q2)
              dimensionR<-c(R1,R2)
              AXd<-AX[,dimensionP]
              BXd<-BX[,dimensionQ]
              CXd<-CX[,dimensionR]
              AYd<-AY[,dimensionP]
              BYd<-BY[,dimensionQ]
              CYd<-CY[,dimensionR]

              # representa las filas segun el primer cubo en la izquierda con las etiquetas segun los colores de los grupos
              # a los que pertenecen
              # representa las columnas del primer cubo en el centro con las etiquetas y vectores desde el origen
              # segun los colores de los grupos a los que pertenecen
              # representa las repeticiones segun el primer cubo en la derecha con las etiquetas y vectores desde el origen

              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sosX,M1=AXd,M2=BXd,M3=CXd,
                     lim1=AXd,lim2=BXd,lim3=CXd,names1=namesf,names2=namesc1,names3=namesr,
                     colores1=coloresf,colores2=coloresc1)
              contfX<-contributions2(AX,dimensionP)
              contcX<-contributions2(BX,dimensionQ)
              contreX<-contributions2(CX,dimensionR)
              contfX[is.nan(contfX)]<-0
              contcX[is.nan(contcX)]<-0
              contreX[is.nan(contreX)]<-0
              mcontfX<-max(contfX)/2
              mcontcX<-max(contcX)/2
              mcontreX<-max(contreX)/2
              lcoloresfX<-coloresf[contfX>=mcontfX]
              lcolorescX<-coloresc1[contcX>=mcontcX]
              lnamesfX<-namesf[contfX>=mcontfX]
              lnamescX<-namesc1[contcX>=mcontcX]
              lnamesrX<-namesr[contreX>=mcontreX]
              lAX<-matrix(AXd[contfX>=mcontfX,],ncol=2)
              lBX<-matrix(BXd[contcX>=mcontcX,],ncol=2)
              lCX<-matrix(CXd[contreX>=mcontreX,],ncol=2)
              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sosX,M1=lAX,M2=lBX,M3=lCX,
                     lim1=AXd,lim2=BXd,lim3=CXd,names1=lnamesfX,names2=lnamescX,names3=lnamesrX,
                     colores1=lcoloresfX,colores2=lcolorescX,titles=TRUE)

              # representa las filas segun el segundo cubo en la izquierda con las etiquetas segun los colores de los grupos
              # a los que pertenecen
              # representa las columnas del segundo cubo en el centro con las etiquetas y vectores desde el origen
              # segun los colores de los grupos a los que pertenecen
              # representa las repeticiones segun el segundo cubo en la derecha con las etiquetas y vectores desde el origen

              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sosY,M1=AYd,M2=BYd,M3=CYd,
                     lim1=AYd,lim2=BYd,lim3=CYd,names1=namesf,names2=namesc2,names3=namesr,
                     colores1=coloresf,colores2=coloresc2)
              contfY<-contributions2(AY,dimensionP)
              contcY<-contributions2(BY,dimensionQ)
              contreY<-contributions2(CY,dimensionR)
              contfY[is.nan(contfY)]<-0
              contcY[is.nan(contcY)]<-0
              contreY[is.nan(contreY)]<-0
              mcontfY<-max(contfY)/2
              mcontcY<-max(contcY)/2
              mcontreY<-max(contreY)/2
              lcoloresfY<-coloresf[contfY>=mcontfY]
              lcolorescY<-coloresc2[contcY>=mcontcY]
              lnamesfY<-namesf[contfY>=mcontfY]
              lnamescY<-namesc2[contcY>=mcontcY]
              lnamesrY<-namesr[contreY>=mcontreY]
              lAY<-matrix(AYd[contfY>=mcontfY,],ncol=2)
              lBY<-matrix(BYd[contcY>=mcontcY,],ncol=2)
              lCY<-matrix(CYd[contreY>=mcontreY,],ncol=2)
              plotmt(dim1=dimensionP,dim2=dimensionQ,dim3=dimensionR,sos=sosY,M1=lAY,M2=lBY,M3=lCY,
                     lim1=AYd,lim2=BYd,lim3=CYd,names1=lnamesfY,names2=lnamescY,names3=lnamesrY,
                     colores1=lcoloresfY,colores2=lcolorescY,titles=TRUE)
            }
          }
        }
      }
    }
  }

  # fin del programa para el COTUCKER3

}
