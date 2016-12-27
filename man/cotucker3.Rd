\name{COTUCKER3}
\alias{
COTUCKER3
PCA
BGA
COIA
PTA
BGCOIA
STATICO
COSTATIS
TUCKER3

read
colores
preproc
pesos
contributions
contributions2
screeplot
compr
plotm
tensorial
desplegar
results
every
plotmt
}
\title{
Co-Tucker3 analysis of two sequences of matrices
}
\description{
Computes the Co-Tucker3 (Co-Inertia Analysis + Tucker3 method) analysis of two three-dimensional arrays with the same first and third dimensions. The package also contains other functions called PCA() (Principal Component Analysis) and BGA() (Between-Groups Analysis), which perform analysis of one matrix, COIA() (Co-Inertia Analysis), which performs analysis of two matrices, PTA() (Partial Triadic Analysis) and TUCKER3(), which perform analysis of a sequence of matrices, and BGCOIA() (Between-Groups Co-Inertia Analysis), STATICO() (STATIS method + Co-Inertia Analysis), COSTATIS() (Co-Inertia Analysis + STATIS method), which also perform analysis of two sequences of matrices.
}
\usage{
COTUCKER3(X, Y, p = NULL, q = NULL, r = NULL, P1 = NULL, P2 = NULL, Q1 = NULL, Q2 = NULL,
          R1 = NULL, R2 = NULL, dimAX = NULL, dimAY = NULL, dimBX = NULL, dimBY = NULL,
          dimCX = NULL, dimCY = NULL, coloresf = NULL, coloresc1 = NULL,
          coloresc2 = NULL, maximo = 5, iter = 100, tol = 10^-8, norm = FALSE,
          contr = FALSE)
PCA(X, dimX = NULL, dimY = NULL, coloresf = NULL, coloresc = NULL, norm = FALSE,
    contr = FALSE)
BGA(X, gruposf, dimX = NULL, dimY = NULL, coloresf = NULL, coloresc = NULL, norm = FALSE,
    contr = FALSE)
COIA(X, Y, dimXx = NULL, dimYx = NULL, dimXy = NULL, dimYy = NULL, dimX = NULL,
     dimY = NULL, coloresf = NULL, coloresc1 = NULL, coloresc2 = NULL, norm = FALSE,
     contr = FALSE, cotucker = FALSE, tcotucker = FALSE)
PTA(X, dimX = NULL, dimY = NULL, coloresf = NULL, coloresc = NULL, norm = FALSE,
    contr = FALSE)
BGCOIA(X, Y, dimXx = NULL, dimYx = NULL, dimXy = NULL, dimYy = NULL, dimX = NULL,
       dimY = NULL, coloresf = NULL, coloresc1 = NULL, coloresc2 = NULL, norm = FALSE,
       contr = FALSE)
STATICO(X, Y, dimX = NULL, dimY = NULL, coloresf = NULL, coloresc1 = NULL,
        coloresc2 = NULL, norm = FALSE, contr = FALSE)
COSTATIS(X, Y, dimX = NULL, dimY = NULL, coloresf = NULL, coloresc1 = NULL,
         coloresc2 = NULL, norm = FALSE, contr = FALSE)
TUCKER3(X, p = NULL, q = NULL, r = NULL, P1 = NULL, P2 = NULL, Q1 = NULL, Q2 = NULL,
        R1 = NULL, R2 = NULL, coloresf = NULL, coloresc = NULL, maximo = 5, iter = 100,
        tol = 10^-8, norm = FALSE, contr = FALSE)

# Auxiliary functions for the previous ones:

read(X, Y = NULL)
colores(filas, columnas, coloresf, coloresc, conf, columnas2 = NULL, coloresc2 = NULL)
preproc(X, norm, cubo = FALSE, capas = FALSE)
pesos(t)
contributions(M, names, title, first = FALSE)
contributions2(M, dimension)
screeplot(d)
compr(dimX, dimY, d, c = NULL, tucker = FALSE)
plotm(dim, d, M1, M2, M3 = NULL, M4 = NULL, M5 = NULL, M6 = NULL, lim1, lim2, names1,
      names2, colores1, colores2, contf, contc, cotucker = FALSE)
tensorial(X, dimension, U)
desplegar(X, dimension)
results(X, title, names = NULL, axis = FALSE, first = FALSE)
every(p, q, r, T, iter, tol, mas = FALSE)
plotmt(dim1, dim2, dim3, sos, M1, M2, M3, lim1, lim2, lim3, names1, names2, names3,
       colores1, colores2, titles = FALSE)
}
\arguments{
\item{X,Y}{the two arrays with the numerical data, with or without labels. By default, if the labels are omitted,
the names 1, 2, ... are used for the first dimension, V1, V2, ... for the second dimension, or R1, R2, ... for the third dimension. They are similarly defined for the rest of the analysis.}
\item{p,q,r}{the number of components for each dimension to be retained. By default, if they are omitted, the analysis for all the combinations is computed.}
\item{P1,P2}{the pair of components for the first dimension to be plotted with the Tucker3 analysis from the Co-Tucker3 for the two arrays separately. By default, if they are omitted, these plots are not obtained.}
\item{Q1,Q2}{the pair of components for the second dimension to be plotted with the Tucker3 analysis from the Co-Tucker3 for the two arrays separately. By default, if they are omitted, these plots are not obtained.}
\item{R1,R2}{the pair of components for the third dimension to be plotted with the Tucker3 analysis from the Co-Tucker3 for the two arrays separately. By default, if they are omitted, these plots are not obtained.}
\item{dimAX,dimAY}{the axes for the first dimension to be plotted with the Co-Inertia analysis from the Co-Tucker3. By default, if they are omitted, these plots are not obtained.}
\item{dimBX,dimBY}{the axes for the second dimension to be plotted with the Co-Inertia analysis from the Co-Tucker3. By default, if they are omitted, these plots are not obtained.}
\item{dimCX,dimCY}{the axes for the third dimension to be plotted with the Co-Inertia analysis from the Co-Tucker3. By default, if they are omitted, these plots are not obtained.}
\item{coloresf}{a vector with the labels of the colours for the first dimension. If it is omitted, black labels are used.}
\item{coloresc1,coloresc2}{vectors with the labels of the colours for the second dimension of the two arrays. If they are omitted, black labels are used.}
\item{coloresc}{vector with the labels of the colours for the second dimension of one array or one matrix for the rest of the analysis. If it is omitted, black labels are used.}
\item{maximo}{the maximum number of components for each dimension to be retained. By default, if it is omitted, 5 is used.}
\item{iter}{the maximum number of iterations computed for each combination of components. By default, if it is omitted, 100 is used.}
\item{tol}{the minimum tolerance. By default, if it is omitted, 1e-08 is used.}
\item{norm}{a logical value indicating whether the data should be normalized by lateral slices. By default, if it is omitted, it is \code{FALSE}.}
\item{contr}{a logical value indicating whether the tables with the contributions for the first, second and third dimensions should be obtained. By default, if it is omitted, it is \code{FALSE}.}
\item{dimX,dimY,dimXx,dimYx,dimXy,dimYy}{the axis to be plotted with the rest of the analysis. By default, if they are omitted, they are not obtained.}
\item{gruposf}{a vector with the labels of the groups that the rows belong to. It has to have the same length as the number of rows.}
\item{cotucker}{a logical value indicating whether the Co-Inertia Analysis to be performed comes from the Co-Tucker3 analysis for the first ot the third dimensions. By default, if it is omitted, it is \code{FALSE}.}
\item{tcotucker}{a logical value indicating whether the Co-Inertia Analysis to be performed comes from the Co-Tucker3 analysis for the second dimension. By default, if it is omitted, it is \code{FALSE}.}

\item{filas, columnas, conf, columnas2}{arguments for the auxiliary function \code{colores}.}
\item{cubo, capas}{arguments for the auxiliary function \code{preproc}.}
\item{t}{argument for the auxiliary function \code{pesos}.}
\item{M, names, title, first}{arguments for the auxiliary functions \code{contributions, contributions2, results}.}
\item{dimension}{argument for the auxiliary functions \code{contributions2, tensorial, desplegar}.}
\item{d}{argument for the auxiliary functions \code{screeplot, compr, plotm}.}
\item{c, tucker}{arguments for the auxiliary function \code{compr}.}
\item{dim, M1, M2, M3, M4, M5, M6, lim1, lim2, names1, names2, colores1, colores2, contf, contc}{arguments for the auxiliary functions \code{plotm, plotmt}.}
\item{U}{argument for the auxiliary function \code{tensorial}.}
\item{axis}{argument for the auxiliary function \code{results}.}
\item{T, mas}{arguments for the auxiliary function \code{every}.}
\item{dim1, dim2, dim3, sos, lim3, names3, titles}{arguments for the auxiliary function \code{plotmt}.}
}
\details{
The number of elements of the first and third dimensions of both arrays must be the same.

The number of labels of the colours for the first and second dimensions of both arrays must be the same as the number of elements of the first and second dimensions of both arrays.

If there is any missing value in one row, this row is omitted.

\code{maximo,iter} must be integers greater than 1.

\code{p,q,r} must be integers greater than or equal to 1, and less than or equal to \code{maximo} and the number of elements of the corresponding dimension.

\code{P1,P2} must be integers greater than or equal to 1, less than or equal to the number of elements of the first dimension, and it must be \code{P1}<\code{P2} (exception: if only one component is retained). The same for \code{Q1,Q2} and \code{R1,R2}.

\code{dimAX,dimAY} must be integers greater than or equal to 1, less than or equal to the number of elements of the first dimension, and it must be \code{dimAX}<\code{dimAY}. The same for \code{dimBX,dimBY} and \code{dimCX,dimCY}.
}
\value{
If a combination of components has not been chosen, the result is a file with two tables: the one with all the combinations, and another with the combinations with a better fit. Moreover, a scree-plot, and a picture with the shape of the arrays are obtained.

If a combination of components has been chosen, the result is a file with three tables: the one with the percentages of fit for each component of each dimension, the core arrays, and the coordinates for each one of the dimensions. If \code{contr} is \code{TRUE}, the contributions of the three dimensions are obtained too.

If the axes (of the two individual Tucker3 analysis) have been chosen to be plotted, the graphics for the three dimensions of each array are obtained.

If the axes (of the two individual Tucker3 analysis) have not been chosen to be plotted, but the axes of the Co-Inertia Analysis have been chosen to be plotted, the three Co-Inertia Analysis are performed and the graphics with each one of the three dimensions according to both arrays are obtained.
}
\examples{
# X<-array(sample(0:1000,size=1408,replace=TRUE)/100,dim=c(22,16,4))
# A<-X[,1:10,]
# B<-X[,11:16,]
# COTUCKER3(A,B)
# COTUCKER3(A,B)
# COTUCKER3(A,B,p=3,q=3,r=2)
# COTUCKER3(A,B,p=3,q=3,r=2,P1=1,P2=2,Q1=1,Q2=2,R1=1,R2=2)
# COTUCKER3(A,B,p=3,q=3,r=2,dimAX=1,dimAY=2)
}
