# PRINCIPAL COMPONENT ANALYSIS

library("ggplot2")
#install.packages("ggfortify")
library("ggfortify")
library("gridExtra")
library("carData")
library("car")
library("factoextra")
library("corrplot")
wineQuality <- read.csv("C:/Users/Lenovo/Downloads/wineQualityReds.csv")
numerical_data <- wineQuality[,2:13]
head(numerical_data)
summary(numerical_data)
res.pca <- prcomp(numerical_data, scale = TRUE)
print(res.pca)
summary(res.pca)
eig.val<-get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca,addlabels = TRUE, col.var="blue")
var <- get_pca_var(res.pca)
var
head(var$cos2)
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
fviz_cos2(res.pca, choice = "var", axes = 1:2)
fviz_pca_var(res.pca,
             col.var = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE)
# Contributions of variables to PC1
a<-fviz_contrib(res.pca, choice = "var", axes = 1)
# Contributions of variables to PC2
b<-fviz_contrib(res.pca, choice = "var", axes = 2)
grid.arrange(a,b, ncol=2, top='Contribution of the variables to the first two PCs')

# FACTOR ANALYSIS
R=cor(numerical_data )
round(R,4)

R.smc <- (1 - 1 / diag(solve(R)))
diag(R) <- R.smc
round(R, 2)

r.eigen <- eigen(R)
r.eigen$values

library(nFactors)
ap <- parallel(subject=nrow(numerical_data),var=ncol(numerical_data),
               rep=100,cent=.05)
nS <- nScree(x=r.eigen$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

tot.prop <- 0
for (i in r.eigen$values) {
  tot.prop <- tot.prop + i / sum(r.eigen$values)
  print(tot.prop)
}

r.lambda <- as.matrix(r.eigen$vectors[,1:3]) %*% diag(sqrt(r.eigen$values[1:3]))
r.lambda

r.h2 <- rowSums(r.lambda^2)
r.u2 <- 1 - r.h2
com <- rowSums(r.lambda^2)^2 / rowSums(r.lambda^4)

cor.pa <- data.frame(cbind(round(r.lambda, 2), round(r.h2, 2), round(r.u2, 3), round(com, 1)))
colnames(cor.pa) <- c('PA1', 'PA2', 'PA3', 'h2', 'u2', 'com')
cor.pa


library(psych)
wq.cor.fa <- fa(numerical_data, nfactors = 3, rotate = 'none', fm = 'pa', max.iter = 1)
wq.cor.fa
wq.cor.fa.v <- fa(numerical_data, nfactors = 3, rotate = 'varimax', fm = 'pa', max.iter = 1)
wq.cor.fa.v
wq.cor.pa <- principal(numerical_data, nfactors = 3, rotate = 'varimax')
wq.cor.pa



# CANONICAL CORRELATION ANALYSIS
pole=read.csv("C:/Users/Lenovo/Downloads/wineQualityCCA.csv",header=TRUE)
acidity<-pole[,1:4]
head(acidity)

ide<-pole[,5:7]
head(ide)

library(GGally)
ggpairs(ide)
ggpairs(acidity)

library(lme4)
library(CCA) #facilitates canonical correlation analysis
library(CCP) #facilitates checking the significance of the canonical variates
#checking the between and within set associations
cormat<-matcor(ide,acidity)
#extracting the within study correlations for set 1 and set 2 and #between set cor
round(cormat$Ycor, 4)
round(cormat$Xcor, 4)
#between set associations
cormat<-matcor(ide,acidity)
round(cormat$XYcor, 4)

#obtaining the canonical correlations
can_cor1=cc(ide,acidity)
can_cor1$cor

#raw canonical coefficients
can_cor1[3:4]

#computes the canonical loadings
can_cor2=comput(ide,acidity,can_cor1)
can_cor2[3:6] #displays the canonical loadings

#test of canonical dimensions
rho=can_cor1$cor
##defining the number of observations, no of variables in first set,
#and number of variables in second set
n=dim(ide)[1]
p=length(ide)
q=length(acidity)
##Calculating the F approximations using different test statistics
#using wilks test statistic
p.asym(rho,n,p,q,tstat="Wilks")

#standardizing the first set of canonical coefficients(ide)
std_coef1<-diag(sqrt(diag(cov(ide))))
std_coef1%*%can_cor1$xcoef

##standardizing the coeficents of the second set (acidity)
std_coef2<-diag(sqrt(diag(cov(acidity))))
std_coef2%*%can_cor1$ycoef


