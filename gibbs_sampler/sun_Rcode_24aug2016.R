setwd("/nas02/home/k/y/kys6")

library(mvtnorm) 
library(MASS)
library(coda)
library(mcmcplots)
library(matrixStats)
library(corrplot)
library(gplots)
library(Rmisc)

#################
# Gibbs sampler #
#################


gibbs <- function(n.sims, s2inv.0, t2inv.0, a, b, y, R_inv, burnin, thin, beta.0, X, sig2binv)
{
  n <- nrow(y)
  p <- ncol(X)
  I.p <- diag(p)
  I.n <- diag(n)
  
  # matrices and vectors to output at end of for loop
  u.draws <- matrix(NA, nrow = (n.sims - burnin)/thin, ncol=n)
  beta.draws <- matrix(NA, nrow = (n.sims - burnin)/thin, ncol=p)
  s2inv.draws <- c()
  t2inv.draws <- c()
  h2.draws <- c()
  
  s2inv.cur <- s2inv.0    #set arbitrary initial values for beta, s2inv, t2inv
  t2inv.cur <- t2inv.0
  beta.cur <- beta.0
  
  # co-diagonalization for faster matrix inverses
  eR <- eigen(R_inv, symmetric=TRUE)  # eigenvalues for spectral decomposition 
  U <- eR$vectors                     # store once for kinship matrix
  d <- eR$values
  eQ <- eigen(t(X)%*%X, symmetric=TRUE)
  W <- eQ$vectors
  f <- eQ$values
  
  # Iterative update functions   
  u.update <- function(X, beta, y, s2inv, t2inv, R_inv, n, U, d)
  {
    K <- ((t2inv*d+s2inv)^-0.5)*t(U)
    Kinv <- U*((t2inv*d+s2inv)^0.5)
    ztest <- rnorm(n)
    V_inv <-t(K)%*%K
    M <- (y-X%*%beta)*s2inv
    u.cur <- t(K)%*%ztest + V_inv%*%M
  }
  
  beta.update <- function(y, s2inv, u, sig2binv, X, p, f, W)
  { 
    K <- ((s2inv*f+sig2binv)^-0.5)*t(W)
    Kinv <- W*((s2inv*f+sig2binv)^0.5)
    V_inv <-t(K)%*%K
    ztest <- rnorm(p)
    M <- (t(X)%*%(y-u))*s2inv
    beta.cur <-  t(K)%*%ztest +V_inv%*%M
  }
  
  s2inv.update <-function(a, b, y, n, u, X, beta)
  {
    b.star <- 0.5*(t(y-X%*%beta-u)%*%(y-X%*%beta-u)) + b
    a.star <- a+(n/2)
    s2inv.cur <- rgamma(1, a.star, rate=b.star)
  }
  
  t2inv.update <-function(u, a, b, R_inv, n)
  {
    b.star <- 0.5*(t(u)%*%R_inv%*%u) + b
    a.star <- a+(n/2)
    t2inv.cur <- rgamma(1, a.star, rate=b.star)
  }
  
  h2.update <-function(t2inv, s2inv)
  {
    h2.cur <- s2inv/(t2inv+s2inv)
  }
  
  for (i in 1:n.sims) 
  {
    u.cur <- u.update(X=X, beta=beta.cur, y=y, R_inv=R_inv, s2inv=s2inv.cur, t2inv=t2inv.cur, n=n,  U=U, d=d)
    beta.cur <- beta.update(y=y, s2inv=s2inv.cur, u=u.cur, sig2binv=sig2binv, X=X,  p=p, W=W, f=f)
    s2inv.cur <- s2inv.update(X=X, beta=beta.cur, a=a, b=b, y=y, n=n, u=u.cur)
    t2inv.cur <- t2inv.update(u=u.cur, a=a, b=b, R_inv=R_inv, n=n)
    h2.cur <- h2.update(t2inv=t2inv.cur, s2inv=s2inv.cur)
    
    if (i > burnin & (i - burnin)%%thin == 0) 
    {
      j=(i - burnin)/thin
      u.draws [j,] <- u.cur
      beta.draws[j,] <- beta.cur
      s2inv.draws [j] <- s2inv.cur
      t2inv.draws [j] <- t2inv.cur
      h2.draws [j] <- h2.cur
      print(j)
    }
  }
  return(list(u.draws=u.draws, s2inv.draws=s2inv.draws, beta.draws=beta.draws,
              t2inv.draws=t2inv.draws, h2.draws=h2.draws)) 
}



######################
# DO 1 and 2 kinship #
######################


kinDO2 <- readRDS("kin_DO2.rds")
kin.nomiss <- kinDO2[c(-90, -181, -243, -195:-199), c(-90, -181, -243, -195:-199)]
kininv.nomiss <- solve(kin.nomiss)
  # DO2: remove samples with missing data

kinDO1 <- readRDS("kin_DO1.rds")
kin1.nomiss <- kinDO1[c(-15, -18, -242, -269, -295, -300, -307, -309, -310, -313), 
                      c(-15, -18, -242, -269, -295, -300, -307, -309, -310, -313)]
  kin1inv.nomiss <- solve(kin1.nomiss)
  # DO1: remove samples with missing data

###########
# DO data #
###########

phen <- read.csv("../DO_Pheno9.csv",header=TRUE)
phen$Diet[phen$Diet == 1] <- 9
phen$Diet[phen$Diet == 2] <- 0
phen$Diet[phen$Diet == 3] <- 1
  #remove diet 1 (control) from dataset so treatment is binary


phen.DO2 <- phen[-1:-316,]
phen.nomiss <- phen.DO2[c(-90, -181, -243, -195:-199),]
phen.DO1<- phen[1:315,]
phen1.nomiss <- phen.DO1[c(-15, -18, -242, -269, -295, -300, -307, -309, -310, -313), ]
  #remove individuals with missing data from phenotype dataset

#for DO1
myvars1 <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
             "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh", 
             "Sex")
#for DO2
myvars2 <-c("XREVS56", "XINTS56", "XRPM56", "XMAX56", "MassWk6", "MassWk12", "MassWk17", "PrMassCh", "Diet")

short_phen1 <- phen1.nomiss[myvars1]
short_phen2 <- phen.nomiss[myvars2]

##############################

#### starting from here, change inputs to reflect DO1 or DO2 depending on the cohort you want analyzed

n.real <-nrow(kin1.nomiss)
x_0.real <- rep(1, n.real)
x_1.real <- as.numeric(short_phen1$Sex) #sex as covariate
x.real <- as.matrix(cbind(x_0.real, x_1.real), ncol=2, byrow=FALSE)
p.real <-ncol(x.real)
b0.real <- rep(1, p.real)


#######################
# Run sampler for DO  #
#######################

sims <-5000
runs<-2
  #runs = number of chains for each set of y
pheno <-10
  #DO1 has 10 phenotypes, DO2 has 8
length <-350 
  #length = (sims-burnin) / thin

h2.end <- matrix(NA, nrow=runs*pheno, ncol=length)
s2inv.end <- matrix(NA, nrow=runs*pheno, ncol=length)
t2inv.end <- matrix(NA, nrow=runs*pheno, ncol=length)             
u.endlist <- list()

for (i in 1:pheno) 
{
  y.real <- as.matrix(short_phen1[,i], nrow=n.real)
  ysc.real <- scale(y.real, scale=FALSE)
  ymean.real <- mean(as.vector(y.real))
  yvar.real <- var(as.vector(y.real))
  
  for (j in 1:runs)
  {
    posterior.cur <- gibbs(n.sims = sims, s2inv.0=1, t2inv.0=1, a=1, b=yvar.real/2, 
                           y=ysc.real, beta.0=b0.real, R_inv=kin1inv.nomiss, burnin=1500, thin=10, X=x.real, sig2binv=1)
    num <- ((i-1)*runs)+j
    s2inv.end[num,]  <- as.vector(posterior.cur$s2inv.draws)
    t2inv.end[num,]  <- as.vector(posterior.cur$t2inv.draws)
    h2.end[num,] <- as.vector(posterior.cur$h2.draws)
    u.endlist [[num]] <- posterior.cur$u.draws
  }
}

saveRDS(s2inv.end, file="s2inv.1_12aug2016.rds")
saveRDS(t2inv.end, file="t2inv.1_12aug2016.rds")
saveRDS(h2.end, file="h2.1_12aug2016.rds")
saveRDS(u.endlist, file="u.1_12aug2016.rds")


## DO2 was saved as parameter.weight_10aug2016.rds files

##################################################
# Heterogeneous Stock kinship for simulated data #
##################################################

kinHS <- readRDS("kinship.Keele02052016.RDS")
kinhs_inv <- solve(kinHS)
kinhssm <- kinHS[1:300, 1:300]
kinshssm_inv <- solve(kinhssm)

########################
# Data for simulations # 
########################

n<-nrow(kinhssm)
p <- 1
x_1 <- rep(1, n)
x <- as.matrix(cbind(x_1), nrows=n)  
b0.0 <- rep(1, p)

##################################
# Run sampler for simulated data #
##################################

sims <-1500
runs<-5
length <-100

h2.test <- matrix(NA, nrow=runs, ncol=length)
#s2inv.test <- matrix(NA, nrow=runs, ncol=length)
#t2inv.test <- matrix(NA, nrow=runs, ncol=length)             
h2.actual <- numeric(runs)
#o2.actual <- numeric(runs)
#t2.actual <- numeric(runs)
#u.testlist <- list()
  
#more important to compare h2

for (i in 1:runs) 
{
  set.seed(abs(round(rnorm(1)*1000, digits=0)))
  
  bet<- rnorm(p)*100
  o2 <- abs(round(rnorm(1)*1000, digits=0)) # randomly generated sigma2 and tau2
  t2 <- abs(round(rnorm(1)*1000, digits=0))
  e <- scale(mvrnorm(1, rep(0, n), diag(n)*(o2)), scale=FALSE)
  u <- scale(mvrnorm(1, rep(0, n), kinhssm*(t2)), scale=FALSE)
  
  h2 <- t2/(o2+t2)
  h2.actual [i] <- h2
  #o2.actual [i] <- o2
  #t2.actual [i] <- t2
  y <- x%*%bet + u + e
  y.sc <- scale(y, scale=FALSE)
  var.y <- colVars(y)
  
  posterior.test <- gibbs(n.sims = sims, s2inv.0=1, t2inv.0=1, a=1, b=var.y/2, 
                          y=y.sc, beta.0=b0.0, R_inv=kinshssm_inv, burnin=500, thin=10, X=x, sig2binv=1)
  
  #s2inv.test [i, ] <- as.vector(posterior.test$s2inv.draws)
  #t2inv.test [i, ] <- as.vector(posterior.test$t2inv.draws)
  h2.test [i, ] <- as.vector(posterior.test$h2.draws)
  #u.testlist[[i]] <- posterior.test$u.draws
}

#plot concordance of estimated and actual h2
plot(rowMeans(h2.test), type='l', main='Comparison of heritability estimates from simulated data', 
     xlab='Run number', ylab="Heritability estimate")
mtext("n=300", 3)
lines(h2.actual, col='red')
legend("topright", inset=0.05, c("Estimated", "Actual"), text.col=c("red", "black"))


#############################################################
# Reverse normal transformation for other phenotypes in DO2 #
#############################################################

phenotype <- c("VO2_kg_mean",	"VO2_kg_mean_day", "VO2_lm_mean",	"VO2_lm_mean_day",	"VO2_mean",	"VO2_mean_day",
"VCO2_kg_mean",	"VCO2_kg_mean_day",	"VCO2_lm_mean",	"VCO2_lm_mean_day",	"VCO2_mean",	"VCO2_mean_day",	"RER_mean",
"RER_mean_day",	"Heat_kg_mean",	"Heat_kg_mean_day",	"Heat_lm_mean",	"Heat_lm_mean_day",	"Heat_mean",	"Heat_mean_day",	
"XT_YT_mean",	"XT_YT_mean_day",	"Drink_mean",	"Drink_mean_day",	"Feed_mean",	"Feed_mean_day",	"TSEDistK_mean",	
"TSEDistK_mean_day",	"TSEDistD_mean",	"TSEDistD_mean_day",	"VO2_kg_mean_night",	"VO2_lm_mean_night",	"VO2_mean_night",
"VCO2_kg_mean_night",	"VCO2_lm_mean_night",	"VCO2_mean_night",	"RER_mean_night",	"Heat_kg_mean_night",	"Heat_lm_mean_night",
"Heat_mean_night",	"XT_YT_mean_night",	"Drink_mean_night",	"Feed_mean_night",	"TSEDistK_mean_night",	"TSEDistD_mean_night",
"O2_mean_Day2",	"O2_mean_day_Day2",	"CO2_mean_Day2",	"CO2_mean_day_Day2",	"VO2_kg_mean_Day2",	"VO2_kg_mean_day_Day2",	
"VO2_lm_mean_Day2",	"VO2_lm_mean_day_Day2",	"VO2_mean_Day2",	"VO2_mean_day_Day2",	"VCO2_kg_mean_Day2",	"VCO2_kg_mean_day_Day2",
"VCO2_lm_mean_Day2",	"VCO2_lm_mean_day_Day2",	"VCO2_mean_Day2",	"VCO2_mean_day_Day2",	"RER_mean_Day2",	"RER_mean_day_Day2",	
"Heat_kg_mean_Day2",	"Heat_kg_mean_day_Day2",	"Heat_lm_mean_Day2",	"Heat_lm_mean_day_Day2",	"Heat_mean_Day2",	
"Heat_mean_day_Day2",	"XT_YT_mean_Day2",	"XT_YT_mean_day_Day2",	"Drink_mean_Day2",	"Drink_mean_day_Day2",	
"Feed_mean_Day2",	"Feed_mean_day_Day2",	"TSEDistK_mean_Day2",	"TSEDistK_mean_day_Day2",	"TSEDistD_mean_Day2",	
"TSEDistD_mean_day_Day2",	"VO2_kg_mean_night_Day2",	"VO2_lm_mean_night_Day2",	"VO2_mean_night_Day2",	"VCO2_kg_mean_night_Day2",
"VCO2_lm_mean_night_Day2",	"VCO2_mean_night_Day2",	"RER_mean_night_Day2",	"Heat_kg_mean_night_Day2",	"Heat_lm_mean_night_Day2",
"Heat_mean_night_Day2",	"XT_YT_mean_night_Day2",	"Drink_mean_night_Day2",	"Feed_mean_night_Day2",	"TSEDistK_mean_night_Day2",	
"TSEDistD_mean_night_Day2", "Diet")


phen.nomiss <- phen.DO2[c(-65,-195:-199),  ]
  #remove individuals with missing data
phen.rint <- phen.nomiss[phenotype]
kin.rint <- kinDO[c(-65,-195:-199),c(-65,-195:-199)]
  #remove individuals with missing data
kininv.rint <- solve(kin.rint)
n.real <- nrow(phen.rint)

# run rank based reverse normal transformation on 94 phenotypes
pheno <-length(phenotype)-1
phen.rinted <- matrix(NA, nrow=nrow(phen.rint), ncol=pheno)

for (i in 1:pheno) 
{
  phen.rinted[,i] <- qnorm((rank(phen.rint[,i])-0.5)/length(phen.rint[,i]))
}

sims <-5000
runs<-1
pheno <-94
length <-350


h2.end <- matrix(NA, nrow=runs*pheno, ncol=length)
s2inv.end <- matrix(NA, nrow=runs*pheno, ncol=length)
t2inv.end <- matrix(NA, nrow=runs*pheno, ncol=length)             
u.endlist <- list()

for (i in 1:pheno) 
{
  y.real <- as.matrix(phen.rinted[,i], nrow=n.real)
  ysc.real <- scale(y.real, scale=FALSE)
  ymean.real <- mean(as.vector(y.real))
  yvar.real <- var(as.vector(y.real))
  
  for (j in 1:runs)
  {
    posterior.cur <- gibbs(n.sims = sims, s2inv.0=1, t2inv.0=1, a=1, b=yvar.real/2, 
                           y=ysc.real, beta.0=b0, R_inv=kininv.rint, burnin=1500, thin=10, X=x.real, sig2binv=1)
    num <- ((i-1)*runs)+j
    s2inv.end[num,]  <- as.vector(posterior.cur$s2inv.draws)
    t2inv.end[num,]  <- as.vector(posterior.cur$t2inv.draws)
    h2.end[num,] <- as.vector(posterior.cur$h2.draws)
    u.endlist [[num]] <- posterior.cur$u.draws
  }
}

# saved for 94 phenotypes run on DO2 with diet as covariate as parameter.rint_10aug2016.rds files


#################
# Data analysis #
#################

# Step 1: concatenate chains together, e.g. 
# DO1 was run in duplicate, so pairs of output are added together

u1.1 <- rbind(u.1[[1]], u.1[[2]])
u1.2 <- rbind(u.1[[3]], u.1[[4]])
u1.3 <- rbind(u.1[[5]], u.1[[6]])
u1.4 <- rbind(u.1[[7]], u.1[[8]])
u1.5 <- rbind(u.1[[9]], u.1[[10]])
u1.6 <- rbind(u.1[[11]], u.1[[12]])
u1.7 <- rbind(u.1[[13]], u.1[[14]])
u1.8 <- rbind(u.1[[15]], u.1[[16]])
u1.9 <- rbind(u.1[[17]], u.1[[18]])
u1.10 <- rbind(u.1[[19]], u.1[[20]])



# Step 2: print highest posterior density regions, traceplots, and density plots for outputs using coda package
# this example code is for DO1 data; do something similar for DO2 and RINT transformed data

for (i in 1:pheno) 
{
  j <- i+(i-1)
  k <- j+1
  l <- i+8
  
  s2inv.mcmc <- as.mcmc.list(lapply(as.data.frame(t(s2inv.1[j:k,])), mcmc))
  t2inv.mcmc <- as.mcmc.list(lapply(as.data.frame(t(t2inv.1[j:k,])), mcmc))
  h2.mcmc <- as.mcmc.list(lapply(as.data.frame(t(h2.1[j:k,])), mcmc))
  
  #s2inv.mcmc <- as.mcmc.list(lapply(as.data.frame(t(s2inv.rint[i,])), mcmc))
  #t2inv.mcmc <- as.mcmc.list(lapply(as.data.frame(t(t2inv.rint[i,])), mcmc))
  #h2.mcmc[[i]] <- as.mcmc.list(lapply(as.data.frame(h2.fin1[i,]), mcmc))
  
  print(i)
  
  print(HPDinterval(h2.mcmc, prob=0.95))
  print(HPDinterval(s2inv.mcmc, prob=0.95))
  print(HPDinterval(t2inv.mcmc, prob=0.95))
  
  denoverplot1(s2inv.mcmc, xlab="Sampling for s2 inv")
  denoverplot1(t2inv.mcmc, xlab="Sampling for t2 inv")
  denoverplot1(h2.mcmc[[i]], xlab="Heritability", ylab="Density", main = "Estimated posterior distribution of heritability")
  
  traceplot(s2inv.mcmc, xlab="Iterations", ylab="Sampling for s2 inv")
  traceplot(t2inv.mcmc, xlab="Iterations", ylab="Sampling for t2 inv")
  traceplot(h2.mcmc, xlab="Iterations", ylab="Sampling for h2", main = "Traceplot of heritability estimates")
}

# Step 3: Getting correlations; this might not be the easiest way to get correlations between vectors

# Correlations between unique u vectors at each sample
for(i in 1:10)
{
  utemp <- u1.9
  distu.9 <- matrix(NA, nrow = 10, ncol=700)
  
  for (j in 1:700)
  {
    distu.9[1,j] <- cor(utemp[j,], u1.1[j,])
    distu.9[2,j] <- cor(utemp[j,], u1.2[j,])
    distu.9[3,j] <- cor(utemp[j,], u1.3[j,])
    distu.9[4,j ] <- cor(utemp[j,], u1.4[j,])
    distu.9[5,j ] <- cor(utemp[j,], u1.5[j,])
    distu.9[6,j ] <- cor(utemp[j,], u1.6[j,])
    distu.9[7,j ] <- cor(utemp[j,], u1.7[j,])
    distu.9[8,j ] <- cor(utemp[j,], u1.8[j,])
    distu.9[9,j ] <- cor(utemp[j,], u1.9[j,])
    distu.9[10,j ] <- cor(utemp[j,], u1.10[j,])
  }
}
    # run through distu.1, distu.2, ... distu.n as necessary to create nxn vectors
    # e.g. distu.9[6, 425] is the correlation between phenotype 9 and 6 at sample 425


# Correlations between mean u for each individual
for(i in 1:10)
{
  print("spearman")
  print(cor(colMeans(allu.list[[1]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[2]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[3]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[4]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[5]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[6]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[7]]), colMeans(allu.list[[i]]), method="spearman"))
  print(cor(colMeans(allu.list[[8]]), colMeans(allu.list[[i]]), method="spearman"))
}

# Visualize correlations
# unique u: pearson correlation heatmap
pearson1u <- as.matrix(read.csv("pearson1_u.csv",header=FALSE))
colnames(pearson1u) <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
                         "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh")
rownames(pearson1u) <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
                         "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh")

# mean u for each individual: spearman correlation heatmap
spearman1 <- as.matrix(read.csv("spearman1_1-10.csv",header=FALSE))
colnames(spearman1) <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
                         "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh")
rownames(spearman1) <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
                         "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh")
corrplot(pearson1u, method="color", col=col3(100))
corrplot(spearman1, method="color", col=col3(100))

