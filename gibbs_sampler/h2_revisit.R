
#Mac
setwd("~/Google Drive/Summer rotation/KSun final docs") 

#PC
setwd("C:/Users/kathie/Google Drive/Summer rotation/KSun final docs/")

#source("../../Dropbox/source/nutridiallel2014/src/plot.hpd.R")

library(dglm)
library(varComp)
library(INLA)
library(rstan)
library(coda)

# DO 1 and 2 kinship #
kindo2 <- readRDS("./data/kin_DO2.rds")
s="/home/tangi/scratch/pompDO2miceonly//pomp2DP."
colnames(kindo2) <- gsub(s,"",colnames(kindo2))
rownames(kindo2) <- colnames(kindo2)

kindo1 <- readRDS("./data/kin_DO1.rds")
s="/home/tangi/scratch/pompDO1miceonly//pomp1DP."
colnames(kindo1) <- gsub(s,"",colnames(kindo1))
rownames(kindo1) <- colnames(kindo1)

kinmat <- list()
kinmat[[1]] <- kindo1
kinmat[[2]] <- kindo2

# DO data #
phen <- read.csv("./data/DO_Pheno9.csv",header=TRUE)
phen$Diet[phen$Diet == 1] <- 9
phen$Diet[phen$Diet == 2] <- 0
phen$Diet[phen$Diet == 3] <- 1
phen$Sex <- ifelse(phen$Sex == 0, 2, 1)
#remove diet 1 (control) from dataset so treatment is binary

#for DO1
allvars <- list()
allvars[[1]] <- c("XRPM56", "XINTS56", "XREVS56", "XMAX56", 
              "MassWk6", "MassWk7", "MassPost", "PrMassCh", "PrFatPost", "PrFatCh", 
              "Sex", "MouseID")
#for DO2
allvars[[2]] <-c("XREVS56", "XINTS56", "XRPM56", "XMAX56", "MassWk6", "MassWk12", "MassWk17", "PrMassCh", "Diet","MouseID")

colnames(phen)[1] <- "MouseID"

compl_phen <- list()
compl_kin <- list()
for(i in 1:2){
  compl_phen[[i]] <- phen[which(phen$DO == i), allvars[[i]]]
  #miss <- which(rowSums(is.na(phen_temp[,allvars[[i]]])) > 0)
  #mouseID_miss <- phen_temp$MouseID[miss]
  #remove <- which(as.numeric(unlist(strsplit(colnames(kinmat[[i]]),"[.]"))[c(F,T,F)]) %in% mouseID_miss)
  remove <- setdiff(compl_phen[[i]]$MouseID, as.integer(unlist(strsplit(colnames(kinmat[[i]]),"[.]"))[c(F,T,F)]))
  if(length(remove) > 0){
    if(remove %in% as.integer(unlist(strsplit(colnames(kinmat[[i]]),"[.]"))[c(F,T,F)])){
      compl_kin[[i]] <- kinmat[[i]][-remove, -remove]
    } else {
      compl_phen[[i]] <- compl_phen[[i]][-which(compl_phen[[i]]$MouseID == remove),]
      compl_kin[[i]] <- kinmat[[i]]
    }
  } else {
    compl_kin[[i]] <- kinmat[[i]]
  }
  #compl_phen[[i]] <- phen_temp[-which(phen_temp$MouseID %in% mouseID_miss),]
}

##########
#  STAN  #
##########
stan_file <- file.path("../../..","Dropbox", "doubleGLM.stan")
#dissector.sm <- stan_model(file=stan_file)
stan_fit <- list()
vceg <- list()

for(j in 1:2){
  covar <- ifelse(j==1, "Sex", "Diet")

  for(i in 1:7){ 
    miss <- which(is.na(compl_phen[[j]][,allvars[[j]][i]]))
    mouseID_miss <- compl_phen[[j]]$MouseID[miss]
    remove <- which(as.numeric(unlist(strsplit(colnames(compl_kin[[j]]),"[.]"))[c(F,T,F)]) %in% mouseID_miss)
    if (length(mouseID_miss) > 0){
      temp_phen <- compl_phen[[j]][-which(compl_phen[[j]]$MouseID %in% mouseID_miss),]
      temp_kin <- compl_kin[[j]][-remove, -remove]
    } else {
      temp_phen <- data.frame(compl_phen[[j]])
      temp_kin <- compl_kin[[j]]
    }
    
    vceg$x <- temp_phen[,covar]
    vceg$R <- temp_kin
    vceg$Rinv <- solve(vceg$R)
    vceg$y <- as.vector(scale(as.numeric(temp_phen[,allvars[[j]][i]])))
    stan_fit <- stan(file = stan_file, control = list(adapt_delta = 0.8), 
                     data = list(num_cov = length(unique(vceg$x)),
                                                   N = nrow(temp_kin),
                                                   phenotype = vceg$y,
                                                   cov = vceg$x,
                                                   R = temp_kin), chains=3, iter = 5000, 
                                                   thin=10, warmup = 1000)
  }
}

#h2_mcmc <- as.mcmc(stan_fit@sim$samples[[1]]$h2)
#mcmc_trace(regex_pars = 'sigma')
rstan::traceplot(stan_fit, pars = 'h2', inc_warmup=F)
rstan::traceplot(stan_fit, pars = 'grand_sig2', inc_warmup=F)
rstan::traceplot(stan_fit, pars = 'cov_vef', inc_warmup=F)
######TO-DO stats::acf(stan_fit@sim$samples[[1]]$h2)


###########
# Varcomp #
###########
tau2 <- list()
sig2 <- list()
h2 <- list()

for(j in 1:2){
  tau2[[j]] <- vector()
  sig2[[j]] <- vector()
  h2[[j]] <- vector()
  covar <- ifelse(j==1, "Sex", "Diet")
  for (i in 1:(length(allvars[[j]])-2)){
    miss <- which(is.na(compl_phen[[j]][,allvars[[j]][i]]))
    mouseID_miss <- compl_phen[[j]]$MouseID[miss]
    remove <- which(as.numeric(unlist(strsplit(colnames(compl_kin[[j]]),"[.]"))[c(F,T,F)]) %in% mouseID_miss)
    if (length(mouseID_miss) > 0){
      temp_phen <- compl_phen[[j]][-which(compl_phen[[j]]$MouseID %in% mouseID_miss),]
      temp_kin <- compl_kin[[j]][-remove, -remove]
    } else {
      temp_phen <- data.frame(compl_phen[[j]])
      temp_kin <- compl_kin[[j]]
    }
    use_form <- formula(paste("Y ~ 1 +", covar))
    temp_phen$Y <- scale(temp_phen[,allvars[[j]][i]])
    test <- varComp(use_form, data=temp_phen, varcov=list("R"=temp_kin))
    tau2[[j]][i] <- test$varComps
    sig2[[j]][i] <- test$sigma2
    h2[[j]][i] <- tau2[[j]][i]/(sig2[[j]][i]+tau2[[j]][i])
  }
  names(h2[[j]]) <- allvars[[j]][1:length(h2[[j]])]
}

########
# INLA #
########
inla.h2.results <- list()
inla.mat <- list()
vceg <- list()
inla.fit <- list()
prior="beta"

for(j in 1:2){
  covar <- ifelse(j==1, "Sex", "Diet")
  inla.h2.results[[paste0("do",j,"_",prior)]] <- list()
  
  for(i in 1:3){   #(length(allvars[[j]])-2)
    miss <- which(is.na(compl_phen[[j]][,allvars[[j]][i]]))
    mouseID_miss <- compl_phen[[j]]$MouseID[miss]
    remove <- which(as.numeric(unlist(strsplit(colnames(compl_kin[[j]]),"[.]"))[c(F,T,F)]) %in% mouseID_miss)
    if (length(mouseID_miss) > 0){
      temp_phen <- compl_phen[[j]][-which(compl_phen[[j]]$MouseID %in% mouseID_miss),]
      temp_kin <- compl_kin[[j]][-remove, -remove]
    } else {
      temp_phen <- data.frame(compl_phen[[j]])
      temp_kin <- compl_kin[[j]]
    }
    
    vceg$x <- temp_phen[,covar]
    vceg$R <- temp_kin
    vceg$Rinv <- solve(vceg$R)
    vceg$y <- scale(as.numeric(temp_phen[,allvars[[j]][i]]))
    a <- ifelse(prior=="gam",0.01, 1)
    b <- ifelse(prior=="gam",0.01, var(vceg$y)/2)
    
    inla.data <- list(prec.R = list(prior="loggamma", param=c(a,b)),
                      prec.e = list(prior="loggamma", param=c(a,b)),
                      inv.R = vceg$Rinv,
                      y = vceg$y,
                      x = vceg$x,
                      idx = 1:length(vceg$y))  
    
    inla.formula <- y ~ x + f(idx, model="generic2", constr=TRUE, Cmatrix=inv.R,
                              hyper=list(theta1=prec.R, theta2 = prec.e))
    
    inla.fit[[paste0("do",j,"_",prior)]][[i]] <- inla(formula=inla.formula, data=inla.data)
    temp <- inla.fit[[paste0("do",j,"_",prior)]][[i]]
    if(!is.null(temp$marginals.hyperpar)){
      hpd.ci <- inla.hpdmarginal(p=0.95, marginal=temp$marginals.hyperpar[["h2 for idx"]])
      inla.h2.results[[paste0("do",j,"_",prior)]][[i]] = cbind(temp$summary.hyperpar["h2 for idx",], hpd.ci)
      rownames(inla.h2.results[[paste0("do",j,"_",prior)]][[i]]) <- allvars[[j]][i]
    }
  }
  inla.mat[[paste0(prior, j)]] <- do.call(rbind, inla.h2.results[[paste0("do",j,"_",prior)]])
}

#########################



