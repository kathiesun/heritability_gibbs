inla.h2.results <- list()
for(i in 1:length(model.formula)){
  vceg <- list()
  vceg$y <- model.frame(formula(model.formula[i]), data=data.list[[data.use[i]]])[,1]
  vceg$R <- K[as.character(subjects.list[[data.use[i]]]), as.character(subjects.list[[data.use[i]]])]
  
  covariates <- model.matrix(formula(model.formula[i]), data=data.list[[data.use[i]]])
  #covariates <- as.vector(covariates[,-remove])
  vceg$x <- covariates
  
  inla.data <- list(prec.R = list(prior="loggamma", param=c(1,1)),
                    prec.e = list(prior="loggamma", param=c(1,1)),
                    inv.R = solve(vceg$R),
                    y = scale(as.numeric(vceg$y)),
                    x = vceg$x,
                    idx = 1:length(vceg$y))  
  
  inla.formula <- y ~ x + f(idx, model="generic2", constr=TRUE, Cmatrix=inv.R,
                            hyper=list(theta1=prec.R, theta2 = prec.e))
  
  inla.fit <- inla(formula=inla.formula, data=inla.data)
  hpd.ci <- inla.hpdmarginal(p=0.95, marginal=inla.fit$marginals.hyperpar[["h2 for idx"]])
  
  inla.h2.results[[i]] = cbind(inla.fit$summary.hyperpar["h2 for idx",], hpd.ci)
  names(inla.h2.results)[i] <- all.vars(formula(model.formula[i]))[1]
}
inla.vec <- unlist(inla.h2.results)
inla.mean <- round(inla.vec[grep("mean", names(inla.vec))], 3)
inla.median <- round(inla.vec[grep(".0.5quant", names(inla.vec), fixed=TRUE)], 3)
inla.mode <- round(inla.vec[grep("mode", names(inla.vec), fixed=TRUE)], 3)

inla.lb <- round(inla.vec[grep("0.025quant", names(inla.vec))], 3)
inla.ub <- round(inla.vec[grep("0.975quant", names(inla.vec))], 3)
inla.quant.ci <- paste("(", inla.lb, ", ", inla.ub, ")", sep="")

inla.hpd.lb <- round(inla.vec[grep("low", names(inla.vec))], 3)
inla.hpd.ub <- round(inla.vec[grep("high", names(inla.vec))], 3)
inla.hpd.ci <- paste("(", inla.hpd.lb, ", ", inla.hpd.ub, ")", sep="")

