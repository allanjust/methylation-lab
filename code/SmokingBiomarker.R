#'###########
#'LASSO on model 2: beta values, fully adjusted
#' select most significant 1000 hits
betas.restricted <- t(betas.clean[rownames(betas.clean)%in%rownames(results2.tab[1:1000,]),])
#' run a lasso and plot the lambda parameter
fits.lasso <- cv.glmnet(betas.restricted, pheno[,1], alpha=1, family="binomial")
plot(fits.lasso)

# we perform variable selection based on the smallest CVM
small.lambda.index <- which(fits.lasso$lambda == fits.lasso$lambda.min)
small.lambda.betas <- fits.lasso$glmnet.fit$beta[,small.lambda.index]
#' select the hits and check how many are not null
selectedVariableIndex<-(abs(small.lambda.betas)>0)
beta_nonzero=small.lambda.betas[selectedVariableIndex]
length(beta_nonzero)
colnames(betas.restricted)[selectedVariableIndex]
#' coefficients of the lasso
Coeff.lasso<-coef(fits.lasso, s = fits.lasso$lambda.min)
#' 
#' #predict smoking
pheno$predicted.val = predict(fits.lasso, betas.restricted, s=0.001, type="response")
plot(predicted.val, pheno[,1])
