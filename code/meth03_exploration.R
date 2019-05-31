#'Load package for Age-Prediction
library(wateRmelon)
DNAmAge<-as.vector(agep(betas.clean))

#'Agreement between chronological age and DNAm-Age

#' Correlation
cor.test(DNAmAge,pheno$AGE)
plot(pheno$AGE,DNAmAge,pch=21,ylab="Age Predicted",
     xlab="Age Reported",cex=1.2, bg=alpha("deepskyblue",0.45),main="Prediction of Age")
legend("topleft",legend=paste0("r=",round(cor(DNAmAge,pheno$AGE),2)),bty="n")
abline(lm(DNAmAge~pheno$AGE),col="red",lw=2)


#' Age Acceleration Residuals
AgeAccelerationResidual<-residuals(lm(DNAmAge~pheno$AGE))
boxplot(AgeAccelerationResidual ~pheno$SMOKE_STATUS, col=c("red","blue"))
wilcox.test(AgeAccelerationResidual ~ pheno$SMOKE_STATUS)
t.test(AgeAccelerationResidual ~ pheno$SMOKE_STATUS)


#' Online Calculator
Horvath<-read.csv("C:/EBC3/Data/MethylationData.output.csv")
cor.test(Horvath$DNAmAge,pheno$AGE)
plot(pheno$AGE,Horvath$DNAmAge,pch=21,ylab="Age Predicted",
     xlab="Age Reported",cex=1.2, bg=alpha("deepskyblue",0.45),main="Prediction of Age")
legend("topleft",legend=paste0("r=",round(cor(Horvath$DNAmAge,pheno$AGE),2)),bty="n")
abline(lm(Horvath$DNAmAge~pheno$AGE),col="red",lw=2)
