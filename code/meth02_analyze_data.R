#'# Analyze methylation data  
#' Using data preprocessed in our script:  
#'  meth01_process_data.R  

#Elena Colicino

#' recall that we have a processed dataset with 15 samples
dim(WB)

suppressPackageStartupMessages({
  library(CpGassoc)
})

# Evaluate any remaining batch effects from design using the principal component analysis
PCobject = prcomp(na.omit(t(beta)), retx = T, center = T, scale. = T);
# Extract the Principal Components from SVD
PCs = PCobject$x;

# Extract the proportion of variability explained by the top R=2 PC.
R = 2
propvar = summary(PCobject)$importance["Proportion of Variance", 1:R]
sum(propvar)*100
cummvar = summary(PCobject)$importance["Cumulative Proportion", 1:R]
cummvar

X = 1  # First PC
Y = 2  # Second PC

# Save Graphic (Plate Effects)
variable = as.factor(a$Plate_ID)
myColors = rainbow(length(table(variable)))
names(myColors) = levels(variable)
# X-Y
plot(PCs[,X], PCs[,Y], xlab = paste("PC", X,"(",round(propvar[X]*100,1),"%)", sep =""),
     ylab = paste("PC", Y, "(",round(propvar[Y]*100,1), "%)", sep =""), pch=16,
     col = myColors)
legend("topleft", c("Plate 1","Plate 2"), pch=c(21,21),bty='n', cex=1,pt.bg = myColors, col=myColors)

## Extract sex prediction from methylation data
## use median total intensity of the X chromosome-mapped probes and
## the median total intensity of the Y-chromosome-mapped probes
## log2(med(X)) − log2(med(Y)) 
## cutoff of -2
## med(X)/med(y)=1

#'# Running an Epigenome Wide Association

#' Here we run an EWAS on sex (as a predictor of methylation)

#' and with adjustment for cell types

#' looking at results

#' qqplot and lambda interpretation

#' volcano plot interpretation

#' looking at a top hit



#' End of script