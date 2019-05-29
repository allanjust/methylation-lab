# The file `beta_matrix.csv` contains methylation data from 3 samples (A,B,C)
# One sample belong to a newborn, the others to nonagenarians.
# Your task is to identify the newborn sample?
# Use the Horvath model to predict the epigenetic age.
# The coefficients of the Horvath model are provided in the file `13059_2013_3156_MOESM3_ESM.csv`.
# Model coefficients are in the second column.
# Multiply the model coefficients with the methylation levels of the corresponding CpG site (353 biomarkers), summarize 
# Use as last step the function below to transform the predictions from the Horvath model to chronological years.
# The task can be completed with base R, but you can of course use any package you want.

inverse_transformation <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }