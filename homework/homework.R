## `beta_matrix.csv` contains methylation data from 3 samples, one newborn and two nonagenarians.
## Which sample is the newborn? Use the Horvath epigenetic age prediction to find out.
## Coefficients of the Horvath model are provided in `13059_2013_3156_MOESM3_ESM.csv`.
## Use as last step the function below to transform the predictions from the Horvath model to chronological years.

inverse_transformation <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }