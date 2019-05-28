library(ewastools)
library(stringi)
library(data.table)
library(magrittr)
library(purrr)

# Download phenotype data (copy URL in browser to see what the requested file looks like)
pheno = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85210&targ=gsm&form=text&view=brief"
pheno = readLines(pheno)

# Split into individual samples
pheno = split(pheno,cumsum(pheno %like% "^\\^SAMPLE = GSM"))

# Extract GSM accessions
names(pheno) = map(pheno,1) %>% stri_match_first(regex="GSM\\d+")

# Parse pheno data
imap(pheno,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> pheno

pheno = rbindlist(pheno)

# Keep only information on sample characteristics and supplementary files
pheno = pheno[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = pheno[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = pheno$value[i] %>% stri_split(fixed=": ")
pheno$variable[i] = map_chr(ch,1)
pheno$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
pheno[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
pheno[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

# Reshape data.table from long to wide format
pheno = dcast(pheno, gsm ~ variable)

# Select and parse the relevant variables
pheno = pheno[,.(gsm,smoker=factor(`subject status`),red,grn)]
pheno = rbind(pheno
	,list(
	 "GSM1185585"
	,"non-smoker"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1185nnn/GSM1185585/suppl/GSM1185585_6285625091_R04C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1185nnn/GSM1185585/suppl/GSM1185585_6285625091_R04C01_Grn.idat.gz")
	,list(
	 "GSM2219539"
	,"smoker"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2219nnn/GSM2219539/suppl/GSM2219539_6222421029_R02C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2219nnn/GSM2219539/suppl/GSM2219539_6222421029_R02C01_Grn.idat.gz")
	 )


# Select the samples
setkey(pheno,"gsm")
pheno = pheno[c(
	## first 15 smokers
 "GSM2260480","GSM2260482","GSM2260485","GSM2260486","GSM2260487"
,"GSM2260488","GSM2260489","GSM2260491","GSM2260493","GSM2260494"
,"GSM2260495","GSM2260496","GSM2260498","GSM2260499","GSM2260500"
	
	## first 15 non-smokers
,"GSM2260481","GSM2260483","GSM2260484","GSM2260490","GSM2260492"
,"GSM2260497","GSM2260501","GSM2260511","GSM2260514","GSM2260516"
,"GSM2260519","GSM2260525","GSM2260528","GSM2260530","GSM2260532"
	
,"GSM2260543" # same person as GSM2260485
,"GSM2260653" # this is the potentially contaminated sample
,"GSM1185585" # unrelated sample from another GSE, granulocytes instead of whole blood
,"GSM2219539" # unrelated sample of lung tissue
,"GSM2260573" # sample for which we'll change sex
)]

dir.create("data",showWarnings=FALSE)

# Download .idat files
map2(pheno$red,"data/" %s+% pheno$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(pheno$grn,"data/" %s+% pheno$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
pheno$red = NULL; pheno$grn = NULL

# Import the methylation data
meth = read_idats("data/" %s+% pheno$gsm)
pheno[,c("X","Y"):=check_sex(meth)]
pheno[,sex:=ifelse(X>1.,"f","m")]

pheno = pheno[,.(gsm,sex,smoker)]
pheno[gsm=="GSM2260573",sex:="f"]

write.csv(pheno,file="data/pheno.csv",row.names=FALSE)
write.csv(pheno[1:30],file="data/pheno_clean.csv",row.names=FALSE)
