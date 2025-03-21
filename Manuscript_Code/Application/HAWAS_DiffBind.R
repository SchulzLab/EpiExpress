library("DiffBind")
library("tidyverse")

out_folder <- "HAWAS_DiffBind/"

run_tag <- "CLL_DiffBindInputNoCofac"
samples <- read.csv(paste0('*DiffBindInput_CLL_samples.txt'), sep='\t')

dbObj <- dba(sampleSheet=samples)

# Count number of reads in peaks.
dbObj <- dba.count(dbObj, summits=FALSE)
dbObj

# Decide on the comparison that will be made.
dbObj <- dba.contrast(dbObj, design="~Factor")
dbObj <- dba.contrast(dbObj, minMembers=2)

# I don't have the slightest idea when this line is needed and when not. 
# Sometimes it causes a crash when it's there, other times it ensures the run goes through.
# A guess is that it's only needed when using ChIP input.
options(mc.cores = 2)  # When any of these fail, good luck!
BiocParallel::register(BiocParallel::SerialParam())  # For the greylisting with the Input to work.

dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)#, bParallel=FALSE)  Comment in when not using the BiocParallel call above.

# Extract results.
deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast=1, th=1)
# Convert to dataframe to change GRanges into bed-coordinates.
out <- as.data.frame(deseq)

write.table(out, file=paste0(out_folder, run_tag, "_diffBind_Results.txt"), sep="\t", quote=F, row.names=F)

