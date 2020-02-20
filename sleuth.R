library(sleuth)
library(data.table)
library(dplyr)

stab <- read.table("sample_covariates.txt", header=TRUE, stringsAsFactors=FALSE, sep = ' ')
print(stab)
#colnames(stab) <- c("sample", "path")
so <- sleuth_prep(stab)
so
#fit a model comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')
so
#fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')
#perform the likelihood ratio test for differential expression between conditions
so<- sleuth_lrt(so, 'reduced', 'full')



#extract the test results from the sleuth object
sleuth_table<- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results (FDR/qval < 0.05)
sleuth_significant<- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
head(sleuth_significant)
head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10)
#write top 10 transcripts to file
write.table(sleuth_significant, file="significant.txt",quote = FALSE,row.names = FALSE)
