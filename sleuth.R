library(sleuth)
library(data.table)
library(dplyr)

#upload the sample covariates file with the path, sample names, and conditions needed
stab <- read.table("sample_covariates.txt", header=TRUE, stringsAsFactors=FALSE, sep = ',')
print(stab)

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
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results (FDR/qval < 0.05)
sleuth_significant<- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
head(sleuth_significant)
#write top transcripts to file
sleuth_significant <- data.frame(sleuth_significant$target_id, sleuth_significant$test_stat, sleuth_significant$pval, sleuth_significant$qval)
colnames(sleuth_significant) <- c("target_id", "test_stat", "pval","qval")
sleuth_significant <- na.omit(sleuth_significant)

write.table(sleuth_significant, file="sleuth_output.txt",quote = FALSE,row.names = FALSE)
