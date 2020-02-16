library(sleuth)
library(data.table)
library(dplyr)

stab <- read.table("sample_covariates.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(stab) <- c("path", "sample")
so <- sleuth_prep(stab)
so

so <- sleuth_lrt(so, 'reduced', 'full')


#extract the test results from the sleuth object
sleuth_table<- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results (FDR/qval < 0.05)
sleuth_significant<- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
head(sleuth_significant)
head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10)
#write top 10 transcripts to file
write.table(sleuth_significant, file="topten.txt",quote = FALSE,row.names = FALSE)


#first extract needed results from kallisto for plotting
#so <- sleuth_prep(stab, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)
