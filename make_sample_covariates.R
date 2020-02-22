
# make the sample covariates file for sleuth
# sleuth object requires a sample covariates file with the sample name and the path to that folder
# the sample covariates file includes the four SRR folders

args <- commandArgs(trailingOnly=TRUE)
print(args)

make_sample_to_covariates <- function(x, x1, x2, x3){
  sample <- c(x, x1, x2, x3)
  path <- c(paste('./', x, '/', sep =''), paste('./', x1, '/', sep =''), paste('./', x2, '/', sep =''), paste('./', x3, '/', sep =''))
  #path <- c(x, x1, x2, x3)
  condition <- c("2dpi", "6dpi", "2dpi", "6dpi")
  data <- data.frame(sample=sample, path=path, condition=condition)
  write.table(data, 'sample_covariates.txt', quote=F, row.names = F, col.names = T, sep = ',')
}

make_sample_to_covariates(args[1], args[2], args[3], args[4])
