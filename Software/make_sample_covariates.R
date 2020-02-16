

args <- commandArgs(trailingOnly=TRUE)
print(args)

make_sample_to_covariates <- function(x, x1, x2, x3){
  sample <- c(x, x1, x2, x3)
  path <- c(paste("./", x, sep =''), paste("./", x1, sep =''), paste("./", x2, sep =''), paste("./", x3, sep =''))
  data <- data.frame(sample=sample, path=path)
  write.table(data, 'sample_covariates.txt', quote=F, row.names = F, col.names = T)
}

make_sample_to_covariates(args[1], args[2], args[3], args[4])
