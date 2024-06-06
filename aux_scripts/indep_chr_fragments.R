# R script to read the GRM binary file and calculate the inverse of the variance of the off-diagonal elements
args=(commandArgs(TRUE));


ReadGRMBin <- function(prefix, AllN = F, size = 4) {
  sum_i <- function(i) {
    return(sum(1:i))
  }
  BinFileName <- paste(prefix, ".grm.bin", sep = "")
  NFileName <- paste(prefix, ".grm.N.bin", sep = "")
  IDFileName <- paste(prefix, ".grm.id", sep = "")
  
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  
  BinFile <- file(BinFileName, "rb")
  grm <- readBin(BinFile, n = n * (n + 1) / 2, what = numeric(0), size = size)
  close(BinFile)
  
  NFile <- file(NFileName, "rb")
  if (AllN == T) {
    N <- readBin(NFile, n = n * (n + 1) / 2, what = numeric(0), size = size)
  } else {
    N <- readBin(NFile, n = 1, what = numeric(0), size = size)
  }
  close(NFile)
  
  i <- sapply(1:n, sum_i)
  return(list(diag = grm[i], off = grm[-i], id = id, N = N))
}

# Function to calculate the inverse of the variance of the off-diagonal elements
calculate_inverse_variance <- function(grm_data) {
  off_diagonal_elements <- grm_data$off
  variance_off_diagonal <- var(off_diagonal_elements)
  inverse_variance_off_diagonal <- 1 / variance_off_diagonal
  return(inverse_variance_off_diagonal)
}

# usage
prefix <- args[1]  # GRM file prefix
grm_data <- ReadGRMBin(prefix)

inverse_variance <- calculate_inverse_variance(grm_data)
print(paste("Inverse of the variance of the off-diagonal elements:", inverse_variance))
print(paste("Bonferroni correction threathold (as a -log10) with the number of independent chromosome fragments:", -log10(0.05/inverse_variance)))
