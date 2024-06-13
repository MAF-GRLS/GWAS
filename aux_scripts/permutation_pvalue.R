args=(commandArgs(TRUE));

## Initialize a list to store permutation p-values
perm_pvalues <- list()

num_permutations <- args[1] ## 1000
# Read permutation p-values
for (perm in 1:num_permutations) {
  for (chr in args[2]) { #for (chr in 1:38)
    perm_file <- paste0(args[3], ".chr", chr, ".", perm, ".mlma")
    perm_data <- read.table(perm_file, header = TRUE)
    perm_pvalues[[paste0("perm_", perm, "_chr_", chr)]] <- perm_data$p
  }
}

# Convert the list of permutation p-values to a matrix
perm_pvalues_matrix <- do.call(cbind, perm_pvalues)

# calc the maximum –log10 (p‐value) of each permutation (i.e. the smallest permuted p-value per permutation)
perm_pvalues_matrix <- -log10(perm_pvalues_matrix)
perm_pvalues_max <- apply(perm_pvalues_matrix, 2, max)

# define the 95th percentile of the ordered recorded values.
perm_pvalues_max_ordered <- sort(perm_pvalues_max)
percentile_95 <- quantile(perm_pvalues_max_ordered, 0.95)
#print(paste("The 95th percentile of the max –log10 (p‐value) of ",args[1], " permutations: ", percentile_95))

cat(percentile_95)
