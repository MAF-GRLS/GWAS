args=(commandArgs(TRUE));

# Load necessary libraries
if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman", repos = "http://cran.us.r-project.org")
}
library(qqman)

# Read observed p-values
observed_pvalues <- read.table(args[2], header = TRUE)
observed_pvalues <- observed_pvalues[, c("SNP", "p")]

# Initialize a list to store permutation p-values
perm_pvalues <- list()

num_permutations <- args[1]
# Read permutation p-values
for (perm in 1:num_permutations) {
  for (chr in 13) { #for (chr in 1:38)
    perm_file <- paste0(args[3], ".chr", chr, ".", perm, ".mlma")
    perm_data <- read.table(perm_file, header = TRUE)
    perm_pvalues[[paste0("perm_", perm, "_chr_", chr)]] <- perm_data$p
  }
}

# Convert the list of permutation p-values to a matrix
perm_pvalues_matrix <- do.call(cbind, perm_pvalues)

# Calculate adjusted p-values
adjusted_pvalues <- observed_pvalues$p
for (i in 1:length(adjusted_pvalues)) {
  adjusted_pvalues[i] <- mean(perm_pvalues_matrix[i, ] <= observed_pvalues$p[i])
}

# Combine SNPs and adjusted p-values into a data frame
adjusted_pvalues_df <- data.frame(SNP = observed_pvalues$SNP, P = observed_pvalues$p, adj_P = adjusted_pvalues)

# Save adjusted p-values
write.table(adjusted_pvalues_df, paste0(args[2], ".adj"), row.names = FALSE, quote = FALSE)

# Try to calculate the inflation rate
adjusted_pvalues_dfx <- adjusted_pvalues_df[which(adjusted_pvalues_df$adj_P!=0),]
chi_squared <- qchisq(adjusted_pvalues_dfx$P, df = 1, lower.tail = FALSE)
chi_squared_adj <- qchisq(adjusted_pvalues_dfx$adj_P, df = 1, lower.tail = FALSE)
x=chi_squared_adj/chi_squared
x=x[x!=0 & !is.infinite(x)]
print(paste("mean(chi_squared_adj/chi_squared):",mean(x)))
print(paste("median(chi_squared_adj/chi_squared):",median(x)))


# Plot Q-Q plot of adjusted p-values
tiff(paste(args[2],".adj.qqplot.tif",sep=""),units="in", width=10, height=10, res=300);
qqman::qq(adjusted_pvalues_df$adj_P, main = "Q-Q Plot of Adjusted P-Values (Permutation Test)")
abline(0, 1, col = "red")
dev.off();


