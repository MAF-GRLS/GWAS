# File: qq_plot_with_ci.R

# Load necessary libraries
#install.packages("qqman", repos = "http://cran.us.r-project.org")
library(qqman)

# Function to create Q-Q plot with 95% confidence interval using qqman
qq_plot_with_ci_qqman <- function(pvalues, conf_level = 0.95) {
  
  # Number of points
  n <- length(pvalues)
  
  # Theoretical quantiles
  theoretical_quantiles <- -log10((1:n) / (n + 1))
  
  # Observed quantiles
  observed_quantiles <- -log10(sort(pvalues))
  
  # Calculate the confidence intervals
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha / 2)
  ci_lower <- qbeta(alpha / 2, 1:n, n:1)
  ci_upper <- qbeta(1 - alpha / 2, 1:n, n:1)
  
  # Use qqman to create the Q-Q plot
  qq(pvalues, main = "Q-Q Plot with 95% Confidence Interval", ylim = c(0, max(-log10(pvalues), na.rm = TRUE) + 0.5))
  
  # Add the confidence interval lines
  lines(-log10((1:n) / (n + 1)), -log10(ci_lower), col = "blue", lty = 2)
  lines(-log10((1:n) / (n + 1)), -log10(ci_upper), col = "blue", lty = 2)
}


## Run for input data after adjustment using λGC
args=(commandArgs(TRUE));
asc <- read.table(args[1],head=TRUE,comment.char="");
asc <- na.omit(asc);
# Calculate chi-squared statistics from p-values
chi_squared <- qchisq(asc$P, df = 1, lower.tail = FALSE)
# Calculate the Genomic Inflation Factor (lambda GC)
lambda_gc <- median(chi_squared) / qchisq(0.5, df = 1)
print(paste("Genomic Inflation Factor (lambda GC):", lambda_gc))
# Adjust p-values using lambda GC
adjusted_chi_squared <- chi_squared / lambda_gc
adjusted_pvalues <- pchisq(adjusted_chi_squared, df = 1, lower.tail = FALSE)

# save new version of the input file 
asc$adjP=adjusted_pvalues
write.table(asc,paste(args[1],".adj",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

tiff(paste(args[2],".qqplot.tif",sep=""),units="in", width=10, height=10, res=300);
qq_plot_with_ci_qqman(adjusted_pvalues);
dev.off();


