## Read input data
args <- commandArgs(TRUE)
asc <- read.table(args[1], head = TRUE, comment.char = "")
asc <- na.omit(asc)

## Adjustment using genome-wide λGC
# Calculate chi-squared statistics from p-values
chi_squared <- qchisq(asc$p, df = 1, lower.tail = FALSE)
# Calculate the Genomic Inflation Factor (lambda GC)
lambda_gc <- median(chi_squared) / qchisq(0.5, df = 1)
print(paste("Genome-wide Genomic Inflation Factor (lambda GC):", lambda_gc))
# Adjust p-values using lambda GC
adjusted_chi_squared <- chi_squared / lambda_gc
adjusted_pvalues <- pchisq(adjusted_chi_squared, df = 1, lower.tail = FALSE)
# save new version of the input file 
asc$GWadjP=adjusted_pvalues

## Adjustment using λGC per chromosome
# Function to calculate lambda GC and adjust p-values per chromosome
adjust_lambda_gc_per_chromosome <- function(data) {
  # Calculate chi-squared statistics from p-values
  chi_squared <- qchisq(data$p, df = 1, lower.tail = FALSE)
  
  # Calculate the Genomic Inflation Factor (lambda GC)
  lambda_gc <- median(chi_squared) / qchisq(0.5, df = 1)
  print(paste("Genomic Inflation Factor (lambda GC) for chromosome", data$Chr[1], ":", lambda_gc))
    
  # Adjust p-values using lambda GC
  adjusted_chi_squared <- chi_squared / lambda_gc
  adjusted_pvalues <- pchisq(adjusted_chi_squared, df = 1, lower.tail = FALSE)
      
  # Add adjusted p-values to the data
  data$adjP <- adjusted_pvalues
  return(data)
}

# Apply the function to each chromosome
asc_adjusted <- do.call(rbind, lapply(split(asc, asc$Chr), adjust_lambda_gc_per_chromosome))

# Save the adjusted data to a file
write.table(asc_adjusted, paste(args[1],".adj",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,  col.names = TRUE)

##################
## Generate QQ plots
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
  qq(pvalues, main = "Q-Q Plot with 95% Confidence Interval", ylim = c(0, max(-log10(pvalues), na.rm = TRUE) + 0.5), cex.axis = 1.8, cex.lab = 1.8, cex = 1.5)

  # Add the confidence interval lines
  lines(-log10((1:n) / (n + 1)), -log10(ci_lower), col = "blue", lty = 2)
  lines(-log10((1:n) / (n + 1)), -log10(ci_upper), col = "blue", lty = 2)
}

tiff(paste(args[2],".unadj.qqplot.tif",sep=""),units="in", width=10, height=10, res=300);
par(mar = c(5, 5, 4, 2) + 0.1);
qq_plot_with_ci_qqman(asc$p);
dev.off();

tiff(paste(args[2],".genWideAdj.qqplot.tif",sep=""),units="in", width=10, height=10, res=300);
par(mar = c(5, 5, 4, 2) + 0.1);
qq_plot_with_ci_qqman(asc$GWadjP);
dev.off();

tiff(paste(args[2],".chrSPAdj.qqplot.tif",sep=""),units="in", width=10, height=10, res=300);
par(mar = c(5, 5, 4, 2) + 0.1);
qq_plot_with_ci_qqman(asc_adjusted$adjP);
dev.off();


