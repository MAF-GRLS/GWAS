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


