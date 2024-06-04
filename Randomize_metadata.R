library(dplyr)

# Load metadata

map_bac <- read.delim('16S_2801_updated_mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

map_vir <- read.delim('virome_cleaned_2801_mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

################Bacteria

df <- map_bac[, -c(1:4)]

# Columns to not randomize
keep.cols <- c("Days","Bronchopulmonary_dysplasia", "Intrapartum_antibiotic_prophylaxis")

columns <- colnames(map_bac)
columns_to_randomize <- columns[!columns %in% keep.cols]

set.seed(123)  # For reproducibility
# Replace values in specified columns with appropriate random values
for (col_name in columns_to_randomize) {
  if (is.numeric(df[[col_name]])) {
    # Numeric columns: Replace with random values from the same distribution
    mean_val <- mean(df[[col_name]], na.rm = TRUE)
    sd_val <- sd(df[[col_name]], na.rm = TRUE)
    df[[col_name]] <- rnorm(n = nrow(df), mean = mean_val, sd = sd_val)
  } else if (is.character(df[[col_name]])) {
    # Character columns: Shuffle the values
    df[[col_name]] <- sample(df[[col_name]])
  }
}

readr::write_tsv(df,"mapping_bacteria_randomized.txt")

################Virome

df <- map_vir[, -c(1:2)]

# Columns to not randomize
keep.cols <- c("Days","Bronchopulmonary_dysplasia", "Intrapartum_antibiotic_prophylaxis")

columns <- colnames(map_bac)
columns_to_randomize <- columns[!columns %in% keep.cols]

set.seed(123)  # For reproducibility
# Replace values in specified columns with appropriate random values
for (col_name in columns_to_randomize) {
  if (is.numeric(df[[col_name]])) {
    # Numeric columns: Replace with random values from the same distribution
    mean_val <- mean(df[[col_name]], na.rm = TRUE)
    sd_val <- sd(df[[col_name]], na.rm = TRUE)
    df[[col_name]] <- rnorm(n = nrow(df), mean = mean_val, sd = sd_val)
  } else if (is.character(df[[col_name]])) {
    # Character columns: Shuffle the values
    df[[col_name]] <- sample(df[[col_name]])
  }
}

readr::write_tsv(df,"mapping_virome_randomized.txt")
