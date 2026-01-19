# Load necessary library
library(tidyverse)
library(here)
library(conflicted)

# Force R to use dplyr for select and filter
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load the annotated data 
# Using read_csv from the readr package (part of tidyverse)
# Read using base R, specifying the first column as row names
data <- read.csv("analysis_ready_count_data.csv", row.names = 1, check.names = FALSE)



