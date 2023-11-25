# bGWAS hogwash script 1.0.1 https://github.com/4C554D43
# 
# 
# Desc: Perform bGWAS utilizing Hogwash on SeqSphere 8.2.0 output data (.csv). 
# Supports Newick phylogenetic tree format, imported through treeio
# Hogwash: https://github.com/katiesaund/hogwash
# treeio:  https://github.com/YuLab-SMU/treeio
# EXPERIMENTAL BUILD. USE AT YOUR OWN RISK
# Version 1.0.0,  19-11-2023: Initial version
# Version 1.0.1,  25-11-2023: code optimization


# ----- CONFIG -----
#directory for input files. Raw string https://en.wikipedia.org/wiki/String_literal#Raw_strings
input_dir = r"(DirectoryPathHere)" 
filename_phenotype = "pheno.csv"
filename_genotype =  "geno.csv"
filename_phylotree = "tree.treefile"

output_dir = r"(DirectoryPathHere)"
output_name = "Hogwash_output" 
create_logfile = TRUE
create_dir = TRUE #create new folder for each run?
# ----- END CONFIG -----

#include and init
if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
if (!requireNamespace("readxl", quietly = TRUE)) {install.packages("readxl")}
if (!requireNamespace("dplyr", quietly = TRUE)) {install.packages("dplyr")}
if (!requireNamespace("hogwash", quietly = TRUE)) {devtools::install_github("4C554D43/hogwash")}
if (!requireNamespace("treeio", quietly = TRUE)) {devtools::install_github("4C554D43/treeio")}
library(hogwash)
library(readxl)
library(dplyr)
library(treeio)
# ----- functions -----
timedate_string <- function()
{
  current_datetime <- Sys.time() 
  formatted_datetime <- format(current_datetime, "%Y-%m-%d_%H-%M-%S")
  return (as.character(formatted_datetime))
}#end timedate
print_timestamped <- function(first, second = NULL, third = NULL, fourth = NULL, fifth = NULL, sixth = NULL)
{
  message <- paste(first, second, third, fourth, fifth, sixth, collapse = " ")
  if ("\n" %in% message) {
    cat(timedate_string(), ": ", message, "\n")
  } else {
    cat(timedate_string(), ": ", message, "\n", sep = "")
  }
}#end print
trim_columns <- function(mat) {
  keep_columns <- c()
  for (col_index in 1:ncol(mat)) 
  {
    if (length(unique(na.omit(mat[, col_index]))) > 1)    # Check if more than one unique non-missing value in column
      keep_columns <- c(keep_columns, col_index)
  }
  removed_columns <- ncol(mat) - length(keep_columns) # Determine the number of columns removed
  result_matrix <- mat[, keep_columns, drop = FALSE]  # Subset the matrix to keep only the selected columns
  #print and return
  print_timestamped("trim_column(): Columns removed: ", removed_columns)
  print_timestamped("trim_column(): Columns kept:", length(keep_columns))
  return(result_matrix)
}#end trim

binary_encode_matrix <- function(mat) { #assign most common value 1, others with 0. Missing assigned 0
  for (col_index in 1:ncol(mat)) 
  {
    col_values <- mat[, col_index]
    value_counts <- table(na.omit(col_values))
    most_common_value <- as.numeric(names(value_counts[which.max(value_counts)]))
    mat[, col_index] <- ifelse(is.na(col_values), 1, ifelse(col_values == most_common_value, 1, 0))
  }
  return(mat)
} #end binary encode
print_phenostats <- function(mat)#prints phenotype input matrix stats
{
  num_true = 0
  num_false = 0
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      # Print the value and its position
      if (mat[i, j] == 1){
        num_true = num_true + 1
      }else{num_false = num_false + 1}
      print_timestamped("Value at position (", i, ",", j, "): ", mat[i, j])
    }
  }
  print_timestamped("Phenotype samples: ", nrow(severity_matrix), "; severe: ", num_true, "; mild: ", num_false)
}#end phenostats

# ----- main -----
#resolve I/O
if (create_dir == TRUE)
{
  output_dir = paste(output_dir, timedate_string(), "_Hogwash")
  dir.create(output_dir)
  if (!dir.exists(output_dir))
  {
    print_timestamped("[ERR] Could not create output directory")
    stop()
  }
}
sink(file.path(output_dir, r"(log.txt)"), split = TRUE, append = TRUE)
Sys.Date()
Sys.time()
R.Version()

path_phenotype = file.path(input_dir, filename_phenotype) #sample severity matrix path) 
path_genotype = file.path(input_dir, filename_genotype) #sample data matrix path
path_phylotree = file.path(input_dir, filename_phylotree)
if (!file.exists((path_phenotype))|| !file.exists((path_genotype))|| !file.exists((path_phylotree)))
{
  print_timestamped("[ERR] An input file (severity, sample data or tree) was not found!")
  stop()
}
phenotype_matrix <- as.matrix(read.csv(path_phenotype, header = TRUE, row.names = 1, sep=';'))
genotype_matrix <- as.matrix(read.csv(path_genotype,header = TRUE, row.names = 1, sep=';'))
phylotree <- read.tree(path_phylotree)


#Trim alleles with same name, and convert to binary for Hogwash
genotype_matrix_trimmed = trim_columns(genotype_matrix)
genotype_matrix_binary = binary_encode_matrix(genotype_matrix_trimmed)
print_timestamped("Phylo tree leaves: ", ape::Ntip(phylotree))


#print stats
print_phenostats(phenotype_matrix)
print_timestamped("Genotype samples: ", nrow(genotype_matrix_binary))
print_timestamped("Genotype loci: ", ncol(genotype_matrix_binary))


cat("Calling hogwash")
#main Hogwash fnc
hogwash(pheno = phenotype_matrix, 
        geno = genotype_matrix_binary,  
        tree = phylotree,      
        tree_type = "phylogram",
        file_name = output_name,
        dir = output_dir,
        perm = 10000,
        fdr = 0.15,
        bootstrap = 0.70,
        group_genotype_key = NULL, 
        grouping_method = "post-ar",
        test = "both")


print_timestamped("Hogwash done")
sink()



