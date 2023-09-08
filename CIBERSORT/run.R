#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")

#"need.txt" input_matrix
#"ref.txt"  input_markers


source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "need.txt", perm=100, QN=TRUE)

