#!/usr/bin/env Rscript

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2); packageVersion("ggplot2")
   }
if (!require("xlsx")) {
   install.packages("xlsx", dependencies = TRUE)
   library(xlsx)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
  }
if (!require("psych")) {
   install.packages("psych", dependencies = TRUE)
   library(psych)
   }
 if (!require("optparse")) {
   install.packages("optparse", dependencies = TRUE)
   library(optparse)
   }

option_list = list(
   make_option(c("-i", "--Input_path"), type="character", default=NULL,
              help="path to fasta.otu.txt", metavar="character"),
   make_option(c("-o", "--Output_path"), type="character", default="Figs",
              help="error corretion model", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$Input_path)){
  print_help(opt_parser)
  stop("Input argument(s) missing.n", call=FALSE)
}

RootPath <- file.path(opt$Input_path) # Change to the path, where the result from DADA2 pipeline or similar pipelines is stored.
FigsPath <- file.path(opt$Output_path) # Path where output will be stored.

print("RootPath")
print(RootPath)

blast_result_file = (list.files(path = FigsPath, pattern ="output_blast_results",full.names=TRUE))

blast_results_db <- read.table( file=blast_result_file) # Results from BLAST analysis above.
colnames(blast_results_db) <- c( "Query ID", "Subject", "Identity percentage",
"Coverage", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
"S.start", "S.end", "Evalue", "Bitscore" ) # Matches ScandiFish and ScandiMarMal build-up

summary(blast_results_db) # Checks if everything works
dim(blast_results_db)
length(unique(blast_results_db$Subject))
length(unique(blast_results_db$"Query ID"))

tax_db = stringr::str_replace_all(blast_results_db$Subject, ";", " / ") # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis
taxonomy = gsub(".*[|]", "", tax_db) # Removes "|" from header, replaces with nothing, to make compatible with downstream analysis

blast_results_db$"Taxonomy" = taxonomy

tax_acc = unique(blast_results_db$Subject)
tax_accno = gsub(".*[|]([^|]+)[|].*", "\\1", tax_acc) # Extracts the accession number from subject id
tax_accno = data.frame(cbind(tax_acc, tax_accno))
colnames(tax_accno) = c("Subject", "Subject accession") # Stores accession number and subject ID in same column.

blast_results_db = left_join(blast_results_db, tax_accno, by = "Subject") # Joins blast_results by accession and subject ID.
blast_results_db$"Subject Taxonomy ID" = blast_results_db$"Subject" # Data formating
blast_results_db$"Source" = "ScandiFish"

blast_results_db_final = data.frame(cbind(blast_results_db$"Query ID",
                                          blast_results_db$"Subject",
                                          blast_results_db$"Subject accession",
                                          blast_results_db$"Subject Taxonomy ID",
                                          blast_results_db$"Identity percentage",
                                          blast_results_db$"Coverage",
                                          blast_results_db$"Evalue",
                                          blast_results_db$"Bitscore",
                                          blast_results_db$"Source",
                                          blast_results_db$"Taxonomy"
                                          ))

colnames(blast_results_db_final) = c("#Query ID",
                                     "#Subject",
                                     "#Subject accession",
                                     "#Subject Taxonomy ID",
                                     "#Identity percentage",
                                     "#Coverage",
                                     "#Evalue",
                                     "#Bitscore",
                                     "#Source",
                                     "#Taxonomy"
                                          )
print(blast_results_db_final)
output_name = paste0(basename(blast_result_file),".lca.tabular")
write.table(blast_results_db_final,file.path(FigsPath, file=output_name), row.names = FALSE, quote = FALSE, sep = "\t") # Stores the BLAST results in a table accepted by the LCA analysis.

