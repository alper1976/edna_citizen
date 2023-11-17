#!/usr/bin/env Rscript

if (!require("dada2")) {
   BiocManager::install("dada2")
   library(dada2); packageVersion("dada2")
   }

if (!require("phyloseq")) {
   BiocManager::install("phyloseq", dependencies = TRUE)
   library(phyloseq); packageVersion("phyloseq")
   }

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
if (!require("limma")) {
   BiocManager::install("limma")
   library(limma)
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
  make_option(c("-i", "--Input_path"), type="character", default="/cluster/projects/nn9745k/02_results/26_lone", 
              help="path to input files", metavar="character")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$Input_path)){
  print_help(opt_parser)
  stop("Input argument(s) missing.n", call=FALSE)
}

RootPath <- file.path(opt$Input_path) # Change to the path, where the result from DADA2 pipeline or similar pipelines is stored.
FigsPath <- file.path(opt$Input_path) # Path where output will be stored.


seqtab_file = (list.files(path = RootPath, pattern = "seqtab_nochim.rds" ,full.names=TRUE, recursive = F))


seqtab = readRDS(file=seqtab_file) # Input is the ASV table.

blast_file = (list.files(path = RootPath, pattern ="output_blast_results$",full.names=TRUE, recursive = F))


blast_results <- read.table(blast_file, fill=TRUE, sep ="\t") # Read Blast result from above.
colnames(blast_results) <- c( "QueryID",  "SubjectID", "Perc.Ident",
"Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
"S.start", "S.end", "E", "Bits" ) #Order Blast results.


blast_results = do.call( rbind,
        lapply( split(blast_results, blast_results[,c("QueryID") ]), function(d){
                                         d[which.max(d$Bits), ] } )
        )

summary(blast_results)
dim(blast_results)


lca_file = (list.files(path = RootPath, pattern ="lca.03_98.tabular$",full.names=TRUE, recursive = F))

lca_table = read.table(lca_file, sep="\t",fill=TRUE) # Read LCA output file
colnames(lca_table) = c("Query","LCA_rank","LCA_taxon","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Method") # Combines both LCA and top blast assignment. 
tax = as.matrix(lca_table) #makes as matrix
rownames(tax) <- tax[,1] # makes first col as rownames
tax = tax[,-1] #removes col which became rownames
tax = as.matrix(tax) # makes as matrix



################## parse seqtab ####################

seqtab_clean = seqtab[, colnames(seqtab) %in% blast_results$QueryID]

seqtab_names =sub('.*\\/', '', rownames(seqtab_clean))
seqtab_names = sub('\\..*', '', (seqtab_names))


rownames(seqtab_clean) = seqtab_names


################# sequence analysis evaluation ######
#track <- read.table(file.path(RootPath, "track_sequence_stats.csv"),
#                                   sep = ",", header = TRUE)

#summary(track)

#track$sample_id = str_extract(track$X, "[^_]+")

#seqtab_samplesum = rowSums(seqtab_clean)


################# Sequence data analysis and visualization ###############
dim(seqtab_clean)
am_otu_n <- otu_table(t(seqtab_clean), taxa_are_rows = TRUE)
dim(am_otu_n)
am_tax <- tax_table(as.matrix(tax))
am_physeq_clean_n <- phyloseq(am_otu_n, am_tax) # Combines ASV/OTU and Taxonomic table.

saveRDS(am_physeq_clean_n, file.path(FigsPath, file="phyloseq_mifish_lca.rds")) ## Phyloseq object to use for further analysis on local computer.


#' @title Summarize taxon composition
#'
#' @description This function takes a phyloseq object and returns phyloseq data (OTUs are now taxonomic ranks)
#' @param phylo_seq_object A phyloseq object with an OTU table and phylogenetic tree.
#' @param taxonomic_rank The taxonomic rank at which the data should be summarized
#' @details
#' Nice for making taxon summaries
#' @return A data.frame with the taxon name the mean, standard deviation as well as min and max values per samples.
#' @keywords phyloseq
#' @export
#' @examples
#' summarize_taxa()

summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    plyr::alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    tapply(counts, taxonomy, sum)
  }
}


#' @title Summarize taxon composition
#'
#' @description This function takes a phyloseq object and returns phyloseq data (OTUs are now taxonomic ranks)
#' @param phylo_seq_object A phyloseq object with an OTU table and phylogenetic tree.
#' @param taxonomic_rank The taxonomic rank at which the data should be summarized
#' @details
#' Nice for making taxon summaries
#' @return A data.frame with the taxon name the mean, standard deviation as well as min and max values per samples.
#' @keywords phyloseq
#' @export
#' @examples
#' phyloseq_summarize_taxa()

phyloseq_summarize_taxa <- function(phylo_seq_object, taxonomic_rank = rank_names(phylo_seq_object)[1], errorIfNULL = TRUE) {
  taxa <- as(phyloseq::tax_table(phylo_seq_object, errorIfNULL)[, taxonomic_rank], 'character')
  sum_tax_table <- summarize_taxa(as(phyloseq::otu_table(phylo_seq_object), 'matrix'), taxa)
  phyloseq::phyloseq(phyloseq::otu_table(sum_tax_table, taxa_are_rows = TRUE),
  	phyloseq::sample_data(phylo_seq_object, FALSE))
}

species_table_lca_n <- otu_table(phyloseq_summarize_taxa(am_physeq_clean_n, 'LCA_taxon'))
#species_table_n[species_table_n < 11] = 0
#species_table_n[species_table_n > 10] = 1

write.csv(t(species_table_lca_n), "species_table_lca_n.csv"))


print("Use the" $phyloseq_name " for downstream")
