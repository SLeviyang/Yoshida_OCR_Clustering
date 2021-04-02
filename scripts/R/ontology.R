library(gprofiler2)
library(dplyr)

GO_enrichment <- function(genes)
{
  
  #if (species=="mouse")
  #  org <- "mmusculus"
  #else
  #  org <- "hsapiens"
  d <- gene_enrichment.GO(genes, organism="mmusculus", significant = F) %>%
    dplyr::select(p_value, intersection_size, term_size, query_size,
                  recall, term_name) 
  return (d)
}

gene_enrichment.GO <- function(genes, organism, significant=T)
{
  d <- gprofiler2::gost(genes, organism=organism, correction_method = "fdr",
                        sources="GO:BP", significant=significant)
  go_d <- dplyr::filter(d$result, grepl("GO:", d$result$term_id)) %>% dplyr::arrange(p_value)
  
  return (go_d)
}

args = commandArgs(trailingOnly=TRUE)
print(args)
 
if (length(args) != 2) {
  print("arguments should be gene_file output_file")
  print(args)
  stop("")
}
gene_file = args[1]
output_file = args[2]

genes <- read.csv(gene_file, header=F) %>% as.matrix %>% as.character()
d <- GO_enrichment(genes)
write.csv(d, output_file, row.names = F)
