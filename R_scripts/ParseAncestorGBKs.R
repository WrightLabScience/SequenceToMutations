# parse ancestor gbk files to get list of gene product names and pgaptmp IDs

# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')

lf <- list.files('pgap/')

gbk <- readLines(con = paste0('pgap/', lf[1]))

# each .gbk file is constructed withL
# HEADER describing the whole file
# EACH CONTIG has:
#     HEADER
#     LIST OF GENES AND THEIR INFORMATION (length, direction, names, etc.)
#     CONTIG sequence (lowercase ATCG's)



contig_start_positions <- grep('^FEATURES +Location/Qualifiers', gbk)
contig_DNAseq_start_positions <- grep('^ORIGIN +', gbk)

contig_lengths <- as.integer(gsub(pattern = ' +source + 1\\.\\.([0-9]+)', 
                                  replacement = '\\1', 
                                  x = gbk[contig_start_positions + 1]))

for (i in seq_along(contig_start_positions)) { # for each contig
   contig <- gbk[contig_start_positions[i]:contig_DNAseq_start_positions[i]]
   gene_positions <- grep(' +gene  +', contig)
   gene_locations <- gsub(' +gene +', '', contig[gene_positions])
   
   df <- data.frame(matrix(NA,
                           nrow = length(gene_positions),
                           ncol = 9,
                           dimnames = list(
                              NULL,
                              c('Contig', 'ContigLength', 'GeneProduct', 'PGAPtmp', 'complement', 'start_pos', 'end_pos', 'start_uncertainty', 'end_uncertainty')
                           )))
   df$Contig <- i
   df$ContigLength <- contig_lengths[i]
   
   heads <- mapply(gene_positions[-length(gene_positions)], 
                   gene_positions[-1] - 1, 
                   FUN = ':')
   for (g in seq_along(heads)) {
      gene <- contig[heads[[g]]]
      gsub()
      df$GeneProduct[g] <- gsub(pattern = ' +/product=\"(.+)\"$', 
                                replacement = '\\1', 
                                x = grep(' +/product=\"(.+)$', gene, value=TRUE))
   }
   
   # information about position, direction, and boundary uncertainty of each gene
   df$complement <- grepl('complement', gene_locations)
   gene_locations <- gsub('complement\\(|\\)', '', gene_locations)
   start_end_pos <- strsplit(gene_locations, split='..', fixed=TRUE)
   start_pos <- sapply(start_end_pos, '[', 1)
   end_pos <- sapply(start_end_pos, '[', 2)
   df$start_uncertainty <- grepl('<', start_pos)
   df$end_uncertainty <- grepl('>', end_pos)
   df$start_pos <- as.integer(gsub('<', '', start_pos))
   df$end_pos <- as.integer(gsub('>', '', end_pos))  
}





















