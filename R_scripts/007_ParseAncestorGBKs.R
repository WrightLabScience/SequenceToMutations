# parse ancestor gbk files to get list of gene product names and pgaptmp IDs

# each .gbk file is constructed withL
# HEADER describing the whole file
# EACH CONTIG has:
#     HEADER
#     LIST OF GENES AND THEIR INFORMATION (length, direction, names, etc.)
#     CONTIG sequence (lowercase ATCG's)


# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')

lf <- list.files('pgap/')
gbkDF_list <- vector('list', length=length(lf))
names(gbkDF_list) <- gsub('\\.gbk', '', lf)

for (l in seq_along(lf)) { # for each annotated and assembled genome (one per ancestor strain)
   gbk <- readLines(con = paste0('pgap/', lf[l]))
   
   MINIMUM_CONTIG_LENGTH <- 500L
   
   contig_start_positions <- grep('^FEATURES +Location/Qualifiers', gbk)
   contig_DNAseq_start_positions <- grep('^ORIGIN +', gbk)
   
   contig_lengths <- as.integer(gsub(pattern = ' +source + 1\\.\\.([0-9]+)', 
                                     replacement = '\\1', 
                                     x = gbk[contig_start_positions + 1]))
   w <- which(contig_lengths > MINIMUM_CONTIG_LENGTH)
   
   contig_start_positions <- contig_start_positions[w]
   contig_DNAseq_start_positions <- contig_DNAseq_start_positions[w]
   contig_lengths <- contig_lengths[w]
   
   df_list <- vector(mode = 'list', length = length(contig_start_positions))
   
   for (i in seq_along(contig_start_positions)) { # for each contig
      cat(i, '\n')
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
      
      product <- grep(' +/product=\"(.+)$', contig, value=TRUE)
      product <- gsub('/|\"', '', product)
      df$GeneProduct <- gsub(pattern = ' +product=(.+)', 
                             replacement = '\\1', 
                             x = product)
      
      pgaptmps <- grep('pgaptmp_[0-9]+', contig, value=TRUE)
      pgaptmps <- gsub('.*(pgaptmp_[0-9]+).*', '\\1', pgaptmps)
      pgaptmps <- unique(pgaptmps)
      if (length(pgaptmps) != nrow(df)) {
         cat('PGAP IDs WRONG !!!')
         break
      }
      df$PGAPtmp <- pgaptmps
      
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
      
      df_list[[i]] <- df
   }
  
   # add this large dataframe to the list of gbk data frames 
   gbkDF_list[[l]] <- do.call(rbind, df_list)
}


save(gbkDF_list, file = 'RdataFiles/ancestor_gbk_dataframes.Rdata')
















