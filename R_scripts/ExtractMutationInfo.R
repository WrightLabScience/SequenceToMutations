# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')
lf <- list.files('mutations_subtracted_ancestor/')

# search patterns to extract mutation information
pattern_mutation_category <- '^mutation_category=(.+)$'
pattern_gene_position <- '^gene_position=(.+)' # general regex for all gene positions. if nonsynomyous
# this next pattern is a little wild. Intergenic mutations give 2 gene positions: 
# 1 relative to the gene on the left and another relative to the gene on the right
# depending on the orientation of the flanking genes, both numbers could be + or -
# but if the genes have the same orientation, one number will + and one will be -, depending on which way the genes face
# WHAT COMPLICATES THIS is when one of the flanking "genes" is not a gene...
# ...but rather the pesky symbol "–" which is not "-" (minus sign) !!!!
# the gene position relative to one of these is just the same symbol: "–" instead of "+56" or "-110", for 2 examples
pattern_gene_position_intergenic <- '^intergenic \\((((–|\\+|-)[0-9]*)/((–|\\+|-)[0-9]*))\\)$'
pattern_gene_position_other <- '^(coding|pseudogene|noncoding) \\(([0-9]+(-[0-9]+)?)/([0-9]+) nt\\)$' # \\1 gives the position (single # if SNP, number range if not SNP), \\3 gives the length of the gene
pattern_position_start <- 'position_start=([0-9]+)'
pattern_position_end <- 'position_end=([0-9]+)'
pattern_contig_info <- 'NODE_([0-9]+)_length_([0-9]+)_cov_.+' # \\1 is contig number, \\2 is contig length
pattern_gene_product <- '^gene_product=(.+)'
pattern_aa_ref_seq <- 'aa_ref_seq=(.+)'
pattern_aa_position <- 'aa_position=(.+)'
pattern_aa_new_seq <- 'aa_new_seq=(.+)'
pattern_nuc_ref_seq <- '^ref_seq=([ATCG])'
pattern_pgaptmp <- 'locus_tag=(.*\\[?pgaptmp_[0-9]+\\]?.*)'
pattern_pgaptmp_multiple <- 'locus_tags_inactivated=(.+)'
pattern_pgaptmp_overlapping <- 'locus_tags_overlapping=(.+)'


mutations_list <- vector('list', length = length(lf))
names(mutations_list) <- gsub('_mutations\\.txt', '', lf)

for (i in seq_along(lf)) { # loop over each file (table of mutations)
   muts <- read.table(file = paste0('mutations_subtracted_ancestor/', lf[i]),
                      header = FALSE,
                      fill = TRUE,
                      col.names = 1:63,
                      sep = '\t',
                      quote = '\"')
   names(muts) <- NULL
   
   df <- data.frame(matrix(ncol = 8,
                           nrow = nrow(muts),
                           dimnames = list(NULL, 
                                           c('Gene', 'PGAPtmp', 'Category', 'Location', 'Mutation', 'GenePosition', 'Contig', 'ContigLength'))))
   
   # mutation can be: 
   #     Mutation                Category Location
   #     3534A>G              if synon in coding
   #     (157/1008)ins1       if small_indel in coding
   #     (81-181/762)del101   if large_deletion in coding      GenePosition would be 81-181
   #     (-89/+291)ins1       if small_indel in intergenic     GenePosition would be -89/+291
   #     delWHOLE_CONTIG      if large_deletion in wholeGene
   
   for (j in seq_len(nrow(muts))) { # loop over each row (mutation)
      df$Contig[j]       <- as.integer(gsub(pattern = pattern_contig_info, 
                                            replacement = '\\1', 
                                            x = muts[j, 4]))
      df$ContigLength[j] <- as.integer(gsub(pattern = pattern_contig_info, 
                                            replacement = '\\2',
                                            x = muts[j, 4]))
      
      df$Gene[j] <- gsub(pattern = pattern_gene_product, 
                         replacement = '\\1', 
                         x = grep(pattern_gene_product, muts[j,], value=TRUE))
      
      df$PGAPtmp[j] <- gsub(pattern = pattern_pgaptmp, 
                            replacement = '\\1', 
                            x = grep(pattern_pgaptmp, muts[j,], value=TRUE))
      
      df$Category[j] <- gsub(pattern = pattern_mutation_category, 
                             replacement = '\\1',
                             x = grep(pattern_mutation_category, muts[j,], value=TRUE))
      
      
      gene_position <- grep(pattern_gene_position, muts[j,], value=TRUE)
      gene_position <- gsub('^gene_position=', '', gene_position)
      
      if (grepl('|', df$Category[j], fixed=TRUE)) { # for overlapping genes, a single mutation will be in both genes, shit is a mess
         df$GenePosition[j] <- gene_position
         cat('WOW!!', df$Category[j], '\n')
         next
      }
      
      # if gene_position is empty, this may be a whole contig deletion
      if (length(gene_position) == 0L) {
         pos_start <- as.integer(gsub(pattern_position_start, '\\1', grep(pattern_position_start, muts[j,], value=TRUE)))
         pos_end <- as.integer(gsub(pattern_position_end, '\\1', grep(pattern_position_end, muts[j,], value=TRUE)))
         df$GenePosition[j] <- paste0(pos_start, '-', pos_end)
         pgaptmps <- gsub(pattern_pgaptmp_multiple, '\\1', grep(pattern_pgaptmp_multiple, muts[j,], value=TRUE))
         overlap_pgaptmps <- grep(pattern_pgaptmp_overlapping, muts[j,], value=TRUE)
         if (length(overlap_pgaptmps) > 0L)
            pgaptmps <- paste0(pgaptmps, ',',
                               gsub(pattern_pgaptmp_overlapping, '\\1', overlap_pgaptmps))
         df$PGAPtmp[j] <- pgaptmps
         df$Location[j] <- 'wholeGene'
         if (pos_start == 1L & pos_end == df$ContigLength[j]) {
            df$Mutation[j] <- 'WHOLE_CONTIG_ABSENT'
         } else {
            df$Mutation[j] <- df$GenePosition[j]
         }
         next
      }
      
      # single position, could be insertion/deletion?
      if (grepl('^[0-9]+$', gene_position)) {
         df$GenePosition[j] <- gene_position
         # df$Location[j] <- 'coding' # not totally positive I can do this yet
         
         # if nonsense or nonsynonymous, write mutation in OLD_AA CODON NEW_AA
         if (df$Category[j] %in% c('snp_nonsense', 'snp_nonsynonymous')) {
            df$Location[j] <- 'coding'
            AA_ref <- gsub(pattern_aa_ref_seq, '\\1', grep(pattern_aa_ref_seq, muts[j,], value=TRUE))
            AA_codon_pos <- gsub(pattern_aa_position, '\\1', grep(pattern_aa_position, muts[j,], value=TRUE))
            AA_new <- gsub(pattern_aa_new_seq, '\\1', grep(pattern_aa_new_seq, muts[j,], value=TRUE))
            df$Mutation[j] <- paste0(AA_ref, AA_codon_pos, AA_new)
            
         # if synonymous
         } else if (df$Category[j] == 'snp_synonymous') {
            df$Location[j] <- 'coding'
            NUC_ref <- gsub(pattern_nuc_ref_seq, '\\1', grep(pattern_nuc_ref_seq, muts[j,], value=TRUE))
            NUC_new <- unlist(muts[j,6]) # column 6 in this case has the mutated nucleotide
            df$Mutation[j] <- paste0(df$GenePosition[j], NUC_ref, '>', NUC_new)
            
         # if a type we haven't handled yet
         } else {
            cat('STOP!!!!', i, j, df$Category[j])
            break
         }
         
         next # we can skip to the next mutation, because we got all the data we needed
      }
      
      # intergenic position
      if (grepl(pattern_gene_position_intergenic, gene_position)) {
         df$Location[j] <- 'intergenic'
         df$GenePosition[j] <- gsub(pattern_gene_position_intergenic, '\\1', gene_position)
         gene_position_for_mutation <- unlist(muts[j,5])
      }
      
      # coding position
      if (grepl(pattern_gene_position_other, gene_position)) {
         df$Location[j] <- 'coding'
         if (df$Category[j] == 'snp_noncoding') df$Location[j] <- 'noncoding'
         df$GenePosition[j] <- gsub(pattern_gene_position_other, '\\2', gene_position)
         gene_position_for_mutation <- paste0(df$GenePosition[j], '/', 
                                              gsub(pattern_gene_position_other, '\\4', gene_position))
      }
      
      # create mutation notations intergenic/coding indels AND snp_intergenics
      if (df$Category[j] %in% c('large_deletion', 'small_indel', 'large_insertion', 'large_amplification', 'large_substitution')) {
         indel_flag <- tolower(muts[j,1])
         indel_length_or_nucleotides <- unlist(muts[j,6])
         if (df$Category[j] == 'large_insertion') { # if this is a large insertion, we don't want to display the literal string of nucleotides that were inserted, but rather the length of the insertion
            indel_length_or_nucleotides <- length(unlist(strsplit(indel_length_or_nucleotides, split='')))
         }
         df$Mutation[j] <- paste0('(', gene_position_for_mutation, ')', indel_flag, indel_length_or_nucleotides)
      } else if (df$Category[j] %in% c('snp_intergenic', 'snp_pseudogene', 'snp_noncoding')) {
         NUC_ref <- gsub(pattern_nuc_ref_seq, '\\1', grep(pattern_nuc_ref_seq, muts[j,], value=TRUE))
         NUC_new <- unlist(muts[j,6]) # column 6 in this case has the mutated nucleotide
         df$Mutation[j] <- paste0(unlist(muts[j,5]), NUC_ref, '>', NUC_new) # column 5 has the contig position of the SNP if it is intergenic
      } else {
         cat('STOP!!!', i, j, df$Category[j], '\n')
         break
      }
   }
   
   mutations_list[[i]] <- df
}

save(mutations_list, file = 'mutations_list_Rdata/mutations_list.Rdata')









