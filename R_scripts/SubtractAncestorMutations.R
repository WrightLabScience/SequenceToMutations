# from breseq output mutations list files
# remove mutations that occurred on contigs less than 500 base pairs ~ shorter contigs have more mutations = probably noise
# 
# ancestor trimmed reads mapped to ancestor genomes will reveal mutations...but why??
# these are likely the result of errors in sequencing, assembly, or something else
# if those mutations show up in the evolved, remove them


# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')


lf <- list.files('breseq_output/', pattern='out$')
strains <- c(3, 4)
MINIMUM_CONTIG_LENGTH <- 500

get_gd_list <- function(fname, mode='') {
   x <- read.delim(file = paste0('breseq_output/', fname, '/data/annotated.gd'),
                   header = FALSE, 
                   skip = 12, 
                   fill = TRUE, 
                   col.names = 1:63)
   
   if (mode == 'ancestor')  w <- which(x$X1 == "RA")[1] - 1
   if (mode == 'evolved') w <- which(x$X1 == "RA" | x$X1 == "MC" | x$X1 == "JC")[1] - 1
   x <- x[1:w, ]
   
   rem <- logical(nrow(x))
   for (i in seq_along(rem)) rem[i] <- any(x[i, ] == 'gene_product=–/–' | x[i, ] == '–')
   rem <- which(rem)
   if (length(rem) > 0) x <- x[-rem, ]
   
   rem <- which(as.integer(gsub('NODE_[0-9]{1,3}_length_([0-9]{1,8})_cov_.+', '\\1', x[,4])) < MINIMUM_CONTIG_LENGTH)
   if (length(rem) > 0) x <- x[-rem, ]
   
   return(x)
}

num_muts <- integer(0L)
all_lineages <- character(0L)

for (s in seq_along(strains)) {
   ancestor_muts <- get_gd_list(fname = paste0('ancestor_', strains[s], '_out'), mode='ancestor')
   
   evolved_fnames <- grep(paste0('^evolved_', strains[s]), lf, value=TRUE)
   lineages <- gsub('^evolved_([0-9]{1,2}[a-z]_[A-Z]{3})_out$', '\\1', evolved_fnames)
   all_lineages <- c(all_lineages, lineages)
   
   for (l in seq_along(evolved_fnames)) {
      cat('Strain:', s, '- evolved lineage:', lineages[l], '\n')
      
      evolved_muts <- get_gd_list(evolved_fnames[l], mode='evolved')
      
      check_cols <- c(1, 4:6, 8:11)
      
      # subtract ancestor --> write to file
      trow <- logical(0L)
      for (i in seq_len(nrow(ancestor_muts))) {
         x <- unlist(ancestor_muts[i, check_cols])
         w <- which(x[1] == evolved_muts$X1
                    & x[2] == evolved_muts$X4 
                    & x[3] == evolved_muts$X5 
                    & x[4] == evolved_muts$X6 
                    & x[5] == evolved_muts$X8 
                    & x[6] == evolved_muts$X9 
                    & x[7] == evolved_muts$X10 
                    & x[8] == evolved_muts$X11)
         trow <- c(trow, w)
      }
      if (length(trow) > 0) evolved_muts <- evolved_muts[-trow, ]
      
      num_muts <- c(num_muts, nrow(evolved_muts))
      
      # subtracted ancestor
      write.table(evolved_muts,
                  file = paste0('mutations_subtracted_ancestor/', lineages[l], '_mutations.txt'),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE,
                  sep = '\t')
   }
}

# distribution of number of mutations
hist(num_muts, breaks=100, ylim=c(0, 25))

df <- data.frame(num_muts, lineage=all_lineages)
df

num_muts_by_strain <- split(df$num_muts, gsub('^([0-9]+).+', '\\1', df$lineage))

# calculate on a strain by strain basis the average number of mutations per evolved lineage, etc.
for (x in seq_along(num_muts_by_strain)) {
   cat('strain:', names(num_muts_by_strain[x]), '\n')
   cat('mean:', mean(num_muts_by_strain[[x]]), '\n')
   cat('median:', median(num_muts_by_strain[[x]]), '\n')
   cat('min:', min(num_muts_by_strain[[x]]), '\n')
   cat('max:', max(num_muts_by_strain[[x]]), '\n')
   cat('\n')
}






