# make the job_map.txt file to run breseq on the grid

# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')

# get the file names of the trimmed reads for evolved lineages AND ANCESTOR STRAINS
lf_evo <- list.files('reads/trimmed_reads/')
names(lf_evo) <- gsub('(ancestor|evolved)_([0-9]+).+', '\\2', lf_evo)

# get the file names of the annotated and assembled ancestor strains
lf_gbk <- list.files('pgap/')
# get the strain name (number in this case) from the annotated genomes
names(lf_gbk) <- gsub('Annot_([0-9]+)\\.gbk', '\\1', lf_gbk)

# behold, the simple beauty of named indexing
df <- data.frame(lf_evo, lf_gbk[names(lf_evo)])

write.table(x = df,
            file = 'OSG_job_files/map_job_map.txt',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
