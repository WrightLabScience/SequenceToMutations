### trimming raw sequencing reads based on quality

# set to path where SequenceToMutations diretory lives on your machine
setwd('~/Desktop/SequenceToMutations/')

# if not already installed on your machine
if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("DECIPHER")

library(DECIPHER)

# you need to have a folder with all of the raw fastq.gz files in the following directory:
# change this directory to the location on your computer where the files are
path_raw_reads <- 'reads/raw_reads/'
lf <- list.files(path_raw_reads, pattern = "fastq") # a list of all of the raw fastq files

# this will create another folder in the same directory as your raw reads directory that will hold the trimmed reads
path_trimmed_reads <- 'reads/trimmed_reads/'
dir.create(path_trimmed_reads)


# adapter sequences
i7 <- DNAStringSet("CAAGCAGAAGACGGCATACGAGAT")
i7rc <- reverseComplement(i7)
i5 <- DNAStringSet("AATGATACGGCGACCACCGAGATCTACAC")
i5rc <- reverseComplement(i5)

# loop over each raw .fastq.gz file
# this script may take a while to run. ~1 minute per .fastq.gz file, depending on their size - your mileage may vary
for (i in seq_along(lf)) {
   print(i)
   
   dna <- readQualityScaledDNAStringSet(paste0(path_raw_reads, lf[i]))  #default quality score is phred
   
   # Trimming DNA based on adapter sequences and on quality
   x <- TrimDNA(DNAStringSet(dna),
                left = c(i7, i5, i7rc, i5rc),
                right = c(i7, i5, i7rc, i5rc),
                quality = dna@quality)
   trimmed <- subseq(dna, start(x), end(x))
   trimmed <- trimmed[width(trimmed) > 0]
   
   # ask for regex help if this line doesn't make sense
   # you want your sample name to be unique for each sample and to correspond to something meaningful about your experiment
   # for example, {organism_name}_{strain_identifier}_{replicate_number}_{etc.}
   if (grepl('^VRSA', lf[i])) { # in my case, this indicated the reads corresponded to an ancestor strain, so a different file name format
      sample_name <- gsub(pattern = '^VRSA([0-9]+).+',
                          replacement = 'ancestor_\\1',
                          x = lf[i])
   } else { # otherwise, the reads corresponded to an evolved lineage = different file name format
      sample_name <- gsub(pattern = '([A-z0-9]+_(BHI|TSB)).+',
                          replacement = 'evolved_\\1',
                          x = lf[i])
   }
   
   writeQualityScaledXStringSet(trimmed, 
                                file = paste0(path_trimmed_reads, sample_name, '_trimmed.fastq.gz'), 
                                compress = TRUE)
}


# next step is assembly, which is run on the grid, so we need a job_map_file
# let's make txt file here
# we only assemble the ancestor
# in my case, with the data I am providing, there are 2 ancestor strains and 4 evolved lineages each

ancestor_trimmed_reads_file_names <- list.files(path_trimmed_reads, pattern='^ancestor')

write.table(ancestor_trimmed_reads_file_names,
            file = 'OSG_job_files/assemble_job_map.txt',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)








