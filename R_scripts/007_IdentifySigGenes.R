# Load necessary data.
setwd('~/Desktop/SequenceToMutations/')
load(file = 'RdataFiles/ancestor_gbk_dataframes.Rdata') # all ancestor assmebly gene and gene product information
load(file = 'RdataFiles/mutations_list_VRSAc50.Rdata') # all evolved lineage mutations tables
mutations_list <- mutations_list[1:2] # Just take VRSA-3 and -4

# Get a vector of gene products that are present in all strains.
gene_prods <- Reduce(intersect, sapply(gbkDF_list, '[[', 'GeneProduct'))

groups <- c('VAN-unexposed', 'VAN-exposed')
strains <- gsub('^Annot_(.+)', '\\1', names(gbkDF_list))

# Create a matrix of counts. Rows are each gene product. Columns are each group.
sig_genes <- matrix(0L,
                    nrow = length(gene_prods),
                    ncol = length(groups),
                    dimnames = list(gene_prods,
                                    groups))

# get list of gene product names that have any mutation in a given lineage.
getGenes <- function(lin, s) {
   pgap_ids <- unique(unlist(lin$PGAPtmp))
   pgap_ids <- unlist(strsplit(pgap_ids, ',')) # 382
   pgap_ids <- unlist(strsplit(pgap_ids, '/')) # 382
   pgap_ids <- unique(pgap_ids) # 217
   pgap_ids <- pgap_ids[pgap_ids != '–']
   if (!all(pgap_ids %in% gbkDF_list[[s]]$PGAPtmp)) {
      print('wtf!')
      break
   }
   
   genes <- gbkDF_list[[s]]$GeneProduct[which(gbkDF_list[[s]]$PGAPtmp %in% pgap_ids)]
   genes <- genes[genes != '–']
   genes <- genes[genes != 'hypothetical protein']
   genes <- genes[genes %in% gene_prods]
   return(genes)
}


for (s in seq_along(strains)) { # Loop over strains
   for (g in seq_along(groups)) { # Loop over groups (VAN+ and VAN-)
      for (l in seq_along(mutations_list[[s]][[g]])) {
         lineage_mut_list <- mutations_list[[s]][[g]][[l]]
         if (is.null(lineage_mut_list))
            next
         
         genes <- getGenes(lineage_mut_list, s)
         sig_genes[genes, groups[g]] <- sig_genes[genes, groups[g]] + 1
      }
   }
}

sig_genes <- data.frame(sig_genes)

num_lins <- 48 # Mine was actually 96 VAN-unexposed/exposed lineages each in the full experiement. Here, I am only using 2 (out of 4) strains, so only 48 lineages.
pvals <- oddsRs <- numeric(nrow(sig_genes))

for (i in seq_len(nrow(sig_genes))) {
   # mutation counts
   ve_mut <- sig_genes$VAN.exposed[i]
   ve_nomut <- num_lins - ve_mut
   vu_mut <- sig_genes$VAN.unexposed[i]
   vu_nomut <- num_lins - vu_mut
   
   mat <- matrix(c(ve_mut, ve_nomut,
                   vu_mut, vu_nomut),
                 ncol=2)
   
   ft <- fisher.test(mat)
   pvals[i] <- ft$p.value
   oddsRs[i] <- ((ve_mut + 1) * (vu_nomut + 1)) / ((ve_nomut + 1) * (vu_mut + 1))
}
rm(i)

# Add columns to data.frame for transformed p-values and odds ratios*
sig_genes$pval <- pvals
sig_genes$OR <- oddsRs


# Show top 6 genes enriched in VAN-exposed group.
head(sig_genes[order(sig_genes$pval), ])

# Bonferroni-corrected p-value=0.01 threshold.
pval_cutoff <- 0.01 / nrow(sig_genes)

# Plotting variables.
col_vec <- c('ns' = '#bbbbbb', 'sig' = '#000000')
pch_vec <- 16
pcex <- 0.75
tcex <- 1

# Function to plot proportion of 
mutsXY <- function(df) {
   plot(NA,
        xlim = c(0, 1),
        ylim = c(0, 1),
        xaxt = 'n', xaxs = 'i', xlab = '',
        yaxt = 'n', yaxs = 'i', ylab = '')
   abline(a=0, b=1, lty=2)
   
   points(x = df$VAN.exposed / num_lins,
          y = df$VAN.unexposed / num_lins,
          xpd = NA,
          pch = pch_vec,
          col = ifelse(df$pval < pval_cutoff, col_vec['sig'], col_vec['ns']),
          bg = ifelse(df$pval < pval_cutoff, col_vec['sig'], col_vec['ns']),
          cex = pcex)
   
   title(xlab = 'proportion in VAN-exposed', line=1.75)
   title(ylab = 'proportion in VAN-unexposed', line=2.25)
   
   axis(side = 1,
        at = seq(0, 1, 0.2),
        labels = seq(0, 1, .2))
   axis(side = 2,
        las = 1,
        at = seq(0, 1, 0.2),
        labels = seq(0, 1, .2))
}

volcano <- function(df, xlimit=6, ylim=c(-0.25,30)) {
   plot(NA, 
        xlim = c(-xlimit, xlimit),
        ylim = ylim,
        yaxs = 'i',
        xlab = '',
        ylab = '',
        yaxt = 'n')
   abline(h = -log10(pval_cutoff), lty = 2)
   abline(v = 0, lty=2)
   
   points(x = log(df$OR),
          y = -log10(df$pval),
          pch = pch_vec,
          cex = pcex,
          col = ifelse(df$pval < pval_cutoff, col_vec['sig'], col_vec['ns']),
          bg = ifelse(df$pval < pval_cutoff, col_vec['sig'], col_vec['ns']))
   
   title(xlab = 'Vancomycin specificity score', line=1.75, xpd=NA)
   title(ylab = expression('-Log'[10]*'(p-value)'), line=1.75)
   
   axis(side = 2, at = seq(0, 25, 5), las = 1)
   
   text(x = xlimit * 1.06,
        y = -log10(pval_cutoff)*1.01 + c(0.5, -0.5),
        adj = c(1, 0.5),
        cex = tcex,
        labels = c('Bonferonni-corrected', 'p-value cutoff'))
}


# Make plots
{
   pdf(file = 'Plots/significant_mutations.pdf', width=10.5, height=5.25)
   par(mfrow=c(1,2), tck=-0.015, mgp=c(1.5,0.5,0), mar=c(3,3,2,1))
   volcano(df=sig_genes)
   legend('bottomright', col=col_vec, legend=names(col_vec), pch=16)
   mutsXY(df=sig_genes)
   dev.off()
}



















