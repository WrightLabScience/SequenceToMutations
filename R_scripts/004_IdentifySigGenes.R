# Load necessary data.
load(file = 'RdataFiles/ancestor_gbk_dataframes.Rdata') # all ancestor assmebly gene and gene product information
load(file = 'RdataFiles/mutations_list.Rdata') # all evolved lineage mutations tables

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

num_lins <- 2 # Mine was actually 96 VAN-unexposed/exposed lineages each
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
sig_genes$pval <- -log10(pvals)
sig_genes$OR <- oddsRs


# Show top 6 genes enriched in VAN-exposed group.
head(sig_genes[order(log(sig_genes$OR), decreasing=T), ])

# Bonferroni-corrected p-value=0.01 threshold.
pval_cutoff <- -log10(0.01 / nrow(sig_genes))

# Plotting variables.
col_vec <- '#bbbbbb'
pch_vec <- 16
pcex <- 0.75
tcex <- 1

# Function to plot proportion of 
mutsXY <- function(df) {
   plot(x = df$VAN.exposed / num_lins,
        y = df$VAN.unexposed / num_lins,
        xlim = c(0, 1),
        ylim = c(0, 1),
        xpd = NA,
        xaxt = 'n', xaxs = 'i', xlab = '',
        yaxt = 'n', yaxs = 'i', ylab = '',
        pch = pch_vec,
        col = col_vec,
        bg = col_vec,
        cex = pcex)
   title(xlab = 'proportion in VAN-exposed', line=1.75)
   title(ylab = 'proportion in VAN-unexposed', line=1.75)
   axis(side = 1,
        at = seq(0, 1, 0.2),
        labels = seq(0, 1, .2))
   axis(side = 2,
        las = 1,
        at = seq(0, 1, 0.2),
        labels = seq(0, 1, .2))
   abline(a=0, b=1)
}

addLegend <- function() {
   x <- -0.63
   y <- 1.01
   xshift <- rep(0.034, 2)
   yshift <- c(0.05, 0.094)
   text(x = x + xshift + 0.025,
        y = y - yshift,
        labels = c('         operon genes', 
                   'other genes'),
        adj = c(0, 0.5),
        cex = tcex,
        xpd = NA)
   
   points(x = x + xshift, 
          y = y - yshift - 0.0025,
          pch = c(17, 19),
          cex = pcex+0.1,
          col = c('black', graycol),
          xpd = NA)
   
   seg_y_pos <- 0.874
   segments(x0 = c(x, x),
            x1 = c(-0.202, x),
            y0 = c(seg_y_pos, seg_y_pos),
            y1 = c(seg_y_pos, 1),
            xpd = NA)
   
   text(x = x + 0.102,
        y = y - yshift[1],
        cex = tcex,
        labels = 'vanA',
        xpd = NA,
        font = 3)
}

volcano <- function(df, xlimit=6, ylim=c(-0.25,30)) {
   plot(x = log(df$OR),
        y = df$pval,
        xlim = c(-xlimit, xlimit),
        ylim = ylim,
        yaxs = 'i',
        xlab = '',
        ylab = '',
        pch = pch_vec,
        cex = pcex,
        bg = col_vec,
        yaxt = 'n',
        col = col_vec)
   
   title(xlab = 'Vancomycin specificity score', line=1.75, xpd=NA)
   title(ylab = expression('-Log'[10]*'(p-value)'), line=1.75)
   
   abline(h = pval_cutoff, lty = 2, lwd=0.8)
   abline(v = 0, lty=3, lwd=1.2, col='#00000066')
   
   axis(side = 2, at = seq(0, 25, 5), las = 1)
   
   text(x = 6.35,
        y = 5.4 + c(0.5, -0.5),
        adj = c(1, 0.5),
        cex = tcex,
        labels = c('Bonferonni-corrected', 'p-value cutoff'))
}

# Make plots
{
   pdf(file = 'Plots/significant_mutations.pdf', width=10, height=5)
   par(mfrow=c(1,2), tck=-0.015, mgp=c(1.5,0.5,0), mar=c(3,3,2,1))
   volcano(df=sig_genes)
   mutsXY(df=sig_genes)
   dev.off()
}



















