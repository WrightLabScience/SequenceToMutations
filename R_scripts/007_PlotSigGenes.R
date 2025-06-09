setwd('~/Desktop/SequenceToMutations/')
load(file='RdataFiles//ddl_info.Rdata') # info on the start and end positions of ddl and surrounding genes for each strain. Will need to get this info for whatever genes you want to plot
load(file='RdataFiles/mutations_list_VRSAc50.Rdata')
mutations_list <- mutations_list[1:2] # Just take VRSA-3 and -4
names(mutations_list[[1]])[4] <- 'TSB-'
names(mutations_list[[2]])[4] <- 'TSB-'

strains <- c(3,4) # Just VRSA-3 and -4 here.
media <- c('BHI','TSB')
groups <- c('+', '-')

numer <- denom <- array(0L, dim=c(length(strains), length(groups), length(media)), dimnames=list(paste0('VRSA-', strains), groups, media))
for (s in seq_along(strains)) {
   for (m in media) {
      for (j in groups) {
         x <- sapply(mutations_list[[s]][[paste0(m, j)]], \(x) any(grepl('D-alanine--D-alanine ligase', x$Gene)))
         numer[s, j, m] <- sum(x)
         denom[s, j, m] <- length(x)
      }
   }
}; rm(s, m, j, x)


cols <- c('UDP-ligase' = "#beaed4", 
          'ligase' = "#386cb0", 
          'rod' = "#03C04a")
del_col <- "#797EF6"
ins_col <- "#FA9C1B"
cex <- 0.84
leg_cex <- 1.1
cex1 <- 0.75
cex2 <- 1.1



pdf(file = 'Plots/DDL_region.pdf',
    bg = "white",
    height = 7.5,
    width = 6.25)
layout(mat = matrix(c(rep(1:8, each=3),
                      rep(0, 6),
                      rep(9:14, each=2)), 
                    ncol=6, 
                    byrow=T),
       heights = c(1,1,1,1,0.8, 1.4, 1.4))
par(oma = c(0, 0, 1, 0))


# Loop over each
for (i in seq_along(strains)) {
   # get position and gene information
   {
      this <- ddl_info[[i]][, c('start position', 'end position')]
      if (i == 4) { # VRSA-10 had a reversed synteny in the ddl region. For consistency across strains, I flipped it here and elsewhere in this script.
         this <- max(this) - this[3:1, ] + 100
      } else {
         this <- this - min(this) + 100
      }
      
      start <-0
      end <- max(this) + 100
      
      colGene <- character(nrow(this))
      products <- character(nrow(this))
      
      for (j in seq_len(nrow(this))) {
         if (rownames(this)[j] == "D-alanine--D-alanine ligase") {
            colGene[j] <- cols["ligase"]
            products[j] <- "ddl"
            
         } else if (grepl("UDP", rownames(this)[j])) {
            colGene[j] <- cols["UDP-ligase"]
            products[j] <- "murF"
            
         } else if (grepl("RodA", rownames(this)[j])) {
            colGene[j] <- cols["rod"]
            products[j] <- "mrdB"
         }
      }
   }
   
   for (s in c(1, 3)) {
      # setup plotting area
      {
         if (s == 1) par(mar = c(0.2, 0.1, 5.6, 2))
         if (s == 3) par(mar = c(0.2, 2,   5.6, 0.1))
         
         xlim <- c(start - 10, end + 10)
         plot(1e9,
              xlim = xlim,
              ylim = c(-1, 1),
              axes = FALSE,
              ann = FALSE)
         
         segments(x0 = start - 25,
                  x1 = end + 25,
                  y0 = 0,
                  lwd = 3)
         
         shift <- c(240, 240, -240)
         rect(xleft = this[, "start position"],
              xright = this[, "end position"] + shift,
              ybottom = -0.7,
              ytop = 0.7,
              border = NA,
              col = colGene)
         
         shift <- c(250, 250, -250)
         for (p in 1:3) {
            polygon(x = c(this[p, "end position"] + shift[p],
                          this[p, "end position"] + shift[p],
                          this[p, "end position"]),
                    y = c(-1.1, 1.1, 0),
                    border = NA,
                    col = colGene[p])
         }
         
         yb <- 1
         yt <- 9
         
         text(x = apply(this, 1, function(x) mean(c(x[1], x[2]))) + shift/4,
              y = 0,
              xpd = NA,
              cex = 1,
              col = "white",
              font = 3,
              labels = products)
         
         
         x0 <- ifelse(s == 1, max(xlim) - 95, min(xlim) + 95)
         
         segments(x0 = x0 - 105,
                  x1 = x0 + 105,
                  y0 = c(mean((24.9 + 1:8) / 4),
                         mean((5 + 1:16) / 4)),
                  xpd = NA)
         
         col_idx <- ifelse(s == 1, 1, 2)
         mtext(text = c(numer[i, col_idx, 2], denom[i, col_idx, 2],
                        numer[i, col_idx, 1], denom[i, col_idx, 1]),
               side = ifelse(s == 1, 4, 2),
               line = -1.2,
               las = 1,
               xpd = NA,
               adj = 0.5,
               cex = cex1,
               at = c(8.075, 6.675, 4.075, 2.675))
         
         
         if (s == 3) {
            mtext(side = 2,
                  las = 1,
                  cex = cex1,
                  line = 1.95,
                  adj = 0.5,
                  text = paste0('VRSA-', strains[i]),
                  font = 2,
                  at = 0)
            mtext(text = c("BHI", "TSB"),
                  side = 2,
                  las = 1,
                  cex = cex1,
                  line = 1.9,
                  adj = 0.5,
                  at = c(mean(c(5, 22))/4, mean(c(34, 25))/4))
         }
         
         axline <- -0.2
         axis(side = 2,
              at = c(5, 22) / 4,
              tck = 0.25,
              line = axline,
              labels = c("", ""),
              xpd = NA)
         axis(side = 2,
              at = c(25, 34) / 4,
              tck = 0.25,
              line = axline,
              labels = c("", ""),
              xpd = NA)
         axis(side = 4,
              at = c(5, 22) / 4,
              tck = 0.25,
              line = axline,
              labels = c("", ""),
              xpd = NA)
         axis(side = 4,
              at = c(25, 34) / 4,
              tck = 0.25,
              line = axline,
              labels = c("", ""),
              xpd = NA)
      }
      
      # per lineage mutations
      for (m in c(s, s+1)) {
         h <- ifelse(m == 1 || m == 3, 0, 20) + 5    
         
         for (l in seq_along(mutations_list[[i]][[m]])) {
            #segments(x0 = -5000, x1 = 5000, y0 = (h+l)/4, xpd = NA, lwd = 0.5)
            
            theseMuts <- mutations_list[[i]][[m]][[l]]
            
            genes <- unlist(theseMuts$Gene)
            w <- which(grepl(rownames(this)[1], genes) | grepl(rownames(this)[2], genes) | grepl(rownames(this)[3], genes))
            
            if (length(w) > 0) {
               for (k in w) {                             # each mutation that is not whole contig...
                  thisMut <- unlist(theseMuts[k, ])
                  
                  if (grepl("ins|del", thisMut["Mutation"])) {
                     pos1 <- integer(2L)
                     
                     if (thisMut["Location"] == "intergenic") {
                        pos1[1] <- as.numeric(thisMut["GenePosition"])
                        len <- as.numeric(gsub(".+(del|ins)([0-9]+)", "\\2", thisMut["Mutation"]))
                        pos1[2] <- pos1[1] + len
                        
                     } else if (grepl("del", thisMut["Mutation"])) { 
                        pos1 <- as.numeric(strsplit(thisMut["GenePosition"], "-")[[1]])
                        
                     } else if (grepl("ins", thisMut["Mutation"])) {
                        pos1[1] <- as.numeric(thisMut["GenePosition"])
                        len <- as.numeric(gsub(".+ins([0-9]+)", "\\1", thisMut["Mutation"]))
                        pos1[2] <- pos1[1] + len
                     }
                     
                     comp <- ifelse(grepl("rod|Rod", thisMut["Gene"]), 1, -1)
                     pos1[1] <- pos1[1] - 3
                     pos1[2] <- pos1[2] + 3
                     pos1 <- pos1 * comp
                     
                     line_col <- ifelse(grepl("del", thisMut["Mutation"]), del_col, ins_col)
                     segments(x0 = this[thisMut["Gene"], "start position"] + pos1[1],
                              x1 = this[thisMut["Gene"], "start position"] + pos1[2],
                              y0 = (h+l) / 4,
                              col = line_col,
                              lwd = 2.5,
                              xpd = NA)
                     
                  } else {
                     pch <- "w"
                     if (thisMut["Category"] == "snp_synonymous")    pch <- 1
                     if (thisMut["Category"] == "snp_nonsynonymous") pch <- 16
                     if (thisMut["Category"] == "snp_pseudogene")    pch <- 4
                     if (thisMut["Category"] == "snp_nonsense")      pch <- 6
                     if (thisMut["Category"] == "snp_intergenic") {   
                        pch <- 0
                        pos <- as.integer(strsplit(thisMut['GenePosition'], '/')[[1]][1])
                        pos <- this[strsplit(thisMut['Gene'], '/')[[1]][1], 1] - pos
                     } else {
                        pos1 <- as.numeric(thisMut["GenePosition"])
                        comp <- ifelse(grepl("rod|Rod", thisMut["Gene"]), 1, -1)
                        pos1 <- pos1 * comp
                        pos <- this[thisMut[1], 1] + pos1
                     }
                     
                     
                     if (pch == 1)
                        print(paste0(unlist(thisMut[1]), ", ", unlist(thisMut[6]), ", ", i, ", ", m, ", ", l, ", ", k))
                     
                     
                     
                     points(x = pos,
                            y = (h+l) / 4,
                            pch = pch,
                            cex = 1,
                            xpd = NA)
                  }
               }
            }
         }
      }
   }
}

dev.off()
