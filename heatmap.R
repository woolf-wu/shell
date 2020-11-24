# Step0 Before starting your project --------------------------------------

#installPackage("data.table")
#
#installPackage("BiocManager")
#
#installPackage("PoisonAlien/maftools")

rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()

Sys.setlocale("LC_ALL","English")
.libPaths( c( "E:/R-packages", "C:/Program Files/R/R-3.6.3/library") )


# Step1 Starting your project -------------------------------------------

library( "pheatmap" )
library("RColorBrewer")


sig_raw <- read.csv(file = "DMR_result.xls",
           header = T, sep = '\t',
           #row.names = 1,
           stringsAsFactors = FALSE,
           comment.char = "#")
mean(sig_raw$C)
mean(sig_raw$KS)
mean(sig_raw$KM)
mean(sig_raw$KH)
dim(sig_raw)
dim(sig_raw[sig_raw$C>0.5,])

for (sample in unique(sig_raw$Sample)) {
  raw <- sig_raw[sig_raw$Sample == sample,]
  anno_raw <- raw[,c(4,1)]
  raw <- raw[,c(4,2,3)]
  rownames(raw) = raw[,1]
  raw <- raw[,-1]
  file = paste(sample, '.pdf', sep = "", collapse = NULL)
  pdf(width = 15, height= 15, file = file)
  pheatmap(t(t(raw)), cluster_cols=F, cellwidth = 20, cellheight = 20,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),)
  dev.off()
}

raw <- sig_raw[sig_raw$C>0.8,]
anno_raw <- raw[,c(1,9)]
raw <- raw[,c(1,2,3,4,5)]
raw <- raw[,c(1,6,7,8)]
rownames(raw) = raw[,1]
raw <- raw[,-1]

annotation_col <- data.frame( Type = factor( anno_raw$type ) )
rownames( annotation_col ) <- anno_raw$gene

unique(anno_raw$type)
col = list( Type = 
              c(intron = '#3ac0c4', 
                exon = '#5e50a1', 
                promoter = '#c276b1',
                TSS = '#79c9f0',
                utr5 = '#fdcab5',
                utr3 = '#fc8a6a',
                TES = '#fcbe6f')
            )

a=scale(raw)
pheatmap( fontsize = 8, fontsize_col = 6,
          color = colorRampPalette(colors = c("#4174b4","white","#b21f48"))(100),
          scale(raw), 
          annotation_row = annotation_col, 
          annotation_colors = col,
          show_rownames = F, show_colnames = T,
          annotation_legend = T, cluster_cols = F)
dev.off()

