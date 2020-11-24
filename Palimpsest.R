# Step0 Before starting your project --------------------------------------
#
inp("BiocManager")
#
dependencies=c('data.table', 'zip', 'XML', 'dplyr', 'NMF', 'GenomicRanges', 
               'reshape2', 'BSgenome.Hsapiens.UCSC.hg19', 'rtracklayer',
               'plyr', 'gtools', 'lsa', 'gplots', 'RCircos', 'plotrix', 'ggplot2',
               'DelayedArray', 'dbplyr', 'curl', 'BSgenome', 'readxl', 'carData',
               'Rgraphviz', 'Biostrings', 'VariantAnnotation', 'Biobase', 'stringi',
               'registry', 'doParallel', 'spgs', 'ggpubr', 'GenomeInfoDbData','pkgmaker')
inp(dependencies)

install.packages("./Palimpsest-master/",repos=NULL,type="source")
install.packages("../BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz", repos = NULL)


rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()


# Step1 Starting your project -------------------------------------------

library("Palimpsest")
library("BSgenome.Hsapiens.UCSC.hg19")

resdir_parent <- "./output/Palimpsest/"
if(!file.exists(resdir_parent))	dir.create(resdir_parent)
resdir <- paste0(resdir_parent,"SBS_denovo/")
if(!file.exists(resdir))	dir.create(resdir)

#load("./vcf.RData")
vcf=read.table(file = "./input/Palimpsest_SNV_INS_DEL.vcf", 
               header= T, fill = T, sep = "\t")

vcf <- annotate_VCF(vcf = vcf, ref_genome = BSgenome.Hsapiens.UCSC.hg19,
                    ref_fasta = "/Genomes/Homo_sapiens_assembly19.fasta",
                    add_ID_cats = FALSE, add_DBS_cats = F)

SBS_input <- palimpsest_input(vcf = vcf, Type = "SBS")
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input,
                                  resdir = resdir)

# Compare the de novo signatures with published COSMIC signatures
compare_tab <- compare_results(reference_sigs = SBS_cosmic, extraction_1 = SBS_denovo_sigs)
compare_tab
readr::write_delim(compare_tab,path = paste0(resdir,"Comparison_table.txt"))

pdf(file.path(resdir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs,SBS_cosmic)
dev.off()

# Define signature colours for plotting
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))

# Calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input,
                                        input_signatures = SBS_denovo_sigs,
                                        threshold = 8,
                                        resdir = resdir, 
                                        signature_colours = SBS_col,
                                        input_vcf = vcf)

pdf(file.path(resdir,"signature_content_plot.pdf"),width=15,height=10)
deconvolution_exposure(signature_colours = SBS_col, signature_contribution = SBS_signatures_exp)
dev.off()


#-------------------------------------------------------------------------------------------------
# 3] Extract with published SBS COSMIC signatures
#-------------------------------------------------------------------------------------------------
resdir <- paste0(resdir_parent,"SBS_COSMIC_Extraction/");if(!file.exists(resdir))	dir.create(resdir) 

# select desired COSMIC SBS reference signatures 
#SBS_liver_names <- 

SBS_liver_names = compare_tab[,1]
SBS_liver_names = c("SBS3","SBS40","SBS6","SBS5","SBS39","SBS1","SBS25","SBS15","SBS42")

## generate colours for new signatures
 for(new_name in c(SBS_liver_names[SBS_liver_names %!in% names(sig_cols)]))
  sig_cols[new_name] <- signature_colour_generator(new_name) 

SBS_liver_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_liver_names,]

# calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, 
                                        input_signatures = SBS_liver_sigs,
                                        threshold = 8, 
                                        signature_colours = sig_cols,
                                        resdir = resdir, 
                                        input_vcf = vcf)

pdf(file.path(resdir,"SBS_signature_content_plot.pdf"),width=13,height=10)
deconvolution_exposure(signature_contribution = SBS_signatures_exp,signature_colours = sig_cols)
dev.off()

#-------------------------------------------------------------------------------------------------
# 6] Assign the probability of each individual mutation being due to each process 
# (Again, this function works for DBS, Indel and SV signatures too)
#-------------------------------------------------------------------------------------------------
# This step is quite computationally-intensive for whole genome data. 
# For this example we restrict the analysis to coding mutations 
vcf.cod <- vcf[(!is.na(vcf$Driver) & vcf$Type=="SNV"),]
vcf.cod <- signature_origins(input = vcf.cod, Type = "SBS",
                             input_signatures = SBS_liver_sigs,
                             signature_contribution = SBS_signatures_exp)
#write.table(vcf.cod, file = "all.xls")

# Estimate and represent the cumulative contribution of signatures to each driver gene
drivers <- c("CHD8","CHEK2",
             "CREBBP","MACF1","MED12","PLCG1","RPL10","TET2","KDM6A","KMT2D","TSC2","HUWE1",
             "KALRN","SYNE1","KMT2C","EP300","DMD","NOTCH2","NOTCH1","TP53")
drivers <- c("FANCA","NTRK1","CREBBP","ELL","DST","KIAA1549","WIF1","ITGA6","ZMYND8","KMT2D","SS18","PLCG1","TAF1","ETV1","BCR","RALGDS","ATIC","TSC2","TP53","FN1","SYNE1","DOCK2","LZTR1","FBN2","ABI1","PCSK7","TSHZ2","NEDD4L","CACNA1D","PTPN13","AHNAK","KMT2A","MDN1","PCM1","TSC1","ATP2B3","AXIN2","ARID1B","ABL1","FBXO11","MED12","SSX2IP","KRAS","NAV3","RAG1","CRTC3","ERBB2IP","ATM","SMO","GOPC","EPPK1","CDKN2A","CBFB","AKAP9","RUNX1","SH3PXD2A","MUC1")
matprob <- matrix(nrow=length(drivers),
                  ncol=length(SBS_liver_names),
                  dimnames=list(drivers, SBS_liver_names))
sig.cols <- paste0(rownames(SBS_liver_sigs),".prob")#grep("prob",colnames(vcf.cod))
for(i in 1:nrow(matprob)){
  g <- rownames(matprob)[i]
  ind <- which(vcf.cod$gene_name==g)
  matprob[i,] <- apply(vcf.cod[ind,sig.cols],2,sum,na.rm=T)
}
barplot(t(matprob),col = sig_cols,border = sig_cols,las=2)
legend("top",names(sig_cols)[names(sig_cols) %in% rownames(SBS_liver_sigs)],fill=sig_cols,ncol=5,
       cex=0.75,bty ="n",inset = c(0,-0.3),xpd = T)

# Compare signature 16 contribution between CTNNB1 mutations and others
library(ggplot2)
gene=c("SYNE1","TP53","KRAS")
for (genename in drivers) {
  vcf.cod$genename <- (vcf.cod$gene_name==genename)
  vcf.cod[is.na(vcf.cod$gene_name),genename] <- FALSE
  #p_value=wilcox.test(c(vcf.cod$SBS6.prob) ~ c(vcf.cod$genename))[3]$p.value
  #if(p_value <= 0.05){
  #  print(genename)
  #  print(p_value)
  #}
}
genename="TP53"
vcf.cod$genename <- (vcf.cod$gene_name==genename)
vcf.cod[is.na(vcf.cod$gene_name),genename] <- FALSE
ggplot(vcf.cod, 
       aes(x = genename,y = SBS1.prob,
           color = genename,
           fill = genename)) +
  geom_violin(alpha = 0.4,adjust = 1)
wilcox.test(c(vcf.cod$SBS1.prob) ~ c(vcf.cod$genename))

# Add signature probability columns to the original vcf
vcf <- merge(vcf,vcf.cod,all=TRUE,sort=FALSE)


#-------------------------------------------------------------------------------------------------
# 7] Clonality analysis
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Clonality");if(!file.exists(resdir)){dir.create(resdir)}

# Load copy number analysis (CNA) and annotation (annot) data
cna_data = read.table(file = "./input/Palimpsest_cnv.xls", 
                      header= T, fill = T, sep = "\t")
#cna_data = load2object("cna_data.RData")
annot = read.table(file = "./input/Palimpsest_purity.list", 
                   header= T, fill = T, sep = "\t")
#annot = load2object("annot_data.RData")

# Calculate the Cancer Cell Fraction (CCF) of each mutation.
vcf_cna <- cnaCCF_annot(vcf= vcf, 
                        annot_data = annot,
                        cna_data = cna_data, 
                        CCF_boundary = 0.95)

# Generate graphical representations of clonality analysis
cnaCCF_plots(vcf=vcf_cna,resdir=resdir)


#-------------------------------------------------------------------------------------------------
# 8] Compare mutational signatures between early clonal and late subclonal mutations in each tumour
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Signatures_early_vs_late/");if(!file.exists(resdir)){dir.create(resdir)}

# Estimate the contribution of each signature to clonal and subclonal mutations in each tumour
vcf.clonal <- vcf_cna[which(vcf_cna$Clonality=="clonal"),]
SBS_input_clonal <- palimpsest_input(vcf = vcf.clonal,Type = "SBS")
signatures_exp_clonal <- deconvolution_fit(input_matrices = SBS_input_clonal,
                                           input_signatures = SBS_liver_sigs, resdir =  resdir,
                                           save_signatures_exp = T)

vcf.subclonal <- vcf_cna[which(vcf_cna$Clonality=="subclonal"),]
SBS_input_subclonal <- palimpsest_input(vcf = vcf.subclonal,Type = "SBS")
signatures_exp_subclonal <- deconvolution_fit(input_matrices = SBS_input_subclonal,
                                              input_signatures = SBS_liver_sigs, resdir =  resdir,
                                              save_signatures_exp = F)

# Generate per tumour comparisons of clonal and subclonal mutations
palimpsest_DissectSigs(vcf=vcf_cna, 
                       signatures_exp_clonal = signatures_exp_clonal,
                       signatures_exp_subclonal = signatures_exp_subclonal,
                       sig_cols = sig_cols,
                       resdir=resdir)

# Generate across the series comparisons of signature assigned to clonal and subclonal mutations
palimpsest_clonalitySigsCompare(clonsig = signatures_exp_clonal$sig_nums,
                                subsig = signatures_exp_subclonal$sig_nums, 
                                msigcol = sig_cols, 
                                resdir = resdir)
dev.off()

write.table(signatures_exp_clonal$sig_nums, file = "clonal.txt")
write.table(signatures_exp_subclonal$sig_nums, file = "subclonal.txt")
#-------------------------------------------------------------------------------------------------
# 9] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"ChromosomeDups_timing/");if(!file.exists(resdir)){dir.create(resdir)}

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf_cna,cna_data = cna_data,cyto=cytoband_hg19)
vcf_cna <- chrom_dup_time$vcf
point.mut.time <- chrom_dup_time$point.mut.time
cna_data <- chrom_dup_time$cna_data

# Visualising timing plots
chrTime_plot(vcf = vcf_cna, point.mut.time = point.mut.time, resdir = resdir,cyto = cytoband_hg19)



#-------------------------------------------------------------------------------------------------
# 11] Visualise the natural history of tumour samples:
#-------------------------------------------------------------------------------------------------

resdir <- file.path(resdir_parent,"Natural_history/");if(!file.exists(resdir)){dir.create(resdir)}

vcf_cna1 = vcf_cna[!vcf_cna$Sample %in% c("P182_5","P192_3"), ]
palimpsest_plotTumorHistories(
  vcf = vcf_cna1, 
  #sv.vcf = SV.vcf,
  cna_data =  cna_data, 
  point.mut.time = point.mut.time, 
  clonsig = signatures_exp_clonal$sig_props, 
  subsig = signatures_exp_subclonal$sig_props, 
  msigcol = sig_cols, 
  #msigcol.sv = SV_cols, 
  resdir = resdir)

write.table(vcf_cna, file = "vcf_cnv.xls")

