# Step0 Before starting your project --------------------------------------

#installPackage("PoisonAlien/maftools")
#installPackage("mclust")

rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()

# Step1 Starting your project -------------------------------------------

maf <- data.table::as.data.table(
  read.csv(file = "./Somatic_mutation_p104_s274.xls",
           header = TRUE, sep = '\t',
           stringsAsFactors = FALSE,
           comment.char = "#"))
maf <- maf[maf$Tumor_Sample_Barcode != "P2018_14791_0", ]
maf <- maf[maf$Tumor_Sample_Barcode != "P2018_759_4", ]
names(maf)[38] = "Protein_Change"
length(unique(maf$Tumor_Sample_Barcode))

unique(maf$Variant_Classification)

phenoData <- read.table("all.group", sep = '\t', 
                        header = T)
names(phenoData)[1] = "Tumor_Sample_Barcode"


library("maftools") 
laml <- read.maf(maf, 
                 clinicalData = phenoData, 
                 removeDuplicatedVariants = F,
                 useAll = T,
                 cnLevel = "all")

col = c('Missense_Mutation' = '#5076ad', 
        'In_Frame_Del' = '#aacd6f', 
        'In_Frame_Ins' = '#018574', 
        'Nonstop_Mutation' = '#9b8dc4',
        'Splice_Site' = '#e3e369',
        'Frame_Shift_Del' = '#ff9800', 
        'Frame_Shift_Ins' = '#78ccbd', 
        'Nonsense_Mutation' = '#c65148',
        'Multi_Hit' = '#f26b71')

col = c('Missense_Mutation' = '#2b83bc', 
        'Nonsense_Mutation' = '#e04f4a',
        'In_Frame_Del' = '#3ac0c4', 
        'In_Frame_Ins' = '#018574', 
        'Nonstop_Mutation' = '#9b8dc4',
        'Splice_Site' = '#6dae8e',
        'Frame_Shift_Del' = '#ffe694', 
        'Frame_Shift_Ins' = '#ff9800', 
        'Multi_Hit' = '#f37258')


pdf("plotmafSummary.pdf", width = 10, height = 7)
plotmafSummary(maf = laml,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE,
               log_scale = FALSE,
               showBarcodes = TRUE,
               textSize = 0.4,
               color = col,
               fs = 1.1,
               titleSize = c(1.2, 0.8))
dev.off()

phecolors = list( Sample = 
                    c("Paired tumor/normal" = '#63b8ff',
                      "Unpaired tumor" = '#000000'),
                  Type = 
                    c("Cohort−1" = '#000080',
                      "Cohort−2" = '#4169e1',
                      "Cohort−3" = '#ff69b4',
                      "Cohort−4" = '#ff0000')
                  
)
phecolors = list( KRAS.mutation = 
                    c(none = '#fafafa',
                      hit = '#bc141a'),
                  TP53.mutation = 
                    c(none = '#fafafa',
                      hit = '#36a12e'),
                  KRAS.type =
                    c(G12D = '#f14432',
                      G12V = '#61409b',
                      none = '#fafafa'),
                  type =
                    c(Paired = '#5ac8dc',
                      Single = '#f48074')
                  
)

top_gene = as.character(unlist(getGeneSummary(laml)[1:10,1]))
#aml_genes_vaf = subsetMaf(maf = laml,
#                          genes = top_gene,
#                          fields = "t_AF",
#                          mafObj = FALSE)[,mean(t_AF, na.rm = TRUE), Hugo_Symbol]
#colnames(aml_genes_vaf)[2] = "VAF"
#head(aml_genes_vaf)'
gene = c("KMT2D","SETD2","NBPF10","NRAS","KRAS","PAX5","FLT3","TP53","NF1","ASXL1","NBPF8","DIAPH1","PRB2","TSFM","MYB","MYC","ETV6","PTPN11","BCORL1","CDKN2A","XBP1","PRSS3","ID3","NBPF11","CAMK2N1","KMT5A","LRP8","PRH2","MAGEA12","SPAG4")
sample_order=c("P2016_26573_3","P2016_28449_0","P2018_1114_2","P2018_12531_7","P2018_16865_7","P2018_3672_13","P2018_6502_2","P2019_173282_2","P2019_8483_2","P2018_4835_7","P2018_14265_10","P2017_3831_0","P2018_14791_0","P2018_4091_1","P2019_4013_3","P2019_11362_7","P2017_9543_2","P2019_14655_3","P2016_24152_2","P2017_5610_3","P2017_13066_2","P2019_847_0","P2018_12659_1","P2018_4710_3","P2018_6545_4","P2019_9871_2","P2018_F1372_6","P2018_19742_0","P2017_10651_1","P2018_759_4","P2017_20985_2","P2018_11720_1","P2017_6822_0","P2018_9982_2","P2016_22804_5","P2018_17522_5","P2018_860_1","P2019_13903_4","P2019_15758_3","P2019_15172_7","P2017_15491_2","P2019_14478_0")
pdf("oncoplot.pdf", width = 10, height = 7)
oncoplot(  maf = laml,
           top = 30,
           #minMut = NULL,
           genes = gene,
           altered = FALSE,
           drawRowBar = TRUE,
           drawColBar = TRUE,
           #leftBarData = aml_genes_vaf,
           #leftBarLims = c(0, 0.5),
           #rightBarData = NULL,
           rightBarLims = c(0, 50),
           #topBarData = NULL,
           logColBar = FALSE,
           includeColBarCN = TRUE,
           clinicalFeatures = c("Type","Sample"),
           annotationColor = phecolors,
           #pathways = NULL,
           #selectedPathways = NULL,
           draw_titv = FALSE,
           showTumorSampleBarcodes = FALSE,
           barcode_mar = 6,
           barcodeSrt = 90,
           gene_mar = 7,
           anno_height = 1,
           legend_height = 3,
           groupAnnotationBySize = FALSE,
           sortByAnnotation = TRUE,
           annotationOrder = c("Cohort−1","Cohort−2","Cohort−3","Cohort−4"),
           sortByMutation = FALSE,
           keepGeneOrder = FALSE,
           #GeneOrderSort = TRUE,
           #sampleOrder = sample_order,
           #additionalFeature = NULL,
           #additionalFeaturePch = 20,
           #additionalFeatureCol = "gray70",
           #additionalFeatureCex = 0.9,
           #genesToIgnore = NULL,
           removeNonMutated = FALSE,
           fill = TRUE,
           #cohortSize = NULL,
           colors = col,
           #bgCol = "#fafafa",
           bgCol = "#cccccc",
           borderCol = "white",
           #annoBorderCol = NA,
           #numericAnnoCol = NULL,
           drawBox = FALSE,
           fontSize = 0.5,
           SampleNamefontSize = 0.85,
           titleFontSize = 1.5,
           legendFontSize = 1,
           annotationFontSize = 1,
           sepwd_genes = 0.6,
           sepwd_samples = 0.25,
           writeMatrix = FALSE,
           colbar_pathway = FALSE,
           showTitle = TRUE,
           titleText = NULL
)
dev.off()


laml.titv = titv(maf = laml, plot = TRUE, useSyn = TRUE)

pdf("plotTiTv.pdf", width = 10, height = 7)
plotTiTv(res = laml.titv,
         showBarcodes = F,
         textSize = 0.8,
         baseFontSize = 1.3,
         axisTextSize = c(1, 1),
         plotNotch = FALSE
         )
dev.off()


library("data.table")
source(file = "lollipopPlot_new.R")
top_gene

for (i in 1:length(top_gene)) {
  file=paste0(top_gene[i], "_lollipopPlot.pdf", collapse = "_")
  pdf(file, width = 16, height = 6)
  lollipopPlot(maf = laml,
               gene = top_gene[i],
               AACol = 'Protein_Change',
               labelPos = "all",
               labPosSize = 0.9,
               showMutationRate = TRUE,
               showDomainLabel = TRUE,
               cBioPortal = FALSE,
               roundedRect = TRUE,
               repel = F,
               collapsePosLabel = TRUE,
               showLegend = TRUE,
               legendTxtSize = 0.8,
               labPosAngle = 30,
               domainLabelSize = 0.8,
               axisTextSize = c(1, 1),
               printCount = TRUE,
               colors = col,
               domainAlpha = 1,
               domainBorderCol = "grey",
               bgBorderCol = "grey",
               labelOnlyUniqueDoamins = FALSE,
               defaultYaxis = TRUE,
               titleSize = c(1, 1),
               pointSize = 1.5)
  dev.off()
}



top_sample = as.character(unlist(getSampleSummary(laml)[,1]))
for (i in 1:length(top_sample)) {
  print(i)
  rainfallPlot(maf = laml,
               tsb = top_sample[i],
               detectChangePoints = TRUE,
               width = 10,
               height = 6,
               fontSize = 1,
               pointSize = 0.6,
               savePlot = TRUE)
}

pdf("tcgaCompare.pdf", width = 10, height = 7)
laml.mutload = tcgaCompare(maf = laml,
                           primarySite = FALSE,
                           col = c("gray70", "black"),
                           bg_col = c("#EDF8B1", "#2C7FB8"),
                           medianCol = "red",
                           decreasing = FALSE,
                           logscale = TRUE,
                           rm_hyper = TRUE,
                           rm_zero = TRUE,
                           cohortName = 'PAAD-42')
dev.off()


plotVaf(maf = laml,
        top = 10,
        orderByMedian = TRUE,
        keepGeneOrder = FALSE,
        flip = FALSE,
        gene_fs = 0.8,
        axis_fs = 0.8,
        showN = T,
        height = 6,
        width = 12,
        fn = "vaf_plot")



all.lesions <- "ALL.all_lesions.conf_99.xls.txt"
amp.genes <- "ALL.amp_genes.conf_99.xls.txt"
del.genes <- "ALL.del_genes.conf_99.xls.txt"
scores.gistic <- "ALL.scores.xls.gistic"
laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, 
                         gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, 
                         gisticScoresFile = scores.gistic)
c_col = c('Amp' = '#a73c52', 
          'Del' = '#018574')

markBands = laml.gistic@cytoband.summary[laml.gistic@cytoband.summary$qvalues<0.01,]$Cytoband
pdf("gisticChromPlot.pdf", width = 14, height = 5)
gisticChromPlot(gistic = laml.gistic,
                fdrCutOff = 0.01,
                color = c_col,
                markBands = markBands,
                ref.build = "hg19",
                cytobandOffset = 0.03,
                txtSize = 0.7,
                cytobandTxtSize = 0.8,
                maf = laml,
                mutGenes = top_gene,
                #y_lims = NULL,
                mutGenesTxtSize = 0.8)
dev.off()

pdf("gisticBubblePlot.pdf", width = 12, height = 10)
gisticBubblePlot(gistic = laml.gistic,
                 color = c_col,
                 markBands = markBands,
                 fdrCutOff = 0.01,
                 log_y = TRUE,
                 txtSize = 0.1)
dev.off()

pdf("gisticOncoPlot.pdf", width = 10, height = 5)
gisticOncoPlot(gistic = laml.gistic,
               #clinicalData = getClinicalData(x = laml),
               #clinicalFeatures = c("T","N","M"),
               #sortByAnnotation = FALSE,
               top = 10,
               #showTumorSampleBarcodes = TRUE,
               gene_mar = 7,
               barcode_mar = 6,
               sepwd_genes = 0.6,
               sepwd_samples = 0.25,
               #sampleOrder = NULL,
               annotationColor = phecolors,
               removeNonAltered = FALSE,
               colors = c_col,
               SampleNamefontSize = 0.85,
               fontSize = 0.55,
               legendFontSize = 1,
               annotationFontSize = 1,
               borderCol = "white",
               bgCol = "#fafafa")
dev.off()


laml <- read.maf(maf, 
                 clinicalData = phenoData, 
                 removeDuplicatedVariants = F,
                 useAll = T,
                 cnLevel = "all",
                 gisticAllLesionsFile = all.lesions, 
                 gisticAmpGenesFile = amp.genes, 
                 gisticDelGenesFile = del.genes, 
                 gisticScoresFile = scores.gistic)

col = c('Missense_Mutation' = '#2b83bc', 
        'Nonsense_Mutation' = '#e04f4a',
        'In_Frame_Del' = '#3ac0c4', 
        'In_Frame_Ins' = '#018574', 
        'Nonstop_Mutation' = '#9b8dc4',
        'Splice_Site' = '#6dae8e',
        'Frame_Shift_Del' = '#ffe694', 
        'Frame_Shift_Ins' = '#ff9800', 
        'Multi_Hit' = '#f37258',
        'Amp' = '#a73c52', 
        'Del' = '#018574')

pdf("gistic_all_OncoPlot.pdf", width = 10, height = 10)
oncoplot(  maf = laml,
           top = 30,
           #minMut = NULL,
           #genes = NULL,
           altered = FALSE,
           drawRowBar = TRUE,
           drawColBar = TRUE,
           #leftBarData = aml_genes_vaf,
           #leftBarLims = c(0, 0.5),
           #rightBarData = NULL,
           #rightBarLims = NULL,
           #topBarData = NULL,
           logColBar = FALSE,
           includeColBarCN = TRUE,
           clinicalFeatures = c("KRAS.type","KRAS.mutation","TP53.mutation","type"),
           annotationColor = phecolors,
           #pathways = NULL,
           #selectedPathways = NULL,
           draw_titv = FALSE,
           showTumorSampleBarcodes = TRUE,
           barcode_mar = 6,
           barcodeSrt = 90,
           gene_mar = 5,
           anno_height = 1,
           legend_height = 3,
           sortByAnnotation = TRUE,
           groupAnnotationBySize = FALSE,
           #annotationOrder = NULL,
           sortByMutation = FALSE,
           keepGeneOrder = FALSE,
           #GeneOrderSort = TRUE,
           #sampleOrder = NULL,
           #additionalFeature = NULL,
           #additionalFeaturePch = 20,
           #additionalFeatureCol = "gray70",
           #additionalFeatureCex = 0.9,
           #genesToIgnore = NULL,
           removeNonMutated = FALSE,
           fill = TRUE,
           #cohortSize = NULL,
           colors = col,
           bgCol = "#fafafa",
           borderCol = "white",
           #annoBorderCol = NA,
           #numericAnnoCol = NULL,
           drawBox = FALSE,
           fontSize = 0.5,
           SampleNamefontSize = 0.85,
           titleFontSize = 1.5,
           legendFontSize = 1,
           annotationFontSize = 1,
           sepwd_genes = 0.6,
           sepwd_samples = 0.25,
           writeMatrix = FALSE,
           colbar_pathway = FALSE,
           showTitle = TRUE,
           titleText = NULL
)
dev.off()


pdf("somaticInteractions.pdf", width = 10, height = 10)
somaticInteractions(
  maf = laml,
  top = 25,
  #genes = NULL,
  pvalue = c(0.05, 0.01),
  returnAll = TRUE,
  #geneOrder = NULL,
  fontSize = 0.5,
  showSigSymbols = FALSE,
  showCounts = TRUE,
  countStats = "all",
  countType = "all",
  countsFontSize = 0.8,
  countsFontColor = "black",
  colPal = "BrBG",
  showSum = TRUE,
  colNC = 9,
  nShiftSymbols = 5,
  sigSymbolsSize = 2,
  sigSymbolsFontSize = 0.9,
  pvSymbols = c(46, 42),
  limitColorBreaks = TRUE
)

dev.off()


pathway = OncogenicPathways(maf = laml, fontSize = 1, panelWidths = c(2, 4, 4))
for (i in pathway$Pathway) {
  file = paste('Pathway_', i, '.pdf', sep = "", collapse = NULL)
  pdf(width = 9.65, height= 5.76, file = file)
  PlotOncogenicPathways(maf = laml,
                        pathways = i,
                        fullPathway = FALSE,
                        removeNonMutated = TRUE,
                        tsgCol = "red",
                        ogCol = "blue",
                        fontSize = 0.6,
                        showTumorSampleBarcodes = TRUE,
                        sampleOrder = NULL,
                        SampleNamefontSize = 0.6)
  dev.off()
}


sample = as.character(unlist(getSampleSummary(laml)[,1]))
sample = sample[-match("P2018_860_1", sample)]
sample = sample[-match("P2018_11720_1", sample)]
sample = sample[-match("P2019_13903_4", sample)]
sample = sample[-match("P2018_F1372_6", sample)]

sample = sample[-match("P2018_860_1", sample)]
sample = sample[-match("P2018_9982_2", sample)]
sample = sample[-match("P2019_13903_4", sample)]
sample = sample[-match("P2018_11720_1", sample)]
sample = sample[-match("P2018_F1372_6", sample)]

for (i in sample) {
  het = inferHeterogeneity(maf = laml, 
                           tsb = i,
                           vafCol = "t_AF",
                           segFile = "CNV.segment.txt",
                           useSyn = TRUE)
  file = paste('Het_', i, '_N.pdf', sep = "", collapse = NULL)
  pdf(width = 9.65, height= 5.76, file = file)
  plotClusters(clusters = het,
               showCNvars = T, 
               tsb = i, 
               genes = "CN_altered")
  dev.off()
}


library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml,
                               prefix = 'chr',
                               add = TRUE,
                               ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                               ignoreChr = NULL,
                               useSyn = TRUE,
                               fn = NULL)
plotApobecDiff(tnm = laml.tnm,
               maf = laml,
               pVal = 1,
               title_size = 1,
               axis_lwd = 1,
               font_size = 1.2)

library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm,
                               nMin = 2,
                               nTry = 6,
                               nrun = 10,
                               parallel = 4,
                               pConstant = NULL,
                               verbose = TRUE,
                               plotBestFitRes = FALSE)
laml.sig = extractSignatures(mat = laml.tnm, n = 4,
                             plotBestFitRes = FALSE,
                             parallel = 4,
                             pConstant = NULL)
laml.og30.cosm = compareSignatures(nmfRes = laml.sig,
                                   sig_db = "legacy")
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities,
                   cluster_rows = FALSE,
                   main = "cosine similarity against validated signatures")
pheatmap::pheatmap(mat = laml.v3.cosm$cosine_similarities,
                   cluster_rows = FALSE,
                   main = "cosine similarity against validated signatures v3")

titv_col = c('C>T' = '#f6695e', 
        'T>C' = '#ff9508',
        'C>A' = '#4dabf5', 
        'C>G' = '#6574c4', 
        'T>G' = '#523810',
        'T>A' = '#70bf73')
maftools::plotSignatures(nmfRes = laml.sig,
                         title_size = 0.8,
                         sig_db = "legacy",
                         contributions = FALSE,
                         color = titv_col,
                         patient_order = NULL,
                         font_size = 1.2,
                         show_title = TRUE,
                         axis_lwd = 2,
                         show_barcodes = FALSE,
                         yaxisLim = 0.3)
maftools::plotSignatures(nmfRes = laml.sig,
                         title_size = 0.8,
                         sig_db = "SBS",
                         contributions = FALSE,
                         color = titv_col,
                         patient_order = NULL,
                         font_size = 1.2,
                         show_title = TRUE,
                         axis_lwd = 2,
                         show_barcodes = FALSE,
                         yaxisLim = 0.3)

library("barplot3d")
sig = laml.sig$signatures[,4]
barplot3d::legoplot3d(contextdata = sig, labels = FALSE,
                      scalexy = 0.01, sixcolors = "sanger",
                      alpha = 0.5, gap = 0.2,
                      theta = 60, phi = 40, gridlines = FALSE,
                      zlabels = TRUE, zsub = "Probability")

laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig,
                              minMut = 5, useCNV = FALSE, fn = NULL)

plotEnrichmentResults(enrich_res = laml.se, pVal = 0.1,
                      cols = NULL,
                      annoFontSize = 0.8,
                      geneFontSize = 0.8,
                      legendFontSize = 0.8,
                      showTitle = TRUE)
