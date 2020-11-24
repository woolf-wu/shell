
# NO RUN ------------------------------------------------------------------

# rm(list = objects( all = TRUE ))
# 
# if (!is.null( dev.list() )) dev.off()
# 
# sysPackages <- (.packages())
# 
# if ( capabilities( "libcurl" ) == T ) {
#   options( download.file.method = "libcurl" )
# }
# 
# options( stringsAsFactors = FALSE )
# 
# Sys.setlocale("LC_ALL","English")
# 
# .libPaths( c( "G:/R-packages", "C:/Program Files/R/R-3.5.2/library") )
# 
# local({
#   options( repos  = "https://mirrors.ustc.edu.cn/CRAN/" )
#   options( BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/" )
# })
# 
# 
# bioPackages <-
#   c(
#     "RCurl", "XML", "RSelenium", ## step1
#     "R.utils", "data.table", ## step2
#     "GEOquery", ## download
#     "FactoMineR", "factoextra", "ggfortify", ## PCA
#     "pheatmap", ## heatmap
#     "ggplot2", ## Volcano plot
#     "limma", "DESeq2", "edgeR", ## DEG
#     "clusterProfiler", "org.Hs.eg.db", "stringr", ## step3
#     "pathview", ## kegg
#     "ggpubr", ## plot
#     "survival", "survminer", ## OS
#     "deconstructSigs", "BSgenome", "BSgenome.Hsapiens.UCSC.hg38"
#   )
# 
# 
# lapply( bioPackages,
#         function(bioPackage) {
#           if ( !require( bioPackage, character.only = T ) ) {
#             CRANpackages <- available.packages()
# 
#             ## install packages by CRAN
#             if ( bioPackage %in% rownames( CRANpackages ) ) {
#               install.packages( bioPackage )
# 
#             }else{
#               ## install packages by bioconductor
#               ## R version >= 3.5 ===> BiocManager
#               if ( as.character( sessionInfo()$R.version$minor ) >= 3.5 ) {
#                 if (!requireNamespace("BiocManager", quietly = TRUE))
#                   install.packages("BiocManager")
#                 BiocManager::install(bioPackage, update = TRUE, ask = FALSE)
# 
#               }else{
#                 ## R version < 3.5 ===> BiocInstaller
#                 if (!requireNamespace("BiocInstaller", quietly = TRUE))
#                   source( "https://bioconductor.org/biocLite.R" )
#                 BiocInstaller::biocLite( bioPackage, ask = FALSE)
#               }
#             }
#           }
#         }
# )




# RUN ---------------------------------------------------------------------

tumor_type <- "PAAD"

suppressMessages(library(RCurl))
suppressMessages(library(XML))
suppressMessages(library(RSelenium))

pars_subpage <- function(type) {
  webElems <- remDr$findElements(using = "css selector", "a")
  sub_type <- unlist(lapply(webElems, function(x) {x$getElementText()}))
  targetElem <- webElems[[grep(type, sub_type)]]
  targetElem$clickElement()
  Sys.sleep(3)
  webElems <- remDr$findElements(using = "css selector", "a")
  sub_url <- unlist(lapply(webElems, function(x) {x$getElementText()}))
  url <- sub_url[grep("http", sub_url)[1]]
  return(url)
}

system(paste('"C:\\Program Files\\Java\\jdk-10.0.2\\bin\\java.exe"',
             '-Dwebdriver.chrome.driver="C:\\Program Files (x86)\\Google\\Chrome\\Application\\chromedriver.exe"',
             '-jar "C:\\Program Files\\Java\\selenium-server-standalone-3.9.1.jar"'),
       wait = FALSE)

url <- 'https://xenabrowser.net/datapages/'
remDr <- remoteDriver(browserName = "chrome",
                      remoteServerAddr = "localhost",
                      port = 4444L)
remDr$open()
remDr$navigate(url)

## GDC
webElems <- remDr$findElements(using = "css selector", "a")
sub_type <- unlist(lapply(webElems, function(x) {x$getElementText()}))
targetElem <- webElems[[grep("GDC Hub", sub_type)]]
targetElem$clickElement()

webElems <- remDr$findElements(using = "css selector", "a")
sub_type <- unlist(lapply(webElems, function(x) {x$getElementText()}))
targetElem <- webElems[[grep(tumor_type, sub_type)]]
targetElem$clickElement()

phenoData_url <- pars_subpage("Phenotype")
remDr$goBack()
raw_data_url <- pars_subpage("Counts")
remDr$close()

save(phenoData_url, raw_data_url, file = "./data/url.Rdata")


# NO RUN ------------------------------------------------------------------

if (!dir.exists('./raw_data/')) {
  dir.create("./raw_data/")
}

if (!dir.exists('./data/')) {
  dir.create("./data/")
}

if (!dir.exists('./fig/')) {
  dir.create("./fig/")
}



# RUN ---------------------------------------------------------------------

load("./data/url.Rdata")
Rdata_file <- paste( './data/TCGA.', tumor_type, '.counts.Rdata', sep = '' )
Rdata_file
if (!file.exists( Rdata_file )) {
  url <- raw_data_url
  url_word <- strsplit( url, split = "/" )
  url_word <- url_word[[1]]
  url_length <- length( url_word )
  gz_file <- url_word[ url_length ]
  destfile <- paste('./raw_data/', gz_file, sep = '')
  destfile
  download.file(url, destfile)
  
  suppressMessages(library(R.utils))
  gunzip(destfile, remove = F)
  
  suppressMessages(library(data.table))
  raw_data <- fread( sub(".gz", "", destfile), sep = '\t', header = T)
  raw_data <- as.data.frame( raw_data )
  raw_data[1:5, 1:5] 
  rownames( raw_data ) <- raw_data[, 1]
  raw_data <- raw_data[, -1]
  raw_data[1:5, 1:5]
  raw_data <- 2^raw_data - 1
  raw_data <- ceiling( raw_data )
  raw_data[1:5, 1:5]
  dim( raw_data )
  save( raw_data, file = Rdata_file )
}else{
  load( Rdata_file )
}


Rdata_file <- paste('./data/TCGA.', tumor_type, '.phenoData.Rdata', sep = '')
Rdata_file
if (!file.exists( Rdata_file )) {
  url <- phenoData_url
  url_word <- strsplit(url, split = "/")
  url_word <- url_word[[1]]
  url_length <- length( url_word )
  gz_file <- url_word[ url_length ]
  destfile <- paste('./raw_data/', gz_file, sep = '')
  destfile
  download.file(url, destfile)
  
  phenoData <- read.table( destfile,
                           header = T,
                           sep = '\t',
                           quote = '' )
  rownames( phenoData ) <- phenoData[ , 1]
  ## name <- gsub(pattern = "-", replacement = ".", name)
  colnames( phenoData )[1] <- "Tumor_Sample_Barcode"
  phenoData[1:5, 1:5]
  save( phenoData, file = Rdata_file )
}else{
  load( Rdata_file )
}




# RUN ---------------------------------------------------------------------

Rdata_file <- paste( './data/TCGA.', tumor_type, '.counts.Rdata', sep = '' )
Rdata_file
load( Rdata_file )


Rdata_file <- paste('./data/TCGA.', tumor_type, '.phenoData.Rdata', sep = '')
Rdata_file
load( Rdata_file )


head(phenoData)

# lymph_node_examined_count
# neoplasm_histologic_grade
# tobacco_smoking_history
# gender.demographic
# race.demographic
# year_of_birth.demographic
# year_of_death.demographic
# age_at_diagnosis.diagnoses
# alcohol_history.exposures
# cigarettes_per_day.exposures
# years_smoked.exposures

## OS.time
head(phenoData$days_to_death.diagnoses)
head(phenoData$days_to_last_follow_up.diagnoses)
phenoData[, 'days_to_death.diagnoses'][is.na(phenoData$days_to_death.diagnoses)] = 0
phenoData[, 'days_to_last_follow_up.diagnoses'][is.na(phenoData$days_to_last_follow_up.diagnoses)] = 0
phenoData$OS.days <- as.numeric(phenoData[, 'days_to_death.diagnoses']) + as.numeric(phenoData[, 'days_to_last_follow_up.diagnoses'])
phenoData$time <- phenoData$OS.days / 30
phenoData <- phenoData[-which(phenoData$time == 0),]
boxplot(phenoData$time)

## status
phenoData$status <- ifelse(phenoData$days_to_death.diagnoses != 0, 0, 1)
table(phenoData$status)

## stage
phenoData$tumor_stage <- ifelse(phenoData$tumor_stage %in% c("stage i", "stage ia", "stage ib"),
                                "early_stage", "later_stage")
table(phenoData$tumor_stage)

## smoke
phenoData$smoked <- ifelse(is.na(phenoData$years_smoked.exposures), 0, phenoData$years_smoked.exposures)
phenoData$smoked_type <- ifelse(phenoData$smoked== 0, "no_smoking", "smoking")
table(phenoData$smoked_type)

## gender
phenoData$gender <- phenoData$gender.demographic
table(phenoData$gender)

library(stringr)
## race
phenoData$race <- str_split(phenoData$race.demographic, ' ', simplify = T)[, 1]
table(phenoData$race)

## bind your clinical data
phenoData <- phenoData[,c("Tumor_Sample_Barcode",
                          "status",
                          "time",
                          "histological_type",
                          "tumor_stage",
                          "smoked_type",
                          "gender",
                          "race")]
colnames(phenoData)[1] <- 'ID'
phenoData$ID <- toupper(phenoData$ID)
dim(phenoData)

## 是否要挑选亚型进行研究
sub_sample <- rownames(phenoData)[phenoData$race == "white"] 
AssayData <- raw_data[, colnames(raw_data) %in% sub_sample]



## 是否要注释基因
dim(AssayData)
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
keytypes(org.Hs.eg.db)
library("stringr")
rownames( AssayData ) <- str_sub(rownames( AssayData ), start = 1, end = 15)
AssayData$ENSEMBL <- rownames( AssayData )
gene_name <- bitr( AssayData$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", 
                   OrgDb = org.Hs.eg.db )
AssayData <- AssayData[gene_name$ENSEMBL, ]
gene_name$max <- apply(AssayData, 1, max)
gene_name <- gene_name[order(gene_name$SYMBOL,
                             gene_name$max,
                             decreasing = T), ]
gene_name <- gene_name[!duplicated(gene_name$SYMBOL), ]
dim( gene_name )




## 是否要挑选lncRNA
{
  gene2type = read.table( './raw_data/gencode.v25lift37.annotation.gtf.gene2type' )
  colnames( gene2type ) = c( "gene", "type" )
  dim( gene2type )
  sort( table( gene2type$type ) )
  gene2type = gene2type[ gene2type[,2] == 'lincRNA', ]
  length( unique( gene2type$gene ) )
  
  gene_name <- gene_name[gene_name$SYMBOL %in% gene2type$gene, ]
  length(gene_name$SYMBOL)
}




AssayData <- AssayData[gene_name$ENSEMBL, ]

rownames(AssayData) <- gene_name$SYMBOL
AssayData <- AssayData[-ncol(AssayData)]

rownames(AssayData) <- sub("\\-", "_", rownames(AssayData))

group_list <- factor(
  ifelse( as.numeric( substr( 
    colnames( AssayData ),14,15)) < 10, 'tumor', 'normal'))

tumor_sample <- colnames(AssayData)[group_list == "tumor"]
normal_sample <- colnames(AssayData)[group_list == "normal"]
paired_sample <- intersect(substring( tumor_sample, 1, 12 ),
                           substring( normal_sample, 1, 12 ))
normal_sample <- normal_sample[substring( normal_sample, 1, 12 ) %in% paired_sample]
AssayData <- AssayData[, c(tumor_sample, normal_sample)]

pheno <- phenoData[phenoData$ID %in% colnames(AssayData), ]
dim(pheno)


group_list <- factor(
  ifelse( as.numeric( substr( 
    colnames( AssayData ),14,15)) < 10, 'tumor', 'normal'))
table(group_list)


## phenodata
pheno$group <- group_list
clinical_info <- pheno
clinical_info$group <- factor(clinical_info$group)

dim(clinical_info)

vars <- colnames(clinical_info)[2:10]

library(tableone)
tableOne <- CreateTableOne(vars, strata = c("group"), data = clinical_info) 
tab_out <- as.data.frame(print(tableOne, md = TRUE))
write.csv(tab_out, "tmp.csv")





## 肿瘤样本和正常样本比较，挑选差异基因
tumor_sample <- tumor_sample[substring( tumor_sample, 1, 12 ) %in% paired_sample]

rep_sample <- tumor_sample[substring(tumor_sample, 1, 12) %in% substring(tumor_sample[substring( tumor_sample, 14, 16 ) %in% "01B"], 1, 12)]

AssayData[1:12, rep_sample]

AssayData <- AssayData[, c(tumor_sample, normal_sample)]

group_list <- factor(
  ifelse( as.numeric( substr( 
    colnames( AssayData ),14,15)) < 10, 'tumor', 'normal'))
table(group_list)

pheno <- pheno[pheno$ID %in% colnames(AssayData), ]
dim(pheno)
dim(AssayData)

draw_heatmap <- function(nrDEG, type){
  library( "pheatmap" )
  nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
  nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
  choose_gene = c( rownames( nrDEG_Z )[1:75], rownames( nrDEG_F )[1:15] )
  choose_matrix = AssayData[ choose_gene, ]
  choose_matrix = t( scale( t( choose_matrix ) ) )
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame( CellType = factor( group_list ))
  rownames( annotation_col ) = colnames( AssayData )
  filename <- paste('./fig/', type, '_heatmap_top100_logFC.png',
                    sep = "", collapse = NULL)
  pheatmap( fontsize = 4, choose_matrix, annotation_col = annotation_col, 
            show_rownames = T, show_colnames = F,
            annotation_legend = T, cluster_cols = F, 
            filename = filename)
}

draw_volcano <- function(nrDEG, type){
  library( "ggplot2" )
  logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
  logFC_cutoff <- 2
  nrDEG$change = as.factor( ifelse( 
    nrDEG$P.Value < 0.01 & abs(nrDEG$logFC) > logFC_cutoff,
    ifelse( nrDEG$logFC > logFC_cutoff, 'UP', 'DOWN' ), 'NOT' ) )
  nrDEGfile <- paste('./data/', type, '_nrDEG_by_logFC.Rdata',
                     sep = "", collapse = NULL)
  save( nrDEG, file = nrDEGfile )
  
  this_tile <- paste0( 
    'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
    '\nThe number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),
    '\nThe number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ) )
  
  volcano = ggplot(data = nrDEG, 
                   aes( x = logFC, y = -log10(P.Value), color = change)) +
    geom_point( alpha = 0.4, size = 1.75 ) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( this_tile ) + 
    theme( plot.title = element_text( size = 15, hjust = 0.5 )) +
    scale_colour_manual( values = c('blue','black','red') )
  print( volcano )
  filename <- paste('./fig/', type, '_volcano_logFC.png',
                    sep = "", collapse = NULL)
  ggsave( volcano, filename = filename )
  dev.off()
}


## DESeq2
library( DESeq2 )
## results of an analysis of differential expression
dds <- DESeqDataSetFromMatrix( countData = AssayData,
                               colData = DataFrame(group_list),
                               design = ~ group_list)

## Differential expression analysis based on the Negative Binomial 
## (a.k.a. Gamma-Poisson) distribution
dds <- DESeq(dds)
resultsNames(dds)
res <- results( dds )
resOrdered <- res[ order(res$padj), ]
head(resOrdered)
nrDEG_DESeq2 <- as.data.frame( resOrdered )
nrDEG_DESeq2 <- na.omit( nrDEG_DESeq2 )
colnames(nrDEG_DESeq2)[2] <- c("logFC") 
colnames(nrDEG_DESeq2)[5] <- c("P.Value") 
draw_heatmap(nrDEG = nrDEG_DESeq2, type = 'DESeq2')
draw_volcano(nrDEG = nrDEG_DESeq2, type = 'DESeq2')


## edgeR
library(edgeR)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Estimate Common Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMCommonDisp(DGElist, design)
## Estimate Trended Dispersion for Negative Binomial GLMs
DGElist <- estimateGLMTrendedDisp(DGElist, design)
## Empirical Bayes Tagwise Dispersions for Negative Binomial GLMs
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

## glmFit fits genewise negative binomial glms, all with the same design matrix 
## but possibly different dispersions, offsets and weights
fit <- glmFit(DGElist, design)
## https://www.biostars.org/p/110861/
## glmLRT conducts likelihood ratio tests for one or more coefficients in the 
## linear model.
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
head(nrDEG_edgeR)
colnames(nrDEG_edgeR)[4] <- c("P.Value") 
draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')


## limma
library(limma)
## A list-based S4 class for storing read counts and associated information 
## from digital gene expression or sequencing technologies.
DGElist <- DGEList( counts = AssayData, group = factor(group_list) )
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
## Calculate Normalization Factors to Align Columns of a Count Matrix
DGElist <- calcNormFactors( DGElist )
DGElist$samples

design <- model.matrix( ~0 + factor(group_list) )
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))

## Transform RNA-Seq Data Ready for Linear Modelling
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
## Fit linear model for each gene given a series of arrays
fit <- lmFit(v, design)

## Construct the contrast matrix corresponding to specified contrasts of a set 
## of parameters.
cont.matrix <- makeContrasts(contrasts = c('tumor-normal'), levels = design)
## Given a linear model fit to microarray data, compute estimated coefficients 
## and standard errors for a given set of contrasts.
fit2 <- contrasts.fit(fit, cont.matrix)
## Empirical Bayes Statistics for Differential Expression
fit2 <- eBayes(fit2)

nrDEG_limma_voom = topTable(fit2, coef = 'tumor-normal', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
draw_heatmap(nrDEG = nrDEG_limma_voom, type = 'limma_voom')
draw_volcano(nrDEG = nrDEG_limma_voom, type = 'limma_voom')


# Step4 Compare three methods ---------------------------------------------

mi <- unique(c(rownames(nrDEG_DESeq2),
               rownames(nrDEG_edgeR),
               rownames(nrDEG_limma_voom)))
lf <- data.frame(DESeq2 = nrDEG_DESeq2[mi, 2],
                 edgeR = nrDEG_edgeR[mi, 1],
                 limma_voom = nrDEG_limma_voom[mi, 1])
cor(na.omit(lf))

library("VennDiagram")

nrDEG_Z = nrDEG_limma_voom[ order( nrDEG_limma_voom$logFC ), ]
nrDEG_F = nrDEG_limma_voom[ order( -nrDEG_limma_voom$logFC ), ]
choose_gene_A = c( rownames( nrDEG_Z )[1:70], rownames( nrDEG_F )[1:30] )

nrDEG_Z = nrDEG_edgeR[ order( nrDEG_edgeR$logFC ), ]
nrDEG_F = nrDEG_edgeR[ order( -nrDEG_edgeR$logFC ), ]
choose_gene_B = c( rownames( nrDEG_Z )[1:70], rownames( nrDEG_F )[1:30] )

nrDEG_Z = nrDEG_DESeq2[ order( nrDEG_DESeq2$logFC ), ]
nrDEG_F = nrDEG_DESeq2[ order( -nrDEG_DESeq2$logFC ), ]
choose_gene_C = c( rownames( nrDEG_Z )[1:70], rownames( nrDEG_F )[1:30] )

## Venn Diagram
venn.plot <- venn.diagram(x = list(A = choose_gene_A, B = choose_gene_B, C = choose_gene_C), 
                          filename = "DIFF.png", height = 450, width = 450,
                          resolution = 300, imagetype = "png", col = "transparent",
                          fill = c("cornflowerblue", "darkorchid1", "red"),
                          alpha = 0.50, cex = 0.45, cat.cex = 0.45)

choose_gene <- intersect(choose_gene_A, choose_gene_B)
choose_gene <- intersect(choose_gene, choose_gene_C)


cbind(nrDEG_limma_voom[choose_gene, c(1,4)],
      nrDEG_edgeR[choose_gene, c(1,4)],
      nrDEG_DESeq2[choose_gene, c(2,5)])

choose_matrix <- AssayData[choose_gene, ]
choose_matrix <- log10(choose_matrix + 0.01)
library( "pheatmap" )

choose_matrix <- t( scale( t( choose_matrix ) ) )
choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2

annotation_col = data.frame( CellType = factor( group_list ))
rownames( annotation_col ) = colnames( AssayData )

pheatmap( fontsize = 8, choose_matrix, annotation_col = annotation_col, 
          show_rownames = T, show_colnames = F,
          annotation_legend = T, cluster_cols = F)














## 是否只保留肿瘤样本
AssayData <- AssayData[, group_list == 'tumor']

pheno <- phenoData[phenoData$ID %in% colnames(AssayData), ]
dim(pheno)
dim(AssayData)




t_exp <- AssayData[apply(AssayData, 1, function(x) sum(x == 0)) < (ncol(AssayData)*0.49), ]
dim(t_exp)

## single
library(survival)
group_data <- apply(t_exp , 1 , function(gene){
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  names(group) <- name
  return(group)
})
group_data <- as.data.frame(group_data, stringsAsFactors = F)
survival_dat <- data.frame(gender = pheno$gender,
                           status = pheno$status,
                           time = pheno$time,
                           smoked_type = pheno$smoked_type,
                           stringsAsFactors = F)
survival_dat <- cbind(group_data, survival_dat)
dim(survival_dat)
covariates <- as.character(colnames(survival_dat))
univ_formulas <- sapply(covariates,
                        function(x){
                          ##print(x)
                          as.formula(paste('Surv(time, status)~', x))
                        })
univ_formulas$status <- NULL
univ_formulas$time <- NULL
length(univ_formulas)
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survival_dat)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 2)
                         beta <- signif(x$coef[1], digits = 2)
                         HR <- signif(x$coef[2], digits = 2)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                         HR.length <- length(x$conf.int[,"lower .95"])
                         if (HR.length > 1) {
                           res <- lapply(1:HR.length, function(x){
                             HR <- paste0(HR, " (", 
                                          HR.confint.lower[x], "-", HR.confint.upper[x], ")")
                             s.res <- c(beta, HR, p.value)
                             return(s.res)
                           })
                           res <- as.data.frame(do.call(cbind, res))
                           names(res) <- as.character(sub('race', "", rownames(x$conf.int)))
                           return(res)
                         }else{
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res <- c(beta, HR, p.value)
                           names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                           return(res)
                         }
                       })
res_single <- as.data.frame(t(do.call(cbind, univ_results)))
table(res_single$p.value <= 0.01)
res_single <- res_single[res_single$p.value <= 0.01, ]
res_single <- res_single[order(res_single$p.value), ]
single_pick <- rownames(res_single)
res_single[single_pick,]

pick_gene <- single_pick

library(formattable)
sign_formatter <- formatter("span", 
                            style = x ~ style(color = ifelse(x > 0, "lightpink", 
                                                             ifelse(x < 0, "lightblue", "black"))))
sign_formatter(c(-1, 0, 1))
formattable( res_single[pick_gene,],
             list(coef = sign_formatter,
                  p.value = color_bar("lightblue")))

pick_gene

## multi
multi_results <- apply(t_exp[pick_gene, ] , 1 , function(gene){
  ## gene <- t_exp[1, ]
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             gender = pheno$gender,
                             status = pheno$status,
                             time = pheno$time,
                             smoked_type = pheno$smoked_type,
                             stringsAsFactors = F)
  res.cox <- coxph(Surv(time, status) ~ smoked_type + gender + group, 
                   data =  survival_dat)
  ## summary(res.cox)
  x <- summary(res.cox)
  p.value <- signif(x$coef['grouplow',5], digits = 2)
  beta <- signif(x$coef['grouplow',1], digits = 2)
  HR <- signif(x$coef['grouplow',2], digits = 2)
  HR.confint.lower <- signif(x$conf.int['grouplow',"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int['grouplow',"upper .95"], 2)
  HR <- paste0(HR, " (", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, p.value)
  names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
  return(res)
})
multi_results <- as.data.frame(t(as.data.frame(multi_results)))
res_multi <- multi_results[order(multi_results$p.value), ]
pick_gene <- rownames(res_multi)
res_multi[pick_gene,]
formattable( res_multi[pick_gene,],
             list(coef = sign_formatter,
                  p.value = color_bar("lightblue")))

## OS paint
library(survminer)
model_exp <- t_exp[pick_gene, ]
dim(model_exp)
x <- matrix(nrow = nrow(model_exp), ncol = 1)
for (i in 1:nrow(model_exp)) {
  gene <- model_exp[i, ]
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
  x[i] <- surv_pvalue(fit)[4]
}

p.v <- as.data.frame(do.call(rbind, x))
rownames(p.v) <- rownames(model_exp)
library(stringr)
p.v$p.value <- str_split(p.v$V1, ' ', simplify = T)[, 3]
p.v
pick_gene <- rownames(p.v)[p.v$p.value < 0.01]
pick_gene

model_exp <- t_exp[pick_gene, ]
os.paint <- apply(model_exp, 1, function(x){
  gene <- unlist(x)
  group <- ifelse(gene >= median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
  ggsurvplot(fit, data = survival_dat,
             ##surv.median.line = "hv",
             legend.title = "Group",
             legend.labs = c("High", "Low"),
             pval = TRUE,
             ##conf.int = TRUE,
             palette = "jco",
             ggtheme = theme_bw(),
             pval.size = 4
  )
})
arrange_ggsurvplots(os.paint, print = TRUE, ncol = 2, nrow = 2)



# step ROC ----------------------------------------------------------------

exp <- AssayData

exp <- exp[apply(exp, 1, function(x) sum(x == 0)) < (ncol(exp)*0.49), ]
pheno <- phenoData[colnames(exp), ]


model_exp <- t(log10(exp[pick_gene,] + 1))
colnames(model_exp) <- pick_gene
dat <- cbind(pheno, model_exp)

## ROC
library("survminer")
colnames(dat)
s <- Surv(time, status) ~ LINC00102 + LINC00524

cox_model <- coxph(s, data = dat )
summary(cox_model, data = dat)
options(scipen = 1)
ggforest(cox_model, data = dat, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)

cox.prob <- predict(cox_model)
library(Hmisc)
with(dat, rcorr.cens(cox.prob, Surv(time, status)))
## 1-(C Index)
library(timeROC)
pheno$fp <- cox.prob
with(pheno,
     ROC <<- timeROC(T = time,
                     delta = status,
                     marker = fp,
                     cause = 1,
                     weighting = "marginal",
                     times = c(30, 60),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC, time = 60, col = "red", add = F)
plot(ROC, time = 30, col = "blue", add = T)

ROC$AUC

legend("top",
       legend = c("t=30, AUC=67.7","t=60, AUC=76.9"),
       ncol = 4,
       cex = 0.8,
       bty = "n",
       col = c("blue","red"),
       lty = 1,lwd = 2) 


confint(ROC)


fp <- predict(cox_model, dat, type = "risk");boxplot(fp)
fp <- predict(cox_model, dat, type = "expected");boxplot(fp)
plot(fp, pheno$days)
fp <- predict(cox_model, dat);boxplot(fp)
basehaz(cox_model) 
library(Hmisc)
options(scipen = 200)
with(dat, rcorr.cens(fp, Surv(time, status)))

library(cowplot)
library(pheatmap)
fp_dat <- data.frame(sample = 1:length(fp), predictor_score = as.numeric(sort(fp)))
sur_dat <- data.frame(sample = 1:length(fp),
                      survival_time = pheno[names(sort(fp)), 'time'],
                      status = pheno[names(sort(fp)), 'status']  ) 
sur_dat$status <- ifelse(sur_dat$status == 0, 'alive', 'death')
colnames(dat)
exp_dat <- dat[names(sort(fp)), 9:12]

plot.point <- ggplot(fp_dat, 
                     aes(x = sample, y = predictor_score)) + 
  geom_point() + xlab("")
print(plot.point)

plot.sur <- ggplot(sur_dat, aes(x = sample, y = survival_time)) + 
  geom_point(aes(col = status)) + xlab("")
print(plot.sur)

mycolors <- colorRampPalette(c("black", "green", "red"), bias = 1.2)(100)
exp_dat <- t(scale(exp_dat))
exp_dat[exp_dat > 1] = 1
exp_dat[exp_dat < -1] = -1
plot.h <- pheatmap(exp_dat, col = mycolors, show_colnames = F, cluster_cols = T)
plot.h <- pheatmap(exp_dat, col = mycolors, show_colnames = F, cluster_cols = F)
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B","C"),
          align = 'v',ncol = 1)



fp <- predict(cox_model, dat);boxplot(fp)
median(fp)
group_list <- ifelse( fp >= median(fp), 'high', 'low' )
table(group_list)
pheno <- pheno[names(group_list), ]
survival_dat <- data.frame(group = group_list,
                           status = pheno$status,
                           time = pheno$time,
                           stringsAsFactors = F)
median(survival_dat[survival_dat$group == 'high', 3])
median(survival_dat[survival_dat$group == 'low', 3])
fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
ggsurvplot(fit, data = survival_dat,
           ##surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("High", "Low"),
           pval = TRUE,
           ##conf.int = TRUE,
           palette = "jco",
           ggtheme = theme_bw())


## gender

exp <- AssayData
exp <- exp[apply(exp, 1, function(x) sum(x == 0)) < (ncol(exp)*0.49), ]
pheno <- phenoData[colnames(exp), ]

table(pheno$smoked_type)

e_phe <- pheno[pheno$smoked_type == "no_smoking", ]
l_phe <- pheno[pheno$smoked_type == "smoking", ]

table(pheno$gender)

l_phe <- pheno[pheno$gender == "female", ]
l_phe <- pheno[pheno$gender == "male", ]

e_exp <- AssayData[, rownames(e_phe)]
l_exp <- AssayData[, rownames(l_phe)]

exp <- e_exp
exp <- l_exp

exp <- exp[apply(exp, 1, function(x) sum(x == 0)) < (ncol(exp)*0.49), ]
pheno <- phenoData[colnames(exp), ]

model_exp <- t(log2(exp[pick_gene,] + 1))
colnames(model_exp) <- pick_gene
dat <- cbind(pheno, model_exp)

dat$gender <- factor(dat$gender)
dat$smoked_type <- factor(dat$smoked_type)

library("survminer")
colnames(dat)
s <- Surv(time, status) ~ LINC00473 + LINC00102 + LINC01290 + LINC00524

cox_model <- coxph(s, data = dat )
summary(cox_model, data = dat)

cox.prob <- predict(cox_model)
library(timeROC)
pheno$fp <- cox.prob
with(pheno,
     ROC <<- timeROC(T = time,
                     delta = status,
                     marker = fp,
                     cause = 1,
                     weighting = "marginal",
                     times = c(48, 50),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC, time = 50, col = "blue", add = F)
plot(ROC, time = 60, col = "red", add = T)
ROC$AUC

legend("topleft",
       legend = c("Shorter_smoking_time\nAUC=74.7","Female\nAUC=92.1"),
       ncol = 4,
       cex = 0.8,
       bty = "n",
       col = c("blue","red"),
       lty = 1,
       lwd = 2) 


confint(ROC)
