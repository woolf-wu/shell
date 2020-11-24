## Author: Wajid Jawaid
## Date: 1 August 2019
## enrichR package: sends data to http://http://amp.pharm.mssm.edu/Enrichr/ for gene enrichment
## in multiple databases.

rm(list = objects( all = TRUE ))

options(stringsAsFactors = FALSE)

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()

library("httr")

getEnrichr <- function(url, ...) {
  tryCatch(
    {
      options(enrichR.live = TRUE)
      x <- GET(url = url, ...)
    },
    error = function(err_msg) {
      message("EnrichR website not responding")
      options(enrichR.live = FALSE)
    },
    finally = function() {
      invisible(x) 
    }
  )
}

options(enrichR.base.address = "http://amp.pharm.mssm.edu/Enrichr/")
options(enrichR.live = TRUE)
packageStartupMessage("Welcome to enrichR\nChecking connection ... ", appendLF = FALSE)
getEnrichr(url = paste0(getOption("enrichR.base.address"), "datasetStatistics"))
if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")



enrichr <- function(genes, dbs) {
  cat("Uploading data to Enrichr... ")
  if (is.vector(genes)) {
    temp <- POST(url = paste0(getOption("enrichR.base.address"), "enrich"),
                 body = list(list = paste(genes, collapse = "\n")))
  } else if (is.data.frame(genes)) {
    temp <- POST(url = paste0(getOption("enrichR.base.address"), "enrich"),
                 body = list(list = paste(paste(genes[,1], genes[,2], sep = ","),
                                      collapse = "\n")))
  }
  getEnrichr(url = paste0(getOption("enrichR.base.address"), "share"))
  cat("Done.\n")
  
  result <- lapply(dbs, function(x) {
    cat("  Querying ", x, "... ", sep = "")
    r <- getEnrichr(url = paste0(getOption("enrichR.base.address"), "export"),
                    query = list(file = "API", backgroundType = x))
    if (!getOption("enrichR.live")) return()
    r <- gsub("&#39;", "'", intToUtf8(r$content))
    tc <- textConnection(r)
    r <- read.table(tc, sep = "\t", header = TRUE, quote = "", comment.char = "")
    close(tc)
    cat("Done.\n")
    return(r)
  })
  
  cat("Parsing results... ")
  names(result) <- dbs
  cat("Done.\n")
  return(result)
}

gene_list <- read.table("./load/WGCNA_20190718_CytoscapeInput.node.anno.sort (2).txt")
gene_list <- gene_list[gene_list$V4 == "mRNA",]

color_list <- unique(gene_list[,2])
head(color_list)
length(color_list)

for (col in color_list) {
  print(which(col %in% color_list))
  print(col)
  pick_gene <- gene_list[gene_list[,2] == col,]
  if (nrow(pick_gene) > 1) {
    dbs <- as.list(c("BioCarta_2016", "GO_Biological_Process_2018", "KEGG_2019_Human"))
    enriched <- enrichr(pick_gene[,3], dbs)
    for (id in 1:length(dbs)) {
      print(id)
      print(dbs[id])
      file_name <- paste(col, dbs[id], "r.csv", sep = "_")
      print(file_name)
      data <- enriched[id]
      write.csv(data, file = file_name)
    }
  }
}
