# Step0 Before starting your project --------------------------------------

dep=c("ggpubr", "magick", "picante", "rgl", "vegan")
inp(dep)


rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()

Sys.setlocale("LC_ALL","English")



# Step1 Starting your project -------------------------------------------

data_raw <- read.csv(file = "./PGA.xls",
           header = T, sep = '\t',
           #row.names = 1,
           stringsAsFactors = FALSE,
           comment.char = "#")
library("ggpubr")

type_x = "sample_Type2"
type_y = "CCF.adj"
my_comparisons <- list( c("NF", "NI"), c("PF", "NF"), c("PI", "NI"), c("PI", "PF") )
ggboxplot(data_raw, x = type_x, y = type_y,
          order = c("NI","NF","PI","PF"),
          color = type_x, palette = c("#00AFBB", "#E7B800", "#FC4E07", "#c00000"),
          add = "none", 
          shape = type_x
          ) + 
  stat_compare_means(comparisons = my_comparisons) 

type_x = "sample_Type1"
type_y = "CCF.adj"
my_comparisons <- list( c("P", "N") )
ggboxplot(data_raw, x = type_x, y = type_y,
          order = c("N","P"),
          color = type_x, palette = c("#00AFBB", "#E7B800"),
          add = "none", 
          shape = type_x
) + 
  stat_compare_means(comparisons = my_comparisons) 



gene_list=c("NOTCH1", "TP53", "TMEM132D", "ANAPC1", "FAT1", "NOTCH2", "LRP1B", "CHEK2", "CDKN2A")

gene="CDKN2A"
for (gene in gene_list) {
  data_raw_1 = data_raw[data_raw$gene_name == gene, ]
  type_x = "Type"
  type_y = "CCF.adj"
  my_comparisons <- list( c("PI", "NI") )
  p <- ggboxplot(data_raw_1, x = type_x, y = type_y,
                 order = c("NI","PI"),
                 color = type_x, palette = c("#00AFBB", "#E7B800"),
                 add = "none", 
                 shape = type_x,
                 title = gene
  ) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",) 
  
  file=paste(gene, '.pdf', sep = "", collapse = NULL)
  ggsave( filename = file, plot = p)
}




x = data_raw[data_raw$type == "N", ]
colnames(x) = c("NF", "NI")

ggplot(data=x, aes(x=NI, y=NF))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=x, method = "pearson") + coord_fixed() +
  scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 5))


x = data_raw[data_raw$type == "P", ]
colnames(x) = c("PF", "PI")

ggplot(data=x, aes(x=PI, y=PF))+geom_point(color="red")+stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=x, method = "pearson") + 
  coord_fixed() + 
  scale_x_continuous(breaks=seq(0, 16, 4), limits=c(0, 16))






data_raw <- read.csv(file = "input/LSK_KO_ch17.depth_stat.xls",
                     header = F, sep = '\t',
                     #row.names = 1,
                     stringsAsFactors = FALSE,
                     comment.char = "#")
library("ggpubr")
library("RColorBrewer")
options(scipen=200)
ggplot(data_raw, aes(x=V2, y=V3), xaxs="d") + 
  geom_bar(stat="identity", color = "#5a5a5a") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(breaks=seq(0, 600000, 100000), limits=c(0, 600000), expand = c(0,0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(rep(2,4),'lines')) +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 0.5) +
  geom_vline(xintercept = 46028374) +
  #geom_vline(xintercept = 46025486) +
  annotate(geom = "rect", xmin = 46025347, xmax = 46025543, ymin=0, ymax=Inf, fill='black', alpha = .08) +
  annotate(geom = "segment", x = 46025486, xend = 46025486, y=530000, yend = 0, color='black', size = 0.6)

46022734-46028374

patient <- read.csv(file = "input/patient_54.list",
                     header = F, sep = '\t',
                     #row.names = 1,
                     stringsAsFactors = FALSE,
                     comment.char = "#")

data_raw <- read.csv(file = "input/sample_108-SBS.xls",
                     header = T, sep = '\t',
                     #row.names = 1,
                     stringsAsFactors = FALSE,
                     comment.char = "#")
library("ggpubr")

group_sample <- c("W00139", "W00204", "W00746", "W00757", "W00790", "W00829", "W00947", "W00980", 
                  "W01075", "W02027", "W02141", "W02240", "W02438", "W02566", "W03075", "W03478", 
                  "W03637", "W05804", "W05853", "W06456", "W08415", "W08632", "W08711", "W10475", 
                  "W10585", "W10793", "W11008", "W11056", "W11320", "W12113", "W12349", "W12412", 
                  "W12452", "W12833", "W13838", "W13900", "W14414", "W14458", "W14494", "W14523", 
                  "W14539", "W15426", "W15603", "W15837", "W15840", "W15988", "W16593", "W16995", 
                  "W17409", "W17947", "W18221", "W20095", "W22379", "W22410")

for (id in group_sample) {
  sample_x = paste0(id, "_1", collapse = "")
  sample_y = paste0(id, "_2", collapse = "")
  
  ditf <- data_raw[data_raw$Sample==sample_x, c("chr", "CCF.adj", "Clonality")]
  colnames(ditf) <- c("chr", "t1.CCF.adj", "t1.Clonality")
  dtb <- data_raw[data_raw$Sample==sample_y,  c("chr", "CCF.adj", "Clonality")]
  colnames(dtb) <- c("chr", "t2.CCF.adj", "t2.Clonality")
  
  x <- merge(ditf, dtb, all = T)
  x$Clonality <- x$t1.Clonality
  x[is.na(x)]<-0
  x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "Clonality"] <- x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "t1.Clonality"]
  x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "Clonality"] <- x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "t2.Clonality"]
  x[x$t1.Clonality=="clonal" & x$t2.Clonality=="subclonal", "Clonality"] <- "t1_clonal_t2_subclonal"
  x[x$t1.Clonality=="subclonal" & x$t2.Clonality=="clonal", "Clonality"] <- "t2_clonal_t1_subclonal"
  
  title_name = patient[patient$V1 == id, 2]
  ggplot(data=x, aes(x=t1.CCF.adj, y=t2.CCF.adj)) + 
    geom_point(aes(color = factor(Clonality))) + 
    coord_fixed() +
    scale_color_manual(values=c('#016eb7','#bf0e10', '#e8ae4a', '#8969aa')) +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
    scale_x_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
    theme(panel.grid.major=element_line(colour=NA),
          legend.title = element_blank()) +
    ggtitle(title_name)
    
  file_path <- paste0(title_name, "_", id, "_ccf.png", collapse = "")
  ggsave( file = file_path , width = 6, height = 6, type = "cairo", dpi = 600) 
}



data_raw <- read.csv(file = "input/sample_12-CCF.xls",
                     header = T, sep = '\t',
                     #row.names = 1,
                     stringsAsFactors = FALSE,
                     comment.char = "#")
library("ggpubr")
library("ggrepel")

group_sample <- c("P1209625", "P1213701", "P1214387", "P1216086", "P1217935", "P53211")

new_x = data.frame()
for (id in group_sample) {
  sample_x = paste0(id, "DITF", collapse = "")
  sample_y = paste0(id, "DTB", collapse = "")
  
  ditf <- data_raw[data_raw$Sample==sample_x, c("chr", "CCF.adj", "Clonality", "Drivrtgene")]
  colnames(ditf) <- c("chr", "DITF.CCF.adj", "DITF.Clonality", "Drivrtgene")
  dtb <- data_raw[data_raw$Sample==sample_y,  c("chr", "CCF.adj", "Clonality", "Drivrtgene")]
  colnames(dtb) <- c("chr", "DTB.CCF.adj", "DTB.Clonality", "Drivrtgene")
  
  x <- merge(ditf, dtb, all = T)
  x$Clonality <- x$DITF.Clonality
  x[is.na(x)]<-0
  x[(x$DITF.Clonality!="0" & x$DTB.Clonality=="0"), "Clonality"] <- x[(x$DITF.Clonality!="0" & x$DTB.Clonality=="0"), "DITF.Clonality"]
  x[(x$DTB.Clonality!="0" & x$DITF.Clonality=="0"), "Clonality"] <- x[(x$DTB.Clonality!="0" & x$DITF.Clonality=="0"), "DTB.Clonality"]
  x[x$DITF.Clonality=="clonal" & x$DTB.Clonality=="subclonal", "Clonality"] <- "DITF_clonal_DTB_subclonal"
  x[x$DITF.Clonality=="subclonal" & x$DTB.Clonality=="clonal", "Clonality"] <- "DTB_clonal_DITF_subclonal"
  
  xx <- x[x$DTB.Clonality == 0,]
  new_x = rbind(xx, new_x)
  file <- paste0(id, ".xls", collapse = "")
  write.table(x, file = file, quote = F, col.names = T, sep = "\t", row.names = F)
  
  ggplot(data=x, aes(x=DITF.CCF.adj, y=DTB.CCF.adj)) + 
    geom_point(aes(color = factor(Clonality))) + 
    geom_text_repel(aes(label=Drivrtgene), size=3, arrow = arrow(length = unit(0.001, "npc")),
                    box.padding = 1.2) +
    coord_fixed() +
    scale_color_manual(values=c('#016eb7', '#e8ae4a', '#8969aa', '#bf0e10')) +
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
    scale_x_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
    theme(panel.grid.major=element_line(colour=NA),
          legend.title = element_blank()) +
    ggtitle(id)
  
  file_path <- paste0(id, "_ccf.pdf", collapse = "")
  ggsave( file = file_path , width = 9.3, height = 6.8) 
  #ggsave( file = file_path , width = 9.3, height = 6.8, type = "cairo", dpi = 600) 
}

write.table(new_x, file = "tmp.xls", quote = F, col.names = T, sep = "\t", row.names = F)



new_x = data.frame()

group_sample = c("W10793", "W00204", "W00790", "W02027", "W17409", "W14523", "W11056", "W10475", "W14414", "W12113", "W13900", "W00139", "W03075", "W14539", "W15603", "W02240", "W02141", "W12412", "W11320", "W08415", "W15426", "W00746", "W02438", "W06456", "W03637", "W08632", "W00829")
group_sample = c("W00757", "W05853", "W12349", "W14458", "W05804", "W14494", "W00947", "W01075", "W08711", "W20095", "W18221", "W13838", "W11008", "W16593", "W22379", "W22410", "W12833", "W16995", "W17947", "W15837", "W00980", "W15840", "W15988", "W03478", "W10585", "W02566", "W12452")

for (id in group_sample) {
  sample_x = paste0(id, "_1", collapse = "")
  sample_y = paste0(id, "_2", collapse = "")
  time = data_raw[data_raw$Sample==sample_x,"Time"][1]
  
  ditf <- data_raw[data_raw$Sample==sample_x, c("chr", "CCF.adj", "Clonality")]
  colnames(ditf) <- c("chr", "t1.CCF.adj", "t1.Clonality")
  dtb <- data_raw[data_raw$Sample==sample_y,  c("chr", "CCF.adj", "Clonality")]
  colnames(dtb) <- c("chr", "t2.CCF.adj", "t2.Clonality")
  
  x <- merge(ditf, dtb, all = T)
  
  x$sample_name = id
  
  x$time = "0"
  
  x[x$time=="0","time"] = time
  
  x$Clonality <- x$t1.Clonality
  x[is.na(x)]<-0
  x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "Clonality"] <- x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "t1.Clonality"]
  x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "Clonality"] <- x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "t2.Clonality"]
  x[x$t1.Clonality=="clonal" & x$t2.Clonality=="subclonal", "Clonality"] <- "t1_clonal_t2_subclonal"
  x[x$t1.Clonality=="subclonal" & x$t2.Clonality=="clonal", "Clonality"] <- "t2_clonal_t1_subclonal"

  new_x = rbind(x, new_x)  
}


library("rgl")
library("caTools")
library("magick")
y=new_x[new_x$t1.CCF.adj!=0,]
new_x=y[y$t2.Clonality!=0,]
t1=new_x$t1.CCF.adj
t2=new_x$t2.CCF.adj
t3=new_x$time

cl<-c(ifelse(new_x$Clonality=="clonal",'#016eb7',
             ifelse(new_x$Clonality=="subclonal",'#bf0e10',
                    ifelse(new_x$Clonality=="t1_clonal_t2_subclonal",'#e8ae4a','#8969aa'))))


open3d(windowRect = c(0, 0, 700, 700))
p1 = plot3d(t2,t1,t3, col = cl)
outputFile = "ccf_3d_N"
movie3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 12, 
        dir = getwd(), movie = outputFile)


new_x = data.frame()

group_sample <- c("W00139", "W00204", "W00746", "W00757", "W00790", "W00829", "W00947", "W00980", 
                  "W01075", "W02027", "W02141", "W02240", "W02438", "W02566", "W03075", "W03478", 
                  "W03637", "W05804", "W05853", "W06456", "W08415", "W08632", "W08711", "W10475", 
                  "W10585", "W10793", "W11008", "W11056", "W11320", "W12113", "W12349", "W12412", 
                  "W12833", "W13838", "W13900", "W14414", "W14458", "W14494", "W14523", 
                  "W14539", "W15426", "W15603", "W15837", "W15840", "W15988", "W16593", "W16995", 
                  "W17409", "W17947", "W18221", "W20095", "W22379", "W22410")

group_sample = c("W10793", "W00204", "W00790", "W02027", "W17409", "W14523", "W11056", "W10475", "W14414", "W12113", "W13900", "W00139", "W03075", "W14539", "W15603", "W02240", "W02141", "W12412", "W11320", "W08415", "W15426", "W00746", "W02438", "W06456", "W03637", "W08632", "W00829")
group_sample = c("W00757", "W05853", "W12349", "W14458", "W05804", "W14494", "W00947", "W01075", "W08711", "W20095", "W18221", "W13838", "W11008", "W16593", "W22379", "W22410", "W12833", "W16995", "W17947", "W15837", "W00980", "W15840", "W15988", "W03478", "W10585", "W02566", "W12452")

for (id in group_sample) {
  sample_x = paste0(id, "_1", collapse = "")
  
  ditf <- data_raw[data_raw$Sample==sample_x, c("chr", "CCF.adj", "Clonality")]
  colnames(ditf) <- c("chr", "t1.CCF.adj", "t1.Clonality")
  
  x = ditf
  
  x$sample_name = id
  
  x$Clonality <- x$t1.Clonality
  x[is.na(x)]<-0
  x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "Clonality"] <- x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "t1.Clonality"]
  x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "Clonality"] <- x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "t2.Clonality"]
  x[x$t1.Clonality=="clonal" & x$t2.Clonality=="subclonal", "Clonality"] <- "t1_clonal_t2_subclonal"
  x[x$t1.Clonality=="subclonal" & x$t2.Clonality=="clonal", "Clonality"] <- "t2_clonal_t1_subclonal"
  
  new_x = rbind(x, new_x)  
}

for (id in group_sample) {

  sample_y = paste0(id, "_2", collapse = "")

  dtb <- data_raw[data_raw$Sample==sample_y,  c("chr", "CCF.adj", "Clonality")]
  colnames(dtb) <- c("chr", "t2.CCF.adj", "t2.Clonality")
  
  x = dtb
  
  x$sample_name = id
  
  x$Clonality <- x$t2.Clonality
  x[is.na(x)]<-0
  x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "Clonality"] <- x[(x$t1.Clonality!="0" & x$t2.Clonality=="0"), "t1.Clonality"]
  x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "Clonality"] <- x[(x$t2.Clonality!="0" & x$t1.Clonality=="0"), "t2.Clonality"]
  x[x$t1.Clonality=="clonal" & x$t2.Clonality=="subclonal", "Clonality"] <- "t1_clonal_t2_subclonal"
  x[x$t1.Clonality=="subclonal" & x$t2.Clonality=="clonal", "Clonality"] <- "t2_clonal_t1_subclonal"
  
  new_x = rbind(x, new_x)  
}

type = c("clonal", "subclonal", "t1_clonal_t2_subclonal","t2_clonal_t1_subclonal")
new_y = data.frame(matrix(type, 4, 1), row.names = 1 )
for (i in unique(new_x$sample_name)) {
  y <- as.data.frame(table(new_x[new_x$sample_name==i,"Clonality"]), row.names = 1)
  colnames(y) <- i
  for (j in type) {
    if (is.na(y[j, 1])) {
      y[j, 1] = "0"
    }else{
      y[j, 1] = sum(new_x[(new_x$sample_name==i) & (new_x$Clonality==j),
                          "t1.CCF.adj"])
    }
  }
  new_y <- cbind(new_y, y)
}


type = c("clonal", "subclonal", "t1_clonal_t2_subclonal","t2_clonal_t1_subclonal")
new_y = data.frame(matrix(type, 4, 1), row.names = 1 )
for (i in unique(new_x$sample_name)) {
  y <- as.data.frame(table(new_x[new_x$sample_name==i,"Clonality"]), row.names = 1)
  colnames(y) <- i
  for (j in type) {
    if (is.na(y[j, 1])) {
      y[j, 1] = "0"
    }
  }
  new_y <- cbind(new_y, y)
}

library("vegan")
library("picante")

otu <- t(new_y)
otu=as.data.frame(lapply(otu,as.numeric))
##物种丰富度 Richness 指数
richness <- rowSums(otu > 0); richness
#或
richness <- estimateR(otu)[1, ]; richness

##Shannon（以下 Shannon 公式的对数底数均设使用 e，在 R 中即表示为 exp(1)）
#Shannon 指数
shannon_index <- diversity(otu, index = 'shannon', base = exp(1)); shannon_index

#Shannon 多样性
shannon_diversity <- exp(1)^shannon_index; shannon_diversity

#Shannon 均匀度（Pielou 均匀度）
pielou <- shannon_index / log(richness, exp(1)); pielou

##Simpson
#Gini-Simpson 指数（我们平时常用的 Simpson 指数即为 Gini-Simpson 指数）
gini_simpson_index <- diversity(otu, index = 'simpson'); gini_simpson_index

#经典 Simpson 指数（使用频率比较低）
simpson_index <- 1 - gini_simpson_index; simpson_index

#Invsimpson 指数（Gini-Simpson 的倒数）
invsimpson_index <- 1 / gini_simpson_index; invsimpson_index
#或
invsimpson_index <- diversity(otu, index = 'invsimpson'); invsimpson_index

#Simpson 多样性
simpson_diversity <- 1 / (1 - gini_simpson_index); simpson_diversity

#Simpson 均匀度（equitability 均匀度）
equitability <- 1 / (richness * (1 - gini_simpson_index)); equitability

##Chao1 & ACE
#Chao1 指数
chao1 <- estimateR(otu)[2, ]; chao1

#ACE 指数
ace <- estimateR(otu)[4, ]; ace

##goods_coverage 指数
goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu); goods_coverage




data_raw$sum=data_raw$SBS1+data_raw$SBS5+data_raw$SBS6

data_raw$SBS1 = data_raw$SBS1/data_raw$sum
data_raw$SBS5 = data_raw$SBS5/data_raw$sum
data_raw$SBS6 = data_raw$SBS6/data_raw$sum

clonsigp <- data_raw[data_raw$clone_type=="clonal", ]
clonsigp <- clonsigp[clonsigp$time=="1", ]
clonsigp <- data_raw[data_raw$clone_type=="clonal", ]
subsigp <- clonsigp[clonsigp$time=="2", ]

subsigp <- data_raw[data_raw$clone_type=="subclonal", ]
clonsigp <- subsigp[subsigp$time=="1", ]
subsigp <- data_raw[data_raw$clone_type=="subclonal", ]
subsigp <- subsigp[subsigp$time=="2", ]

labs <- c("first", "second")
pdf(file.path("Clonal_vs_sublclonal_signature_proportions.pdf"), 
    width = 10, height = 3.5)
bp <- plot(-10, xlim = c(0, 5 * 2 + 4), ylim = c(-0.1, 1.2), axes = F, bty = "n", xlab = "", ylab = "Proportion of mutations", 
           main = "Clonal vs.Subclonal Signature Exposures")
for (j in 1:3) {
  sig <- colnames(clonsigp)[j+3]
  ind <- which(clonsigp[, sig] > 0 | subsigp[, sig] > 0)
  tmp <- data.frame(clon = clonsigp[ind, sig], sub <- subsigp[ind, sig])
  for (i in 1:nrow(tmp)) {
    points(((j - 1) * 2.5):((j - 1) * 2.5 + 1), tmp[i, ], type = "b", pch = c(19), col = "red")
    axis(1, at = ((j - 1) * 2.5):((j - 1) * 2.5 + 1), 
         labels = labs, pos = 0, las = 2, cex.axis = 1.5, 
         tick = FALSE, cex.axis = 1)
  }
  text(mean(((j - 1) * 2.5):((j - 1) * 2.5 + 1)), 1.2, 
       sub("Signature.", "", sig), col = "red")
}
axis(2, at = c(0, 0.5, 1), las = 2)
abline(h = 0)
dev.off()
