rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()


# installPackage("IRanges")
# installPackage("limma")
# installPackage("gridBase")
# 
# library(devtools)
# install_github('hdng/clonevol')
# install_github("genome/bmm")
# install_github("genome/sciClone")
# install_github("chrisamiller/fishplot")
# installPackage("packcircles")
# installPackage("igraph")
library(sciClone)
library(clonevol)
library(fishplot)

sample_list <- c("P1209625", "P1213701", "P1216086", "P1217935", "P53211")

sample_list <- c("glmBM")


for (sample in sample_list) {

#sample = "W00757"

sample_1 <- paste( sample, '_2_1', sep = '' )
sample_2 <- paste( sample, '_3_1', sep = '' )


Rdata_file_1 <- paste( './input/', sample_1, '_sciClone.Rdata', sep = '' )
Rdata_file_2 <- paste( './input/', sample_2, '_sciClone.Rdata', sep = '' )

load(file = Rdata_file_1)
v1 = vaf1
c1 = cnv1
a1 = anno1

load(file = Rdata_file_2)
v2 = vaf1
c2 = cnv1
a2 = anno1


a=merge(v1, v2, by=c("chr", "pos"), all=TRUE, sort=TRUE)
a[is.na(a)] = 0
v1=a[,c(1,2,3,4,5)]
v2=a[,c(1,2,6,7,8)]

colnames(v1) <- c("chr", "pos", "ref_reads", "var_reads", "vaf")
colnames(v2) <- c("chr", "pos", "ref_reads", "var_reads", "vaf")

a=merge(a1, a2, all=TRUE, sort=TRUE)
a1=a
a2=a


names = c(sample_1, sample_2)

#------------------------------------
#3d clustering using three samples:

sc = sciClone(vafs=list(v1,v2),
              copyNumberCalls=list(c1,c2),
              annotation=list(a1,a2),
              minimumDepth = 20,
              sampleNames = names,
              useSexChrs=TRUE, 
              maximumClusters=2,
              doClusteringAlongMargins=FALSE)

file <- paste( './output/', sample_1, '_vs_', sample_2, '_', 'clusters.xls', sep = '' )
writeClusterTable(sc, file)
file <- paste( './output/', sample_1, '_vs_', sample_2, '_', '1d.pdf', sep = '' )
sc.plot1d(sc, file)
file <- paste( './output/', sample_1, '_vs_', sample_2, '_', '2d.pdf', sep = '' )
sc.plot2d(sc, file)

#sc.plot3d(sc, sc@sampleNames, size=700, outputFile="glmBM.3d.gif")


## prepare clonevol input
vafs= data.frame(cluster=sc@vafs.merged$cluster,
                 vaf_1=sc@vafs.merged[,5],
                 vaf_2=sc@vafs.merged[,11],
                 stringsAsFactors=F)
#write.table(vafs, "temp.xls")

#vafs=read.table("temp.xls", header=T)

vafs = vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]
#vafs[2:3] = vafs[2:3] * 100
names(vafs)[2:length(names(vafs))] = names


## run clonevol
res = infer.clonal.models(variants=vafs, cluster.col.name="cluster", vaf.col.names=names,
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          cluster.center="mean", num.boots=1000, founding.cluster=1,
                          min.cluster.vaf=0.01, sum.p=0.001, alpha=0.1, random.seed=63108)
# new clonevol
res = convert.consensus.tree.clone.to.branch(res, branch.scale='sqrt')

plot.clonal.models(res, box.plot=TRUE, fancy.boxplot=TRUE, cell.plot=TRUE,
                   out.format="pdf", overwrite.output=TRUE, scale.monoclonal.cell.frac=TRUE,
                   cell.frac.ci=TRUE, tree.node.shape="circle", tree.node.size=40,
                   tree.node.text.size=0.65, width=8, height=5, out.dir=".")
dev.off()

f = generateFishplotInputs(results=res)
fishes = createFishPlotObjects(f)

fishes[[1]]@frac.table

file <- paste( './output/', sample_1, sample_2, 'fish.pdf', sep = '_' )
pdf(file, width=8, height=5)
par(mar = par()$mar + c(0,0,3,3))
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  print(fish)
  fishPlot(fish,shape="spline", title.btm=sample, cex.title=1.2,
           vlines = seq(1, length(names)), col.vline = "white",
           #bg.col = c("#F1F2F2","#F1F2F2","#F1F2F2"),
           border = 0.1,
           vlab=names, 
           pad.left=0.25)
  drawLegend(fish)
}
par(xpd = T)
legend("bottomright", 
       inset=c(.65,-.2), 
       pch=16, bty="n", 
       col=f$clonevol.clone.colors[2], 
       text.col = f$clonevol.clone.colors[1],
       legend = paste0(round(f$cell.fractions[[1]][2,1], 2),"%"))
legend("bottomright", 
       inset=c(-.1,-.2), 
       pch=16, bty="n", 
       col=f$clonevol.clone.colors[2], 
       text.col = f$clonevol.clone.colors,
       legend = paste0(round(f$cell.fractions[[1]][2,2], 2),"%"))

dev <- dev.off()


}



