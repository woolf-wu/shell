l1 = read.table("./lancet.xls", header = T)
s1 = read.table("./Strelka2_manta.xls", header = T)
m1 = read.table("./muTect_Strelka.xls", header = T)

colnames(l1) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")
colnames(s1) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")
colnames(m1) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")

data=rbind(l1,s1,m1)
data=data[data$Tumor_Sample_Barcode!="Tumor_Sample_Barcode",]
data$t_AF = as.numeric(data$t_AF)
data=data[,c(3,2,4)]

type_x = "Tumor_Sample_Barcode"
type_y = "t_AF"
type_z = "Type"
color = c("#8ebc94", "#6aafd8", "#ec6a5c")
ggboxplot(data,
            x = type_z, y = type_y,
            fill = type_z,
            color = type_z,
            palette = color,
            facet.by = type_x,
            short.panel.labs = TRUE,
            bxp.errorbar = TRUE,
            bxp.errorbar.width = 0.6,
            xlab = FALSE,
            add = "boxplot",
            add.params = list(color = "black",
                              linetype = "solid")) + 
  theme(axis.title=element_text(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=14,vjust= 0.8),
        legend.position = "top")  + 
  stat_compare_means(label = "p.format", paired = TRUE)












l1 = read.table("./lancet.xls", header = T)
s1 = read.table("./Strelka2_manta.xls", header = T)
m1 = read.table("./muTect_Strelka.xls", header = T)

data=merge(l1,s1,by=c("Chromosome_Start_position", "Tumor_Sample_Barcode"), all=F, sort=TRUE)
data=merge(data,m1,by=c("Chromosome_Start_position", "Tumor_Sample_Barcode"), all=F, sort=TRUE)

v1=data[,c(1,2,3,4)]
v2=data[,c(1,2,5,6)]
v3=data[,c(1,2,7,8)]


colnames(v1) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")
colnames(v2) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")
colnames(v3) = c("Chromosome_Start_position","Tumor_Sample_Barcode","t_AF","Type")


data=rbind(v1,v2,v3)

data=data[data$Tumor_Sample_Barcode!="Tumor_Sample_Barcode",]
#data=data[data$Tumor_Sample_Barcode=="B1",]
data$t_AF = as.numeric(data$t_AF)

type_x = "Tumor_Sample_Barcode"
type_y = "t_AF"
type_z = "Type"
color = c("#8ebc94", "#6aafd8", "#ec6a5c")
ggboxplot(data,
          x = type_z, y = type_y,
          fill = type_z,
          color = type_z,
          palette = color,
          facet.by = type_x,
          short.panel.labs = TRUE,
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.6,
          xlab = FALSE,
          add = "boxplot",
          add.params = list(color = "black",
                            linetype = "solid")) + 
  theme(axis.title=element_text(),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=14,vjust= 0.8),
        legend.position = "top")  + 
  stat_compare_means(label = "p.format", paired = TRUE)  + 
  scale_y_continuous(limits=c(0, 1.5))

library(gridExtra)

for (sample in c("S1", "A1", "L1", "L2", "B1", "L3", "H3", "H1", "H2")) {
  print(sample)
  a=data[data$Tumor_Sample_Barcode==sample,]
  a = a[,c(4,3)]
  a$Type = as.factor(a$Type)
  pic=paste(sample)
  pic=ggplot(a, aes(x=t_AF, fill=Type)) + 
    geom_histogram(aes(y=..density..),
                   alpha=0.5, 
                   position="identity")+
    geom_density(alpha=.2)  + 
    facet_grid(Type ~ .) +
    theme(legend.position='none') +
    scale_fill_manual(values=color)
  print(pic)
}








sample_list=c("A1", "B1", "H1", "H2", "H3", "L1", "L2", "L3", "S1")

sample= sample_list[1]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_1 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 

sample= sample_list[2]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_2 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 

sample= sample_list[3]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_3 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 



sample= sample_list[4]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_4 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 


sample= sample_list[5]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_5 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 


sample= sample_list[6]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_6 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 


sample= sample_list[7]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_7 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 


sample= sample_list[8]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_8 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 


sample= sample_list[9]
print(sample)
a=data[data$Tumor_Sample_Barcode==sample,]
a = a[,c(4,3)]
a$Type = as.factor(a$Type)
file=paste("pic", i, sep = "_")
pic_9 <- ggplot(a, aes(x=t_AF, fill=Type)) + 
  geom_histogram(aes(y=..density..),
                 alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2)  + 
  facet_grid(Type ~ .) +
  ggtitle(sample) +
  theme(legend.position='none',
        plot.title=element_text(hjust = 0.5)) +
  scale_fill_manual(values=color) 

grid.arrange(pic_1,pic_2,pic_3,pic_4,pic_5,pic_6,pic_7,pic_8,pic_9,nrow=3,ncol=3)
