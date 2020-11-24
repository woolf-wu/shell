# Step0 Before starting your project --------------------------------------

rm(list = objects( all = TRUE ))

if (!is.null( dev.list() )) dev.off()

sysPackages <- (.packages())

options()$download.file.method

.libPaths()

Sys.setlocale("LC_ALL","English")
.libPaths( c( "E:/R-packages", "C:/Program Files/R/R-3.6.3/library") )


# Step1 Starting your project -------------------------------------------

geneData <- read.table(file = "input/cytoband.xls",
                  header = TRUE, sep = '\t',
                  row.names = 1,
                  stringsAsFactors = FALSE,
                  comment.char = "#")
phenoData <- read.table(file = "input/sample_107_type.list",
                        header = TRUE, sep = '\t',
                        stringsAsFactors = FALSE,
                        comment.char = "#")

dim(phenoData)
dim(geneData)

a <- geneData[, phenoData[phenoData$tumor_stage == "NI",1]]
b <- as.data.frame(apply(a, 1, sum))
c <- geneData[, phenoData[phenoData$tumor_stage == "PI",1]]
d <- as.data.frame(apply(c, 1, sum))
e <- cbind(b,d)
e$gene <- rownames(e)
colnames(e) = c("NI_count", "PI_count", "gene")

# step5 COX analysis ------------------------------------------------------
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model

pheno <- phenoData
t_exp <- geneData[,pheno$ID]

dim(t_exp)
dim(pheno)

## single
library(survival)
group_data <- apply(t_exp , 1 , function(gene){
  name <- rownames(gene)
  gene <- unlist(gene)
  group <- ifelse(gene == "0", '0-none', '1-hit')
  names(group) <- name
  return(group)
})
group_data <- as.data.frame(group_data, stringsAsFactors = F)
survival_dat <- data.frame(status = pheno$status,
                           time = pheno$time,
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
options(scipen = 100)
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
res_single$gene <- rownames(res_single)
res <- merge(e, res_single)
write.table(res, file = "CNV_DEL_cox.xls", quote = F, sep = "\t", row.names = F)

table(res_single$p.value <= 0.05)
res_single <- res_single[res_single$p.value <= 0.05, ]
res_single <- res_single[order(res_single$p.value), ]
single_pick <- rownames(res_single)
res_single[single_pick,]

library(formattable)
sign_formatter <- formatter("span", 
                            style = x ~ style(color = ifelse(x > 0, "lightpink", 
                                                             ifelse(x < 0, "lightblue", "black"))))
sign_formatter(c(-1, 0, 1))
formattable( res_single[single_pick,],
             list(coef = sign_formatter,
                  p.value = color_bar("lightblue")))

pick_gene <- single_pick
pick_gene

## multi
multi_results <- apply(t_exp[pick_gene, ] , 1 , function(gene){
  #gene <- t_exp[1, ]
  gene <- unlist(gene)
  group <- ifelse(gene == "0", '0-none', '1-hit')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  res.cox <- coxph(Surv(time, status) ~ group, 
                   data =  survival_dat)
  ## summary(res.cox)
  x <- summary(res.cox)
  p.value <- signif(x$coef['group1-hit',5], digits = 2)
  beta <- signif(x$coef['group1-hit',1], digits = 2)
  HR <- signif(x$coef['group1-hit',2], digits = 2)
  HR.confint.lower <- signif(x$conf.int['group1-hit',"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int['group1-hit',"upper .95"], 2)
  HR <- paste0(HR, " (", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, p.value)
  names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
  return(res)
})
multi_results <- as.data.frame(t(as.data.frame(multi_results)))
table(multi_results$p.value <= 0.05)
res_multi <- multi_results[multi_results$p.value <= 0.05, ]
res_multi <- res_multi[order(res_multi$p.value), ]
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
  group <- ifelse(gene == "0", 'none', 'hit')
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
pick_gene <- rownames(p.v)[p.v$p.value < 0.05]
pick_gene

model_exp <- t_exp[pick_gene, ]
os.paint <- apply(model_exp, 1, function(x){
  gene <- unlist(x)
  group <- ifelse(gene == "0", '0-none', '1-hit')
  survival_dat <- data.frame(group = group,
                             status = pheno$status,
                             time = pheno$time,
                             stringsAsFactors = F)
  fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
  ggsurvplot(fit, data = survival_dat,
             surv.median.line = "hv",
             legend.title = "Group",
             legend.labs = c("hit", "none"),
             pval = TRUE,
             conf.int = TRUE,
             palette = "jco",
             ggtheme = theme_bw(),
             ylab = "progression-free probability"
  )
})

arrange_ggsurvplots(os.paint[1:6], print = TRUE, ncol = 3, nrow = 2)
arrange_ggsurvplots(os.paint[7:12], print = TRUE, ncol = 3, nrow = 2)
arrange_ggsurvplots(os.paint[13:14], print = TRUE, ncol = 3, nrow = 2)

coef_dat <- as.data.frame(multi_results[pick_gene, ])

exp <- as.data.frame(cbind("(", coef_dat$coef, "* coef_exp$", rownames(coef_dat), ")+"))
exp <- apply(exp, 1, function(x){
  Reduce('paste', x)
  })
expression <- sub("\\+$", "", Reduce('paste0', exp, 'Risk score='))
expression <- gsub("\\$ ", "$", expression)
expression <- gsub("\\(", " (", expression)
expression <- gsub("\\)", " ) ", expression)
expression


# step ROC ----------------------------------------------------------------

exp <- t_exp

model_exp <- t(exp[pick_gene,])
colnames(model_exp) <- pick_gene
dat <- cbind(pheno, model_exp)

dat$gender <- factor(dat$gender)
dat$tumor_stage <- factor(dat$tumor_stage)


## ROC
library("survminer")
colnames(dat)
s <- Surv(time, status) ~ FAM75A6 + ANKRD20A2 + ANKRD20A3 + AQP7P3 + FAM74A2 + FAM74A4 + FAM75A5 + FAM75A7 + FAM95B1 + FOXD4L2 + FOXD4L4 + LOC286297 + LOC642929 + LOC643648

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
                     times = c(36, 60),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC, time = 60, col = "red", add = F)
plot(ROC, time = 60, col = "blue", add = T)
legend("top",
       legend = c("Testing Datasets","TCGA Datasets"),
       ncol = 4,
       cex = 0.8,
       bty = "n",
       col = c("red","blue"),
       lty = 1,lwd = 2) 


ROC$AUC
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
exp_dat <- dat[names(sort(fp)), 8:12]

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





# step USER GROUP ------------------------------------------------------------

exp <- t_exp
exp <- v_exp
exp <- AssayData

exp <- exp[apply(exp, 1, function(x) sum(x == 0)) < (ncol(exp)*0.49), ]
pheno <- phenoData[colnames(exp), ]
model_exp <- exp[pick_gene, ]
coef_exp <- as.data.frame(t(model_exp))

coef_exp$Risk_score=( 1 * coef_exp$LINC01554 )+( -0.77 * coef_exp$MIR4435_2HG )+( 0.74 * coef_exp$LINC01537 )+( -0.75 * coef_exp$LINC01139 )+( -0.78 * coef_exp$MAPKAPK5_AS1 )

median(coef_exp$Risk_score) 

coef_exp$group <- ifelse(coef_exp$Risk_score <= -499.85, 'high', 'low')
table(coef_exp$group)

coef_exp <- coef_exp[rownames(pheno),]
survival_dat <- data.frame(group = coef_exp$group,
                           status = pheno$status,
                           time = pheno$time,
                           age = pheno$age,
                           gender = pheno$gender,
                           stage = pheno$tumor_stage,
                           stringsAsFactors = F)
median(survival_dat[survival_dat$group == 'high', 3])
median(survival_dat[survival_dat$group == 'low', 3])
fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
a <- ggsurvplot(fit, data = survival_dat,
           surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("High", "Low"),
           pval = TRUE,
           conf.int = TRUE,
           palette = "jco",
           ggtheme = theme_bw()
)

list_su <- list()
list_su[[3]] <- a
arrange_ggsurvplots(list_su, print = TRUE, ncol = 3, nrow = 1)

## STAGE

exp <- t_exp
exp <- exp[apply(exp, 1, function(x) sum(x == 0)) < (ncol(exp)*0.49), ]
pheno <- phenoData[colnames(exp), ]

table(pheno$tumor_stage)

e_phe <- pheno[pheno$tumor_stage == "early_stage", ]
l_phe <- pheno[pheno$tumor_stage == "later_stage", ]

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
dat$tumor_stage <- factor(dat$tumor_stage)

library("survminer")
colnames(dat)
s <- Surv(time, status) ~ LINC01554 + MAPKAPK5_AS1 + MIR4435_2HG + LINC01537 + LINC01139

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
                     times = c(36, 60),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC, time = 60, col = "blue", add = F)
plot(ROC, time = 60, col = "red", add = F)
legend("top",
       legend = c("Early subgroup","Late subgroup"),
       ncol = 4,
       cex = 0.8,
       bty = "n",
       col = c("blue","red"),
       lty = 1,lwd = 2) 

ROC$AUC
confint(ROC)


model_exp <- exp[pick_gene,] + 1
coef_exp <- as.data.frame(t(model_exp))

coef_exp$Risk_score=( 1 * coef_exp$LINC01554 )+( -0.77 * coef_exp$MIR4435_2HG )+( 0.74 * coef_exp$LINC01537 )+( -0.75 * coef_exp$LINC01139 )+( -0.78 * coef_exp$MAPKAPK5_AS1 )

median(coef_exp$Risk_score) 

coef_exp$group <- ifelse(coef_exp$Risk_score <= -453.43, 'high', 'low')
table(coef_exp$group)

coef_exp <- coef_exp[rownames(pheno),]
survival_dat <- data.frame(group = coef_exp$group,
                           status = pheno$status,
                           time = pheno$time,
                           age = pheno$age,
                           gender = pheno$gender,
                           stage = pheno$tumor_stage,
                           stringsAsFactors = F)
median(survival_dat[survival_dat$group == 'high', 3])
median(survival_dat[survival_dat$group == 'low', 3])
fit <- survfit(Surv(time, status) ~ group, data = survival_dat)
ggsurvplot(fit, data = survival_dat,
           surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("High", "Low"),
           pval = TRUE,
           conf.int = TRUE,
           palette = "jco",
           ggtheme = theme_bw()
)


## single
univ_formulas <- sapply(colnames(survival_dat),
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
res_single


## multi

res.cox <- coxph(Surv(time, status) ~ group + age + gender + stage, 
                 data =  survival_dat)
## summary(res.cox)
x <- summary(res.cox)
p.value <- signif(x$coef[,5], digits = 2)
beta <- signif(x$coef[,1], digits = 2)
HR <- signif(x$coef[,2], digits = 2)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
HR <- paste0(HR, " (", 
             HR.confint.lower, "-", HR.confint.upper, ")")
res_multi <- as.data.frame(cbind(beta, HR, p.value))
colnames(res_multi) <- c("coef", "HR (95% CI for HR)", "p.value")
res_multi

t_cox <- cbind(res_single, res_multi)
v_cox <- cbind(res_single, res_multi)
a_cox <- cbind(res_single, res_multi)

t_cox
v_cox
a_cox

all_cox <- rbind(t_cox, v_cox, a_cox)
all_cox1 <- rbind(t_cox, v_cox, a_cox)
all_cox


