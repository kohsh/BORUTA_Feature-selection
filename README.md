 # ✒️ BORUTA_Feature-selection

Our Guthub repo include pipeline for running random forest-based algorithm on [SRA](https://www.ncbi.nlm.nih.gov/sra) fastq files using [Boruta](https://www.jstatsoft.org/article/view/v036i11) tool. The pipeline is written in an [R](https://github.com/rstudio/rstudio) programming language.

## ⚙️ Technologies & Tools

![](https://img.shields.io/badge/Code-RScript-informational?style=flat&logo=<#FF6000>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-Rstudio-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-GitHub-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)
![](https://img.shields.io/badge/Tools-SRAtoolkit-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)

```R

#!/usr/bin/Rscript

#Load packages
library(tidyverse)
library(Boruta)


set.seed(123)
setwd("/home/data")

# read in reads file
data <- read.delim("input.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)

# Filtering lowely express genes
group <- paste(c(rep("AD",200), rep("Ctrl",400)))
y <- DGEList(counts= data, group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] 
y <- calcNormFactors(y,method = "TMM") # TMM normalization method
reads_norm <- cpm(y, normalized.lib.sizes = T)

# log transform the reads to make the counts closer to a normal distribution
reads_norm <- log2(reads_norm+1)
datExpr=t(reads_norm)
dim(datExpr)


# Now we read in the physiological trait data
datTraits <- read.delim("meta_df.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)
datTraits <- dplyr::select(datTraits, c('Sample', 'Project_Numeric', 'Strategy_Numeric', 'Condition', 'Age', 'Sex'))
head(datTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$Sample
# show that row names agree
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly.

# select covariate columns for machine learning
var <- c('Sample', 'Condition', 'Age', 'Sex')
ml_covar <- datTraits[,names(datTraits) %in% var]

# combine genes and covariates into one object 
covar_genes <- cbind(datExpr, ml_covar)

# save full cohort with covariate columns for machine learning regression
saveRDS(covar_genes, "ml_regression_input_full_cohort.rds")
```

### Input data format

|   | AD01 | AD02 | AD03 | Ctrl01 | Ctrl02 |Ctrl03 |
| ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |
| ENST00000456328.2  | 7.85718  | 3.52895 | 2.15625 |0 | 0| 0.466484 |
| ENST00000450305.2  | 1.2187  | 1.39157 |  3.47069  | 0.0013632| 0.0355961| 0.0544532|
 
```R

# balanced sub cohort for machine learning classification
set.seed(123)
mldata <- covar_genes

mldata <- mldata[complete.cases(mldata),] #Return a logical vector indicating which cases are complete, i.e., have no missing values.
mldata <- droplevels(mldata) #The function droplevels is used to drop unused levels from a factor or, more commonly, from factors in a data frame.
contrast_labels <- mldata$Condition

rm_cols <- var <- c('Sample', 'Condition', 'Age', 'Sex')
mldata <- mldata[,!(names(mldata) %in% rm_cols)]

saveRDS(mldata,"data_case_vs_ctrl_ml_class.rds")
saveRDS(contrast_labels,"labels_case_vs_ctrl_ml_class.rds")

# input data set up - the labels and data files are produced by the "prepare_for_m_learn" script
labels <- readRDS("labels_case_vs_ctrl_ml_class.rds")
data <- readRDS("data_case_vs_ctrl_ml_class.rds")

# Case VS Control
data_labels <- cbind(data, labels)
boruta_model=Boruta(labels ~ ., data=data_labels, doTrace = 2 , maxRuns = 20000)
print(boruta_model)
saveRDS(boruta_model,file="boruta_model.rds")
final_boruta_model=TentativeRoughFix(boruta_model)
Case_VS_Ctrl_model_sel_feat=getSelectedAttributes(final_boruta_model)
saveRDS(Case_VS_Ctrl_model_sel_feat,file="Case_VS_Ctrl_model_sel_feat.rds")
write.csv(Case_VS_Ctrl_model_sel_feat,file="Case_VS_Ctrl_model_selelected_feat.csv")

k <-lapply(1:ncol(boruta_model$ImpHistory),function(i)
  boruta_model$ImpHistory[is.finite(boruta_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_model$ImpHistory)
boruta_model_feature_rank <- sort(sapply(k,median))
write.csv(boruta_model_feature_rank,file="boruta_model_feature_rank.csv")
```

### Onput data format

| ensembl_transcript_id | Z.score | 
| ------------- | ------------- |
| ENST00000624169.1  | 1.032175543  |
| ENST00000534336.1  | 1.042689192  |



```R

# Stack-plot for max-Zscores
setwd(path.expand("/home/STACK-plot") )
ENR <- read.delim("boruta_model_feature_rank.txt", header=TRUE, fill = TRUE)
ggplot(ENR, aes(x=hgnc_symbol, y=Z_Score, fill=Contrast)) + geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(ENR))) +
theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
ggsave("AD-Ctrl.eps")
```

### Z.score stack-plot

![AD-vs-Ctrl](https://user-images.githubusercontent.com/46553150/218000536-adf416bf-cc60-4f4e-9ffa-3abb9ebfacca.png)


