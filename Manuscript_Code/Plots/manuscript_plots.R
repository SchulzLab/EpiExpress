library(ggplot2)
library("rjson")
library(ppcor)
library(pheatmap)
library(ggforce)
library(ggdist)
library(RColorBrewer)
library(gridExtra)
library(dplyr)

jsonCol <- fromJSON(file= "~/ColourDict.json")
font_config <- list(family= "Helvetica", size= 7)
colors <- unlist(jsonCol$models)

kept_genes <- read.table("~/kept_genes.txt", sep= "\t", header= T)[, 2]
gene_feat <- read.table("~/GeneFeatureMat.txt.gz", sep= "\t", comment.char= "!", header= T)
hits <- which(gene_feat$Ensembl.ID %in% kept_genes)
gene_feat <- gene_feat[hits, ]

biotype_cnt <- table(gene_feat$Biotype_general)
biotype_df <- data.frame(biotype= names(biotype_cnt), frequency= round(as.numeric(biotype_cnt) / nrow(gene_feat) * 100, 2))
myPalette <- brewer.pal(5, "Set2")

na_hits <- which(is.na(gene_feat$X.ENCODE_peCREs) == T)
if(length(na_hits)){
        gene_feat <- gene_feat[-na_hits, ]
}
gene_feat_trimmed <- gene_feat[hits, c("X.TSS", "X.Transcripts", "Gene_length", "Exons_length", "Gene_density", "Expression_ubiquitousness", "Fraction_nonzero_in_train", "Fraction_nonzero_in_test", "Distance_TAD_border", "X.ENCODE_peCREs")]

gene_feat_trimmed <- na.omit(gene_feat_trimmed)
mat <- matrix(NA, ncol= ncol(gene_feat_trimmed), nrow= ncol(gene_feat_trimmed)); for(i in seq(ncol(gene_feat_trimmed))){mat[i, ] <- sapply(seq(ncol(gene_feat_trimmed)), function(j)cor(gene_feat_trimmed[, i], gene_feat_trimmed[, j], method= "s"))}
colnames(mat) <- rownames(mat) <- colnames(gene_feat_trimmed)

### Build the dataframe for storing the model performance
kept_genes_df <- data.frame(gene= kept_genes)

## NN ##
nn <- read.csv("~/MLP_Best_seeds_manuscript.csv") #updated MLP performance
nn_original <- nn
colnames(nn)[which(colnames(nn) == "MSE_S_test")] <- "NN_MSE_S"
colnames(nn)[which(colnames(nn) == "Cor_test")] <- "NN_Cor_test"
##############

## guidedRF ##
all_cor_res <- merge(kept_genes_df, nn, by.x= "gene", by.y= "Name", all.x= T)
peaks_RF <- read.csv("~/RF_1MB_all_H38_Final.csv")
colnames(peaks_RF) <- c("gRF_MSE_test", "gRF_Cor_test", "gRF_MSEs_back_test", "Name")
all_cor_res <- merge(all_cor_res, peaks_RF, by.x= "gene", by.y= "Name", all.x= T)
##############

## unguidedRF ##
binned_RF_cor <- read.csv("~/binned_activity_w1MB_100bs_RF_models_cor_summary.csv", sep= "\t", header= T)[, -1]
binned_RF_err <- read.csv("~/binned_activity_w1MB_100bs_RF_models_err_summary.csv", sep= "\t", header= T)[, -1]
binned_RF <- merge(binned_RF_cor, binned_RF_err, by.x= "Ensembl_ID", by.y= "Ensembl_ID")
colnames(binned_RF)[2:5] <- sapply(colnames(binned_RF)[2:5], function(str)paste0("ugRF_", str))
all_cor_res <- merge(all_cor_res, binned_RF, by.x= "gene", by.y= "Ensembl_ID", all.x= T)
##############

## STITCHIT ##
stitchit <- read.table("~/performance_test_partition0_all_genes.tsv",  header= T)
colnames(stitchit)[2:3] <- c("STITCHIT_cor", "STITCHIT_MSE")
all_cor_res <- merge(all_cor_res, stitchit, by.x= "gene", by.y= "GeneID", all.x= T)
##############

cnn_old <- read.csv("~/CNN_validation_df_warmseeds_90_percent_training.txt", sep= "\t", header= T)[, -1]
cnn <- read.csv("~/CNN_validation_df_warmseeds_90_percent_training_sigmoid.txt", sep= "\t", header= T)[, -1]

#MSE
CREMLP_kept <- nn[which(nn$Name  %in% kept_genes), ]$NN_MSE_S
CRERF_kept <- peaks_RF[which(peaks_RF$Name %in% kept_genes), ]$gRF_MSEs_back_test
stitchit_kept <- stitchit[which(stitchit$GeneID %in% kept_genes), ]$STITCHIT_MSE
binnedRF_kept <- binned_RF_err[which(binned_RF_err$Ensembl_ID %in% kept_genes), ]$test_err
cnn_kept <- cnn[which(cnn$gene %in% kept_genes), ]$test_err

#cor
CREMLP_kept_cor <- nn[which(nn$Name %in% kept_genes), ]$NN_Cor_test
CRERF_kept_cor <- peaks_RF[which(peaks_RF$Name %in% kept_genes), ]$gRF_Cor_test
stitchit_kept_cor <- stitchit[which(stitchit$GeneID %in% kept_genes), ]$STITCHIT_cor
binnedRF_kept_cor <- binned_RF_cor[which(binned_RF_cor$Ensembl_ID %in% kept_genes), ]$test_cor
cnn_kept_cor <- cnn[which(cnn$gene %in% kept_genes), ]$test_cor

#MSE
model_names <- c("CRE-MLP", "CRE-RF", "Binned-RF", "STITCHIT", "Binned-CNN")
mse_list <- list(CREMLP_kept, CRERF_kept, binnedRF_kept, stitchit_kept, cnn_kept)
mse_df <- NULL;
for(i in seq(length(mse_list))){
        mse_df <- rbind(mse_df, data.frame(model= model_names[i], MSE= mse_list[[i]]))
}

#cor
model_names <- c("CRE-MLP", "CRE-RF", "Binned-RF", "STITCHIT", "Binned-CNN")
cor_list <- list(CREMLP_kept_cor, CRERF_kept_cor, binnedRF_kept_cor, stitchit_kept_cor, cnn_kept_cor)
cor_df_all <- NULL;
for(i in seq(length(cor_list))){
        cor_df_all <- rbind(cor_df_all, data.frame(model= model_names[i], cor= cor_list[[i]]))
}

cor_df_all$cor[which(is.na(cor_df_all$cor)) ] <- 0

mse_df$model <- factor(mse_df$model, levels= names(colors))
cor_df_all$model <- factor(cor_df_all$model, levels= names(colors))

mse_df$MSE <- as.numeric(format(round(mse_df$MSE, 3)))
cor_df_all$cor <- as.numeric(format(round(cor_df_all$cor, 2)))

########################## Suppl. Fig. 3a ################################
pdf("~/model_specific_gene_set_pearson.pdf", width=3.5, height=1.6)
print(ggplot(cor_df_all, aes(x= model, y= cor, fill= model)) + geom_violin() + scale_fill_manual(values= colors) + theme_bw() +  ylim(-0.3,1) + 
        ylab("Pearson correlation") + theme(legend.position="none", axis.text=element_text(family= font_config$family, size= font_config$size),axis.title=element_text(family= font_config$family, size= font_config$size), text=element_text(family= font_config$family, size= font_config$size), panel.grid.major= element_line(color = "white"), panel.grid.minor= element_line(color = "white")) + 
        geom_boxplot(width=.1, outlier.shape = NA))
dev.off()

model_names <- c("CRE-MLP", "CRE-RF", "Binned-RF", "STITCHIT", "Binned-CNN")
number_list <- c(length(CREMLP_kept), length(CRERF_kept), length(binnedRF_kept), length(stitchit_kept), length(cnn_kept))
number_df <- data.frame(model_names, number_list)
colnames(number_df) <- c("model", "numberGenes")
number_df$model <- factor(number_df$model, levels= names(colors))

########################## Suppl. Fig. 3b ################################
pdf("~/1MB_total_number_models_barplots.pdf", width=1.75, height=1.8)
print(ggplot(number_df, aes(x= model, y= numberGenes, fill= model)) + geom_bar(stat= "identity") + scale_y_continuous(breaks = c(5000, 10000, 15000, 20000, 25000)) + scale_fill_manual(values= colors)+ theme_bw() + xlab(NULL) + ylab("Number of genes") + theme(text=element_text(family= font_config$family, size= font_config$size),  axis.text.x = element_text(angle = 90), legend.position = "none", panel.grid.major = element_line(color = "white"), panel.grid.minor = element_line(color = "white")))
dev.off()

colnames(cnn_old)[which(colnames(cnn_old) == "test_cor")] <- "CNNrelu_test_cor"
colnames(cnn_old)[which(colnames(cnn_old) == "test_err")] <- "CNNrelu_test_err"
colnames(cnn)[which(colnames(cnn) == "test_cor")] <- "CNN_test_cor"
colnames(cnn)[which(colnames(cnn) == "test_err")] <- "CNN_test_err"
all_cor_res_CNNold <- merge(all_cor_res, cnn_old, by.x= "gene", by.y= "gene", all.x= T)
all_cor_res <- merge(all_cor_res, cnn, by.x= "gene", by.y= "gene", all.x= T)

best2_df <- merge(all_cor_res[, c("gene", "gRF_Cor_test", "CNN_test_cor", "gRF_MSEs_back_test", "CNN_test_err")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)

## replace all NA's with zero
for(j in seq(ncol(all_cor_res))){
        na_hits <- which(is.na(all_cor_res[, j]));
        if(length(na_hits)){
                if(grepl("MSE", colnames(all_cor_res)[j], fixed= T) || grepl("err", colnames(all_cor_res)[j], fixed= T)){
                        print("Skipping the MSE columns")
                }
                else{
                        all_cor_res[na_hits, j] <- 0
                }
        }
}
##############
for(j in seq(ncol(all_cor_res))){
        if(grepl("MSE", colnames(all_cor_res)[j], fixed= T) || grepl("err", colnames(all_cor_res)[j], fixed= T)){
                err_na_hits <- which(is.na(all_cor_res[, j]))
                maxs <- vector("numeric", length(err_na_hits))
                if(length(err_na_hits)){
                        for(k in seq(length(err_na_hits))){
                                maxs[k] <- max(all_cor_res[err_na_hits[k], c("NN_MSE_S", "gRF_MSEs_back_test", "ugRF_test_err", "STITCHIT_MSE", "CNN_test_err")], na.rm= T)
                                if(is.na(maxs[k])){
                                        print(all_cor_res[, "gene"][err_na_hits[k]])
                                        print(err_na_hits[k])
                                        maxs[k] <- 50
                                }else{
                                        maxs[k] <- maxs[k] + 1
                                }
                        }
                        all_cor_res[err_na_hits, j] <- maxs
                }
        }
}
##############

cor_mat <- all_cor_res[, c("gene", "NN_Cor_test", "gRF_Cor_test", "ugRF_test_cor", "STITCHIT_cor", "CNN_test_cor")]
err_mat <- all_cor_res[, c("gene", "NN_MSE_S", "gRF_MSEs_back_test", "ugRF_test_err", "STITCHIT_MSE", "CNN_test_err")]

colnames(cor_mat)[colnames(cor_mat) == "NN_Cor_test"] <- "CRE-MLP"
colnames(cor_mat)[colnames(cor_mat) == "gRF_Cor_test"] <- "CRE-RF"
colnames(cor_mat)[colnames(cor_mat) == "ugRF_test_cor"] <- "Binned-RF"
colnames(cor_mat)[colnames(cor_mat) == "STITCHIT_cor"] <- "STITCHIT"
colnames(cor_mat)[colnames(cor_mat) == "CNN_test_cor"] <- "Binned-CNN"

colnames(err_mat)[colnames(err_mat) == "NN_MSE_S"] <- "CRE-MLP"
colnames(err_mat)[colnames(err_mat) == "gRF_MSEs_back_test"] <- "CRE-RF"
colnames(err_mat)[colnames(err_mat) == "ugRF_test_err"] <- "Binned-RF"
colnames(err_mat)[colnames(err_mat) == "STITCHIT_MSE"] <- "STITCHIT"
colnames(err_mat)[colnames(err_mat) == "CNN_test_err"] <- "Binned-CNN"
##############

#order methods
cor_mat<- cor_mat[ , c("gene", "CRE-RF","CRE-MLP","STITCHIT", "Binned-RF", "Binned-CNN")]
err_mat<- err_mat[ , c("gene", "CRE-RF","CRE-MLP","STITCHIT", "Binned-RF", "Binned-CNN")]

########################## Suppl. Fig. 3d ################################
my_line <- function(x,y,...){
        my_col= rep("blue", times= length(x))
        neg_idx <- which(x < 0 | y < 0)
        my_col[neg_idx] <- "grey"
        points(x,y, col= my_col, ...)
        abline(a = 0, b = 1, col= "red", ...)
}

pdf("~/1MB_model_scatter_pairs_MSE.pdf")
par(family= font_config$family)
pairs(err_mat[, -1], lower.panel= my_line, upper.panel= my_line, pch= 20)
dev.off()

###############

### cor
all_models_list_cor <- list()
for(i in seq(2, ncol(cor_mat))){ ## Starts from 2 to skip the gene ID in the dataframe
        name <- colnames(cor_mat)[i]
        all_models_list_cor[[name]] <- cor_mat[, i]
}
### MSE
all_models_list_err <- list()
for(i in seq(2, ncol(err_mat))){ ## Starts from 2 to skip the gene ID in the dataframe
        name <- colnames(err_mat)[i]
        all_models_list_err[[name]] <- err_mat[, i]
}

get_com_mat <- function(all_models_list, metric){
        com_mat <- matrix(NA, nrow= length(all_models_list), ncol= length(all_models_list))
        rownames(com_mat) <- colnames(com_mat) <- names(all_models_list)
        onevRest_res <- vector("numeric", length(all_models_list))
        onevRest_res_genes <- list()
        for(i in seq(length(all_models_list))){
                x <- all_models_list[[i]]
                diffs <- lapply(all_models_list[seq(length(all_models_list))[-i]], function(l) sign(x - l))
                diffs_mat <- do.call(cbind, diffs)
                if(metric == "cor"){
                        strictly_sup_genes <- which(rowSums(diffs_mat) == (length(all_models_list) - 1))
                        strictly_sup <- length(which(rowSums(diffs_mat) == (length(all_models_list) - 1)))
                }else{
                        strictly_sup_genes <- which(rowSums(diffs_mat) == -(length(all_models_list) - 1))
                        strictly_sup <- length(which(rowSums(diffs_mat) == -(length(all_models_list) - 1)))
                }
                onevRest_res[i] <- strictly_sup
                onevRest_res_genes[[i]] <- strictly_sup_genes
                names(onevRest_res_genes)[i] <- names(all_models_list)[i]
                for(j in seq(length(all_models_list))){
                        if(i == j){
                                next;
                        }
                        y <- all_models_list[[j]]
                        if(metric == "cor"){
                                print(paste("percentage of", names(all_models_list)[i], ">", names(all_models_list)[j]))
                                print(com_mat[i, j] <- length(which(x > y & x != y)) / length(x) * 100)
                        }else{
                                print(paste("percentage of", names(all_models_list)[i], "<", names(all_models_list)[j]))
                                print(com_mat[i, j] <- length(which(x < y & x != y)) / length(x) * 100)
                        }
                }
        }
        names(onevRest_res) <- names(onevRest_res_genes) <- names(all_models_list)
        return(list(com_mat= com_mat, onevRest_res= onevRest_res, onevRest_res_genes= onevRest_res_genes))
}

com_mat_cor <- get_com_mat(all_models_list_cor, "cor")
com_mat_err <- get_com_mat(all_models_list_err, "err")

df <- data.frame(model= names(com_mat_cor$onevRest_res), NumOfGenes= com_mat_cor$onevRest_res)
df$model[df$model == "NN_Cor_test"] <- "CRE-MLP"
df$model[df$model == "gRF_Cor_test"] <- "CRE-RF"
df$model[df$model == "ugRF_test_cor"] <- "Binned-RF"
df$model[df$model == "STITCHIT_cor"] <- "STITCHIT"
df$model[df$model == "CNN_test_cor"] <- "Binned-CNN"
df$model <- factor(df$model, levels= names(colors))

########################## Main Fig. 2b ################################
pdf("~/1MB_strictly_superior_barplots.pdf", width=1.75, height=1.8)
print(ggplot(df, aes(x= model, y= NumOfGenes, fill= model)) + geom_bar(stat= "identity") +scale_fill_manual(values= colors)+ theme_bw() + xlab(NULL) + ylab("Number of superior genes") + theme(text=element_text(family= font_config$family, size= font_config$size),  axis.text.x = element_text(angle = 90), legend.position = "none", panel.grid.major = element_line(color = "white"), panel.grid.minor = element_line(color = "white")))
dev.off()

colnames(com_mat_err$com_mat) <- c("CRE-MLP", "CRE-RF", "Binned-RF", "STITCHIT", "Binned-CNN")
rownames(com_mat_err$com_mat) <- colnames(com_mat_err$com_mat)

colnames(com_mat_cor$com_mat) <- c("CRE-MLP", "CRE-RF", "Binned-RF", "STITCHIT", "Binned-CNN")
rownames(com_mat_cor$com_mat) <- colnames(com_mat_err$com_mat)

com_mat_err$com_mat <- com_mat_err$com_mat[, c("CRE-RF", "CRE-MLP", "STITCHIT", "Binned-RF", "Binned-CNN")]
com_mat_err$com_mat <- com_mat_err$com_mat[c("CRE-RF", "CRE-MLP", "STITCHIT", "Binned-RF", "Binned-CNN"), ]

com_mat_cor$com_mat <- com_mat_cor$com_mat[, c("CRE-RF", "CRE-MLP", "STITCHIT", "Binned-RF", "Binned-CNN")]
com_mat_cor$com_mat <- com_mat_cor$com_mat[c("CRE-RF", "CRE-MLP", "STITCHIT", "Binned-RF", "Binned-CNN"), ]


kakapo_palette <- colorRampPalette(c("#56B4E9", "#FEFEFE", "#D55E00"))(100)
########################## Main Fig. 2c ################################
pdf("~/1MB_Percentage_models_heatmap_MSE_kakapo.pdf", height=4, width=4.5)
pheatmap::pheatmap(com_mat_err$com_mat, display_numbers= T, color = kakapo_palette, number_color= "black", cluster_rows= F, cluster_cols= F, fontfamily= font_config$family, legend = TRUE, cellheight= 40, cellwidth= 40)
dev.off()

########################## Suppl. Fig. 3c ################################
pdf("~/1MB_Percentage_models_heatmap_cor_kakapo.pdf", height=4, width=4.5)
pheatmap::pheatmap(com_mat_cor$com_mat, display_numbers= T, color = kakapo_palette, number_color= "black", cluster_rows= F, cluster_cols= F, fontfamily= font_config$family, legend = TRUE, cellheight= 40, cellwidth= 40)
dev.off()

cnn_gf <- merge(all_cor_res[, c("gene", "validation_loss", "training_cor", "CNN_test_cor", "training_err", "CNN_test_err")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)
stitchit_gf <- merge(all_cor_res[, c("gene", "STITCHIT_cor", "STITCHIT_MSE")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)
nn_gf <- merge(all_cor_res[, c("gene", "NN_MSE_S", "NN_Cor_test")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)
pRF_gf <- merge(all_cor_res[, c("gene", "gRF_MSE_test", "gRF_MSEs_back_test", "gRF_Cor_test")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)
bRF_gf <- merge(all_cor_res[, c("gene", "ugRF_training_cor", "ugRF_test_cor", "ugRF_training_err", "ugRF_test_err")], gene_feat, by.x= "gene", by.y= "Ensembl.ID", all.x= T)

features= c("gene", "X.TSS", "Gene_length", "Exons_length", "Gene_density", "Expression_ubiquitousness", "Fraction_nonzero_in_test", "Distance_TAD_border", "longest_5prime_UTR", "longest_3prime_UTR", "X.ENCODE_peCREs") #, "X.Transcripts")
structural_features= c("X.TSS", "Gene_length", "Exons_length", "longest_5prime_UTR", "longest_3prime_UTR")
genomic_features1 <- c("Expression_ubiquitousness", "Fraction_nonzero_in_test")
genomic_features2 <- c("Gene_density", "X.ENCODE_peCREs")
biotypes <- unique(cnn_gf$Biotype_general)

### There are NA's introduced, mainly because of the peCREs and I need to remove them before computing the pcor, otherwise it crashes.
na_hits <- which(is.na(cnn_gf$X.ENCODE_peCREs) == T)
####
cnn_gf_trimmed <- cnn_gf[-na_hits, c("CNN_test_cor", features)]
cnn_gf_trimmed_err <- cnn_gf[-na_hits, c("CNN_test_err", features)]

cnn_gf_bt <- cnn_gf[, c("CNN_test_cor", "Biotype_general")]
cnn_gf_bt <- cbind(cnn_gf_bt,rep("Binned-CNN", nrow(cnn_gf_bt)))
colnames(cnn_gf_bt) <- c("cor","Biotype_general", "method") 
cnn_gf_bt_err <- cnn_gf[, c("CNN_test_err", "Biotype_general")]

####
stitchit_gf_trimmed <- stitchit_gf[-na_hits, c("STITCHIT_cor", features)]
stitchit_gf_trimmed_err <- stitchit_gf[-na_hits, c("STITCHIT_MSE", features)]

stitchit_gf_bt <- stitchit_gf[, c("STITCHIT_cor", "Biotype_general")]
stitchit_gf_bt <- cbind(stitchit_gf_bt, rep("STITCHIT", nrow(stitchit_gf_bt)))
colnames(stitchit_gf_bt) <- c("cor","Biotype_general", "method")
stitchit_gf_bt_err <- stitchit_gf[, c("STITCHIT_MSE", "Biotype_general")]

####
nn_gf_trimmed <- nn_gf[-na_hits, c("NN_Cor_test", features)]
nn_gf_trimmed_err <- nn_gf[-na_hits, c("NN_MSE_S", features)]

nn_gf_bt <- nn_gf[, c("NN_Cor_test", "Biotype_general")]
nn_gf_bt <- cbind(nn_gf_bt, rep("CRE-MLP", nrow(nn_gf_bt)))
colnames(nn_gf_bt) <- c("cor","Biotype_general", "method")
nn_gf_bt_err <- nn_gf[, c("NN_MSE_S", "Biotype_general")]

####
pRF_gf_trimmed <- pRF_gf[-na_hits, c("gRF_Cor_test", features)]
pRF_gf_trimmed_err <- pRF_gf[-na_hits, c("gRF_MSEs_back_test", features)]

pRF_gf_bt <- pRF_gf[, c("gRF_Cor_test", "Biotype_general")]
pRF_gf_bt <- cbind(pRF_gf_bt, rep("CRE-RF", nrow(pRF_gf_bt)))
colnames(pRF_gf_bt) <- c("cor","Biotype_general", "method")
pRF_gf_bt_err <- pRF_gf[, c("gRF_MSEs_back_test", "Biotype_general")]

####
bRF_gf_trimmed <- bRF_gf[-na_hits, c("ugRF_test_cor", features)]
bRF_gf_trimmed_err <- bRF_gf[-na_hits, c("ugRF_test_err", features)]

bRF_gf_bt <- bRF_gf[, c("ugRF_test_cor", "Biotype_general")]
bRF_gf_bt <- cbind(bRF_gf_bt, rep("Binned-RF", nrow(bRF_gf_bt)))
colnames(bRF_gf_bt) <- c("cor","Biotype_general", "method")
bRF_gf_bt_err <- bRF_gf[, c("ugRF_test_err", "Biotype_general")]

all_gf_bt_list <- list(cnn_gf_bt, nn_gf_bt, stitchit_gf_bt, pRF_gf_bt, bRF_gf_bt)
all_gf_bt_df  <- do.call("rbind", all_gf_bt_list)
all_gf_bt_df_filtered <- all_gf_bt_df[-(which(is.na(all_gf_bt_df$Biotype_general))),]
all_gf_bt_df_filtered <- all_gf_bt_df_filtered[-(which(all_gf_bt_df_filtered$Biotype_general=="TEC")),]

all_gf_bt_df_filtered$method <- factor(all_gf_bt_df_filtered$method, levels= names(colors))

bt_palette <- c("#04cbe3","#0c944c", "#0c6cfb")
names(bt_palette) <- c("protein_coding", "non-coding RNA", "pseudogene")

########################## Suppl. Fig. 4b ################################
pdf("~/model_performance_biotype_boxlot_pearson.pdf")
print(ggplot(all_gf_bt_df_filtered, aes(x= method, y= cor, fill= Biotype_general)) + geom_boxplot(draw_quantiles = c(0.5), 
        outlier.shape= NA) + scale_fill_manual(values = bt_palette) + 
        scale_y_continuous(limits = c(max(quantile(all_gf_bt_df_filtered$cor,0.25)-1.5*IQR(all_gf_bt_df_filtered$cor),min(all_gf_bt_df$cor)), 
        min(quantile(all_gf_bt_df_filtered$cor,0.75)+1.5*IQR(all_gf_bt_df_filtered$cor), max(all_gf_bt_df_filtered$cor)))) + ggtitle("Peformance stratified by biotype") + theme(panel.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),legend.key.size = unit(4.5, "lines"), legend.key = element_rect(fill = "white", color = NA), 
        text=element_text(family= font_config$family, size=7)))

dev.off()

### partial correlation introduction fig. ###
# 1. variance genomic features
# 2. scatter-plot genomic feature and performance

all_gf_trimmed <- c(CER_NN= list(nn_gf_trimmed), CER_RF= list(pRF_gf_trimmed), binned_RF= list(bRF_gf_trimmed), STITCHIT= list(stitchit_gf_trimmed), binned_CNN= list(cnn_gf_trimmed))

schabrackentapir.grau.dunkel <- "#596C89"
schabrackentapir.grau.hell <- "#D3D3D5"

x <- hist(as.numeric(all_gf_trimmed$CER_NN$"X.ENCODE_peCREs"), breaks=20)
cre_mlp_cres_df <- data.frame(
  x = rep(x$breaks, each=2), 
  y = c(0, rep(x$counts, each = 2), 0))

########################## Main Fig. 3a ################################
pdf("/projects/apog/work/CNN_bin/NN_STITCHIT_results/updated_plots/num_CREs_histogram_light.pdf", width=3.5, height=3)
print(ggplot(cre_mlp_cres_df, aes(x,y)) +
        geom_polygon(fill=schabrackentapir.grau.hell) + geom_line(col=schabrackentapir.grau.dunkel) + labs(x = "Number CREs", y = "Count") + theme_classic() + theme(text=element_text(family= font_config$family, size= font_config$size)))
dev.off()        

########################## Main Fig. 3b ################################
cor_CRE_RF_numbCREs <- cor(all_gf_trimmed$CER_RF[,1], all_gf_trimmed$CER_RF[,12], method = "pearson")
pdf("~/CRE_RF_pearson_numbCREs_scatter_corr_loess_colors_adapted_size.pdf", width=3.5, height=3)
ggplot(all_gf_trimmed$CER_RF, aes(X.ENCODE_peCREs, gRF_Cor_test)) +
  geom_point(size = 0.5, color = schabrackentapir.grau.hell, alpha = 0.3) + geom_smooth(method=loess, color="black", fill=schabrackentapir.grau.dunkel, span=0.4) +
  labs(
    x = "Number of CREs",
    y = "CRE-RF Pearson"
  ) +
  theme_classic()
dev.off()

#########################
cnn_pcor_sup <- lapply(seq(length(com_mat_cor$onevRest_res_genes)), function(i)pcor(all_gf_trimmed[[i]][which(all_gf_trimmed[[i]]$gene %in% cor_mat[com_mat_cor$onevRest_res_genes[[i]], ]$gene), -2], method= "spearman"))

cnn_pcor_sup_mat <- NULL
cnn_pcor_sup_pval_mat <- NULL
for(i in seq(length(cnn_pcor_sup))){
        cnn_pcor_sup_mat <- rbind(cnn_pcor_sup_mat, cnn_pcor_sup[[i]]$estimate[1, -1])
        cnn_pcor_sup_pval_mat <- rbind(cnn_pcor_sup_pval_mat, cnn_pcor_sup[[i]]$p.value[1, -1])
}

pval_cutoff <- .05
cnn_pcor_sup_mat_sig <- cnn_pcor_sup_mat
cnn_pcor_sup_mat_sig[cnn_pcor_sup_pval_mat > pval_cutoff] <- 0

rownames(cnn_pcor_sup_mat) <- rownames(cnn_pcor_sup_pval_mat) <- names(cnn_pcor_sup)
gf_strat_sup <- lapply(com_mat_cor$onevRest_res_genes, function(l) return(list(superior= gene_feat[gene_feat$Ensembl.ID %in% cor_mat[l, ]$gene,], rest= gene_feat[-which((gene_feat$Ensembl.ID %in% cor_mat[l, ]$gene)), ])))

df_all <- list()
for(ft in features[-1]){
        df <- NULL
        for(i in seq(length(gf_strat_sup))){
                df <- rbind(df, data.frame(model= names(gf_strat_sup)[i], value= gf_strat_sup[[i]]$superior[, colnames(gf_strat_sup[[i]]$superior) == ft], gene_set= "superior", genomic_feature= ft))
                df <- rbind(df, data.frame(model= names(gf_strat_sup)[i], value= gf_strat_sup[[i]]$rest[, colnames(gf_strat_sup[[i]]$rest) == ft], gene_set= "rest", genomic_feature= ft))
        }
        df_all[[ft]] <- df
}
tmp <- do.call(rbind, df_all)
tmp$model[which(tmp$model == "CNN_test_cor")] <- "Binned-CNN"
tmp$model[which(tmp$model == "ugRF_test_cor")] <- "Binned-RF"
tmp$model[which(tmp$model == "STITCHIT_cor")] <- "STITCHIT"
tmp$model[which(tmp$model == "gRF_Cor_test")] <- "CRE-RF"
tmp$model[which(tmp$model == "NN_Cor_test")] <- "CRE-MLP"
tmp$model <- factor(tmp$model, levels= names(colors))

selected_feats <- c(1, 6, 10)
df_all_selected <- do.call(rbind, df_all[selected_feats])

### obtain medians from boxplots (excluding outliers) ###
# for nonzero fraction and #ENCODE CREs
for (feature in c("X.ENCODE_peCREs", "Fraction_nonzero_in_test", "longest_5prime_UTR")){
        for (method in names(gf_strat_sup)){
        print(feature)
        print(method)
        method_df <-  df_all[[feature]][which(df_all[[feature]]$model==method),]
        method_df_trimmed <- method_df[which(method_df$value <= min(quantile(method_df$value,0.75)+1.5*IQR(method_df$value)) & method_df$value >= max(quantile(method_df$value,0.25)-1.5*IQR(method_df$value),min(method_df$value))), ]
        method_df_superior <- method_df_trimmed[which(method_df_trimmed$gene_set == "superior"), ]
        method_df_rest <- method_df_trimmed[which(method_df_trimmed$gene_set == "rest"), ]
        print(median(method_df_superior$value))
        print(median(method_df_rest$value))
        }
}

portogiesische.galeere.lila <- "#9E42A1" 
kea.grün <-  "#89C934"
dumbo.squid <- "#b27d74"

########################## Main Fig. 3d-e, Suppl. Fig. 4d-i ################################
pdf("~/superior_genes_gene_features_boxplot.pdf")

names(df_all)[which(names(df_all) == "Fraction_nonzero_in_test")] <- "Nonzero fraction"
names(df_all)[which(names(df_all) == "X.ENCODE_peCREs")] <- "Number ENCODE CREs"
names(df_all)[which(names(df_all) == "longest_5prime_UTR")] <- "Longest 5' UTR"

#flip to zero fraction
df_all[["Nonzero fraction"]][["value"]] <- 1 - df_all[["Nonzero fraction"]][["value"]]

p_val_list <- vector("list", length(df_all))
#calculate significance -> unpaired wilkoxon-rank test (Mann-Whitney U)
run_wilcox_unpaired <- function(df) {
        df %>%
        group_by(model) %>%
        summarise(
        genomic_feature = unique(genomic_feature),  
        p_value = wilcox.test(value[gene_set == "superior"], value[gene_set == "rest"])$p.value,
        significance = ifelse(p_value < 0.05, "s", "ns"))        
}

for(i in seq(length(df_all))){ #for each genomic feature 
        df_all[[i]]$model[which(df_all[[i]]$model == "CNN_test_cor")] <- "Binned-CNN"
        df_all[[i]]$model[which(df_all[[i]]$model == "ugRF_test_cor")] <- "Binned-RF"
        df_all[[i]]$model[which(df_all[[i]]$model == "STITCHIT_cor")] <- "STITCHIT"
        df_all[[i]]$model[which(df_all[[i]]$model == "gRF_Cor_test")] <- "CRE-RF"
        df_all[[i]]$model[which(df_all[[i]]$model == "NN_Cor_test")] <- "CRE-MLP"
        df_all[[i]]$model <- factor(df_all[[i]]$model, levels= names(colors))

        #calculate significance
        p_val_list[[i]] <- run_wilcox_unpaired(df_all[[i]])
        
	print(ggplot(df_all[[i]], aes(x= model, y= value, fill= gene_set)) + geom_boxplot(draw_quantiles = c(0.5), outlier.shape= NA) + scale_fill_manual(values=c(kea.grün, dumbo.squid)) + #scale_fill_brewer(palette="Accent") + 
                scale_y_continuous(limits = c(max(quantile(df_all[[i]]$value,0.25)-1.5*IQR(df_all[[i]]$value),min(df_all[[i]]$value)), min(quantile(df_all[[i]]$value,0.75)+1.5*IQR(df_all[[i]]$value), max(df_all[[i]]$value)))) + ggtitle(names(df_all)[i]) + theme(panel.background = element_rect(fill = "white", color = NA), panel.border = element_rect(color = "black", fill = NA, size = 1),legend.key.size = unit(4.5, "lines"), legend.key = element_rect(fill = "white", color = NA), text=element_text(family= font_config$family, size=7)))#, legend.text = element_text(size = 5)))
}
dev.off()

#save p-values to file
pval_df <- do.call(rbind, p_val_list)
write.table(pval_df, "~/unpaired_mann_whitney_u_pvals_superior_remaining_per_method.tsv", sep = "\t", col.names=TRUE, row.names = FALSE, quote = FALSE)

gene_feat_final <- gene_feat[hits,c("X.TSS", "Gene_length", "Exons_length","longest_5prime_UTR", "longest_3prime_UTR" , "Expression_ubiquitousness", "Fraction_nonzero_in_test", "Gene_density", "X.ENCODE_peCREs"),]

na_hits <- which(is.na(gene_feat_final$X.ENCODE_peCREs) == T)
if(length(na_hits)){
        gene_feat_final <- gene_feat_final[-na_hits, ]
}

#number test samples = 0 -> flip feature values 
gene_feat_final[["Fraction_nonzero_in_test"]]  <- 1 - gene_feat_final[["Fraction_nonzero_in_test"]] 

gene_feat_pearson_df <- cor(gene_feat_final)
gene_feat_spearman_df <- cor(gene_feat_final, method = "spearman")

library(reshape2)
gene_feat_pearson_long <- melt(gene_feat_pearson_df)
gene_feat_spearman_long <- melt(gene_feat_spearman_df)
colnames(gene_feat_pearson_long) <- c("feat1","feat2","pearson")
colnames(gene_feat_spearman_long) <- c("feat1","feat2","spearman")

target <- c("X.TSS", "Gene_length", "Exons_length", "longest_5prime_UTR", "longest_3prime_UTR", "Expression_ubiquitousness", "Fraction_nonzero_in_test", "Gene_density", "X.ENCODE_peCREs")
gene_feat_pearson_long$feat1 <- factor(gene_feat_pearson_long$feat1, levels = rev(target))
gene_feat_spearman_long$feat1 <- factor(gene_feat_spearman_long$feat1, levels = rev(target))

########################## Suppl. Fig. 4a ################################
pdf("~/genomic_features_cors_heatmap.pdf", height=7.5, width=7)
print(ggplot(gene_feat_pearson_long, aes(x= feat2, y= feat1, fill= pearson)) + geom_tile(color= "black") + scale_fill_gradient2(low = "#075AFF", mid = "#FFFFFF", high = "#FF0000") + geom_text(aes(label= round(pearson, 3)), color= "black", size= 4) + theme(panel.background = element_blank(), legend.position="top", axis.text.x= element_text(angle= 45, hjust= 1, vjust= 1))) #+ ggtitle("pcor between model test correlation and genomic features"))
print(ggplot(gene_feat_spearman_long, aes(x= feat2, y= feat1, fill= spearman)) + geom_tile(color= "black") + scale_fill_gradient2(low = "#075AFF", mid = "#FFFFFF", high = "#FF0000") + geom_text(aes(label= round(spearman, 3)), color= "black", size= 4) + theme(panel.background = element_blank(), legend.position="top", axis.text.x= element_text(angle= 45, hjust= 1, vjust= 1))) #+ ggtitle("pcor between model test correlation and genomic features"))
dev.off()


library(psych)
## convert Fraction_nonzero_in_test to Fraction_zero_in_test
cnn_gf_trimmed[, "Fraction_nonzero_in_test"] <- 1 - cnn_gf_trimmed[, "Fraction_nonzero_in_test"]
bRF_gf_trimmed[, "Fraction_nonzero_in_test"] <- 1 - bRF_gf_trimmed[, "Fraction_nonzero_in_test"]
pRF_gf_trimmed[, "Fraction_nonzero_in_test"] <- 1 - pRF_gf_trimmed[, "Fraction_nonzero_in_test"]
nn_gf_trimmed[, "Fraction_nonzero_in_test"] <- 1 - nn_gf_trimmed[, "Fraction_nonzero_in_test"]
stitchit_gf_trimmed[, "Fraction_nonzero_in_test"] <- 1 - stitchit_gf_trimmed[, "Fraction_nonzero_in_test"]

## convert expression breadth to 1- breadth
for(ft in structural_features){print(ft); print(partial.r(cnn_gf_trimmed, c("CNN_test_cor", ft), c(genomic_features1, genomic_features2)))}
get_pcor <- function(data, model_col, model, ft, ft1, ft2){
	pcor <-  partial.r(data, c(model_col, ft), c(ft1, ft2), method= "spearman")
	pcor_p <-  corr.p(pcor, nrow(data) - length(c(ft1, ft2)))
	df_res <- data.frame(pcor= pcor[1, 2], method= "spearman", feature= ft, model= model, pvalue= pcor_p$p[1, 2])
	return(df_res)
}
pcor_multivar <- NULL;
for(ft in structural_features){
	print(ft);
	pcor_multivar <- rbind(pcor_multivar, get_pcor(data= cnn_gf_trimmed, model_col= "CNN_test_cor", model= "Binned-CNN", ft= ft, ft1= genomic_features1, ft2= genomic_features2))
	pcor_multivar <- rbind(pcor_multivar, get_pcor(data= bRF_gf_trimmed, model_col= "ugRF_test_cor", model= "Binned-RF", ft= ft, ft1= genomic_features1, ft2= genomic_features2))
	pcor_multivar <- rbind(pcor_multivar, get_pcor(data= pRF_gf_trimmed, model_col= "gRF_Cor_test", model= "CRE-RF", ft= ft, ft1= genomic_features1, ft2= genomic_features2))
	pcor_multivar <- rbind(pcor_multivar, get_pcor(data= nn_gf_trimmed, model_col= "NN_Cor_test", model= "CRE-MLP", ft= ft, ft1= genomic_features1, ft2= genomic_features2))
	pcor_multivar <- rbind(pcor_multivar, get_pcor(data= stitchit_gf_trimmed, model_col= "STITCHIT_cor", model= "STITCHIT", ft= ft, ft1= genomic_features1, ft2= genomic_features2))
}
pcor_multivar_g1 <- NULL;
for(ft in genomic_features1){
	print(ft);
	pcor_multivar_g1 <- rbind(pcor_multivar_g1, get_pcor(data= cnn_gf_trimmed, model_col= "CNN_test_cor", model= "Binned-CNN", ft= ft, ft1= structural_features, ft2= genomic_features2))
	pcor_multivar_g1 <- rbind(pcor_multivar_g1, get_pcor(data= bRF_gf_trimmed, model_col= "ugRF_test_cor", model= "Binned-RF", ft= ft, ft1= structural_features, ft2= genomic_features2))
	pcor_multivar_g1 <- rbind(pcor_multivar_g1, get_pcor(data= pRF_gf_trimmed, model_col= "gRF_Cor_test", model= "CRE-RF", ft= ft, ft1= structural_features, ft2= genomic_features2))
	pcor_multivar_g1 <- rbind(pcor_multivar_g1, get_pcor(data= nn_gf_trimmed, model_col= "NN_Cor_test", model= "CRE-MLP", ft= ft, ft1= structural_features, ft2= genomic_features2))
	pcor_multivar_g1 <- rbind(pcor_multivar_g1, get_pcor(data= stitchit_gf_trimmed, model_col= "STITCHIT_cor", model= "STITCHIT", ft= ft, ft1= structural_features, ft2= genomic_features2))
}
pcor_multivar_g2 <- NULL;
for(ft in genomic_features2){
	print(ft);
	pcor_multivar_g2 <- rbind(pcor_multivar_g2, get_pcor(data= cnn_gf_trimmed, model_col= "CNN_test_cor", model= "Binned-CNN", ft= ft, ft1= structural_features, ft2= genomic_features1))
	pcor_multivar_g2 <- rbind(pcor_multivar_g2, get_pcor(data= bRF_gf_trimmed, model_col= "ugRF_test_cor", model= "Binned-RF", ft= ft, ft1= structural_features, ft2= genomic_features1))
	pcor_multivar_g2 <- rbind(pcor_multivar_g2, get_pcor(data= pRF_gf_trimmed, model_col= "gRF_Cor_test", model= "CRE-RF", ft= ft, ft1= structural_features, ft2= genomic_features1))
	pcor_multivar_g2 <- rbind(pcor_multivar_g2, get_pcor(data= nn_gf_trimmed, model_col= "NN_Cor_test", model= "CRE-MLP", ft= ft, ft1= structural_features, ft2= genomic_features1))
	pcor_multivar_g2 <- rbind(pcor_multivar_g2, get_pcor(data= stitchit_gf_trimmed, model_col= "STITCHIT_cor", model= "STITCHIT", ft= ft, ft1= structural_features, ft2= genomic_features1))
}

pcor_multivar_all <- rbind(pcor_multivar, pcor_multivar_g1, pcor_multivar_g2)
pcor_multivar_all$feature[pcor_multivar_all$feature == "Fraction_nonzero_in_test"] <- "Fraction_zero_in_test"

target <- c("Number of TSSs", "Gene length", "Exon length", "longest 5' UTR", "longest 3' UTR", "Expression breadth", "Test sample expr = 0 (%)", "Number of genes", "Number of CREs")

## rename features for plotting
pcor_multivar_all$feature[pcor_multivar_all$feature == "X.TSS"] <- "Number of TSSs"
pcor_multivar_all$feature[pcor_multivar_all$feature == "Gene_length"] <- "Gene length"
pcor_multivar_all$feature[pcor_multivar_all$feature == "Exons_length"] <- "Exon length"
pcor_multivar_all$feature[pcor_multivar_all$feature == "longest_5prime_UTR"] <- "longest 5' UTR"
pcor_multivar_all$feature[pcor_multivar_all$feature == "longest_3prime_UTR"] <- "longest 3' UTR"
pcor_multivar_all$feature[pcor_multivar_all$feature == "Fraction_zero_in_test"] <- "Test sample expr = 0 (%)"
pcor_multivar_all$feature[pcor_multivar_all$feature == "Gene_density"] <- "Number of genes"
pcor_multivar_all$feature[pcor_multivar_all$feature == "Expression_ubiquitousness"] <- "Expression breadth"
pcor_multivar_all$feature[pcor_multivar_all$feature == "X.ENCODE_peCREs"] <- "Number of CREs"

pcor_multivar_all$feature <- factor(pcor_multivar_all$feature, levels = target)
pcor_multivar_all <- pcor_multivar_all[order(pcor_multivar_all$feature), ]

pcor_multivar_all$model <- factor(pcor_multivar_all$model, levels= rev(names(colors)))

plush.crested.jay.blue <- "#2572E4"
mantis.shrimp.red <- "#D53610"

########################## Main Fig. 1c ################################
pdf("~/partial_cor_multi_feature_categories.pdf")
print(ggplot(pcor_multivar_all, aes(x= feature, y= model, fill= pcor), color= "black") + geom_tile(color= "black") + scale_fill_gradient2(low = plush.crested.jay.blue, mid = "#FFFFFF", high = mantis.shrimp.red) + geom_text(aes(label= round(pcor, 2)), color= "black", size= 4) + theme(panel.background = element_blank(), legend.position="top", axis.text.x= element_text(color= "black", angle= 45, hjust= 1, vjust= 1), axis.text.y= element_text(color= "black"))) 
print(ggplot(pcor_multivar, aes(x= feature, y= model, fill= pcor)) + geom_tile(color= "black") + scale_fill_gradient2(low = plush.crested.jay.blue, mid = "#FFFFFF", high = mantis.shrimp.red) + geom_text(aes(label= round(pcor, 3)), color= "black", size= 4) + theme(legend.position="top", axis.text.x= element_text(angle= 45, hjust= 1, vjust= 1)) + ggtitle("pcor(model_cor, structral_feat | genomic_feat)"))
print(ggplot(pcor_multivar_g1, aes(x= feature, y= model, fill= pcor)) + geom_tile(color= "black") + scale_fill_gradient2(low = plush.crested.jay.blue, mid = "#FFFFFF", high = mantis.shrimp.red ) + geom_text(aes(label= round(pcor, 3)), color= "black", size= 4) + theme(legend.position="top", axis.text.x= element_text(angle= 45, hjust= 1, vjust= 1)) + ggtitle("pcor(model_cor, genomic_feat1 | (structral_feat, genomic_feat2))"))
print(ggplot(pcor_multivar_g2, aes(x= feature, y= model, fill= pcor)) + geom_tile(color= "black") + scale_fill_gradient2(low = plush.crested.jay.blue, mid = "#FFFFFF", high = mantis.shrimp.red ) + geom_text(aes(label= round(pcor, 3)), color= "black", size= 4) + theme(legend.position="top", axis.text.x= element_text(angle= 45, hjust= 1, vjust= 1)) + ggtitle("pcor(model_cor, genomic_feat2 | (structral_feat, genomic_feat1))"))
dev.off()

########################## Suppl. Fig. 1a: CIRCOS PLOT ################################

library(circlize)
library(viridis)

samples <- read.csv("~/IHEC_metadata_harmonization.v1.2.extended.csv")
data <- read.table("~/ENSG00000150687_w1MB_100bs_BinnedActivity.txt.gz", header= T)
hit_idx <- which(samples$epirr_id_without_version %in% data$Sample)
samples <- samples[hit_idx, ]
data_gist <- data[, c(1, 10002)]
data_sample <- merge(data_gist, samples, by.x= "Sample", by.y= "epirr_id_without_version")

organ <- data_sample$harmonized_sample_organ_system_order_AnetaMikulasova
ontology <-  data_sample$harmonized_sample_ontology_intermediate
cnt <- table(ontology)
cnt <- sort(cnt)
sectors <- names(cnt)

s1 <- factor(sectors, levels= sectors)
ranges <- sapply(cnt, function(i) c(0, i))

mycol <- c(viridis(length(unique(organ)) / 2 + 2, begin= 0)[seq(6)], magma(length(unique(organ)) / 2+2, begin= 0)[seq(6)])
organ_col <- as.matrix(mycol);
rownames(organ_col)= unique(organ);
mycol_all <- vector("character", length(sectors));
for(s in seq(length(sectors))){
        hits <- which(ontology == sectors[s]);
        mycol_all[s] <- organ_col[which(rownames(organ_col) == organ[hits][1])]
}

organ_col_df <- data.frame(organ_col)

pdf("~/sample_ontology_circlized_legend.pdf")

circos.clear(); circos.par$circle.margin <- 1.9; circos.par$gap.after <- 6; circos.par(cell.padding = c(0.00, 0.0,0.0, 0.00)); circos.initialize(sectors= s1, xlim= t(ranges)); circos.trackPlotRegion(factors = s1, y= seq(length(s1)), panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, cex= .9, adj = c(0, 0.5), col= mycol_all[CELL_META$sector.numeric.index]); circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 20, col = mycol_all[CELL_META$sector.numeric.index])}, bg.col= NA, track.height= .1, bg.border= NA)
legend("topright", legend = rownames(organ_col_df), fill = organ_col_df[,1], cex = .5, bty = "n", )

dev.off()
