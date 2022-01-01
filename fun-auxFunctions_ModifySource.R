## ==================================================================================== ##
# transcriptAchilles Shiny App for predicting transcript biomarkers.
# Copyright (C) 2018 Fernando Carazo
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# You may contact the author of this code, Fernando Carazo, at <fcarazo@tecnun.es>
## ==================================================================================== ##

### Aux Functions ###


# ggplot of Essentiality score of a gene (gn) against expression value of a transcript (tr)
plotAtaris <- function (tr, gn, iCL, getBM, CCLE_ataris, CCLE_tpm_log2, CCLE_sample_info, method = "pearson", n = NULL, text = TRUE, subtype = FALSE, hline = 0){
  
  CCLE_tpm_log2 <- as.matrix(CCLE_tpm_log2)
  CCLE_ataris <- as.matrix(CCLE_ataris)
  
  A <- data.frame(log2_TPM = CCLE_tpm_log2[tr, iCL], 
                  Essentiality = CCLE_ataris[gn, iCL], 
                  Primary.Site = CCLE_sample_info$Primary.Site[iCL],
                  Subtype = CCLE_sample_info$Subtype[iCL])
  
  pv <- cor.test(CCLE_tpm_log2[tr, iCL], CCLE_ataris[gn, iCL],method = method)$p.value
  cor <- cor(CCLE_tpm_log2[tr, iCL], CCLE_ataris[gn, iCL],method = method)
  
  info_t <- getBM[which(getBM$Transcript_ID==tr),]
  
  
  
  if(subtype){
    a <- (ggplot(A, aes(log2_TPM, Essentiality, color = Subtype)) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            geom_point(size =3) +
            ggtitle(sprintf("Essentiality: %s  |  Biomarker: %s", gn, info_t$Transcrip_name)) +
            geom_text(aes(-Inf, Inf, label= sprintf("P.v = %.2e \n cor = %.2f", pv, cor)),
                      vjust = "inward",
                      hjust = "inward", color =2, size = 3.5))
  }else{
    a <- (ggplot(A, aes(log2_TPM, Essentiality, color = Primary.Site)) + 
            geom_point(size =3) +
            theme(plot.title = element_text(hjust = 0.5)) +
            ggtitle(sprintf("Essentiality: %s  |  Biomarker: %s", gn, info_t$Transcrip_name)) +
            geom_text(aes(-Inf, Inf, label= sprintf("P.v = %.2e \n cor = %.2f", pv, cor)),
                      vjust = "inward",
                      hjust = "inward", color =2, size = 3.5))
  }
  
  if(text){a <- a + geom_text(aes(label=Primary.Site, color = Primary.Site), size = 3, hjust=0, vjust=-1) }
  # print(a)
  
  if(hline != 0){
    a <- a + geom_hline(yintercept = hline, col = 1, linetype = 2)
  }
  
  return(a)
}



corPv_all <- function (CCLE_ataris, CCLE_tpm_log2, nCors = 1000, method = "pearson"){
  pb <- progress_bar$new(
    format = "  calculating [:bar] :percent progress: :eta",
    total = nCors, clear = FALSE, width= 80)
  
  cors <- cor(t(CCLE_ataris), t(CCLE_tpm_log2) , method = method)  
  cors <- data.frame(Essential_gn = rownames(cors), cors)
  cors <- gather(cors, Transcript, cors, -Essential_gn)
  ord2 <- order(abs(cors$cors), decreasing = T)
  cors <- cors[ord2,]
  rownames(cors) <- NULL
  
  pv <- rep(NA, nrow(cors))
  for (p in 1:nCors){
    pb$tick()
    TR <- cors$Transcript[p]
    GN <- cors$Essential_gn[p]
    pv[p] <- cor.test(CCLE_ataris[GN, ], CCLE_tpm_log2[TR, ], method = method)$p.value
  }
  
  cors$pv <- pv
  # cors$adj.pv <- p.adjust(pv, "BH")
  return(cors)
}

corPv <- function (CCLE_ataris, CCLE_tpm_log2, nCors = 1000, method = "pearson"){
  pb <- progress_bar$new(
    format = "  calculating [:bar] :percent progress: :eta",
    total = nCors, clear = FALSE, width= 80)
  
  cors <- cor(t(CCLE_ataris), t(CCLE_tpm_log2) , method = method)  
  cors <- data.frame(Essential_gn = rownames(cors), cors)
  cors <- gather(cors, Transcript, cors, -Essential_gn)
  ord2 <- order(abs(cors$cors), decreasing = T)
  cors <- cors[ord2[1:nCors],]
  rownames(cors) <- NULL
  
  pv <- rep(NA, nCors)
  for (p in 1:nCors){
    pb$tick()
    TR <- cors$Transcript[p]
    GN <- cors$Essential_gn[p]
    pv[p] <- cor.test(CCLE_ataris[GN, ], CCLE_tpm_log2[TR, ], method = method)$p.value
  }
  
  cors$pv <- pv
  cors$adj.pv <- p.adjust(pv, "BH")
  return(cors)
}


# OLD FUNCTION
# plotAtaris_transcripts <- function(tr, gn, getBM, iCL, method = "pearson", n=NULL, text = TRUE, subtype = FALSE){
#   
#   # first plot: correlations ###########################################
#   gg1 <- plotAtaris(tr = tr, 
#                   gn = gn, getBM = getBM,
#                   iCL = iCL,
#                   method = method, n, text, subtype)
#   
#   trN <- getBM$Transcrip_name[which(getBM$Transcript_ID == tr)]
#   
#   gene <- getBM$Gene_name[which(getBM$Transcrip_name == trN)]
#   Trs <- which(getBM$Gene_name == gene)
#   ataris_gn <- gn
#   
#   ord <- order(CCLE_ataris[ataris_gn, iCL], decreasing = T)
#   a <- data.frame(CCLE_tpm_log2[Trs,iCL])
#   
#   if((ncol(a)==1)){
#     a <- data.frame(a[ord,])
#     rownames(a) <- make.unique(CCLE_sample_info$Primary.Site[iCL])
#     colnames(a) <- getBM$Transcrip_name[Trs]
#     gg3 <- ggmatplot(a, highlight = c(gene, trN)) + 
#       ylab("log2_tpm") + 
#       xlab("Samples") + 
#       ggtitle(sprintf("%s: Transcript sel: %s", gene, trN))
#     
#   }
#   
#   if((ncol(a)>1)){
#     
#     colnames(a) <- make.unique(CCLE_sample_info$Primary.Site[iCL])
#     rownames(a) <- getBM$Transcrip_name[Trs]
#     gnE <- log2(colSums(2^a-1)+1)
#     a <- rbind(gnE,a)
#     rownames(a)[1] <- gene
#     
#     gg3 <- ggmatplot(t(a[ord]), highlight = c(gene, trN)) + 
#       ylab("log2_tpm") + 
#       xlab("Samples") + 
#       ggtitle(sprintf("%s: Transcript sel: %s", gene, trN))
#     
#   }
#   
#   b <- data.frame(CCLE_ataris[ataris_gn, iCL])
#   colnames(b) <- "Ataris"
#   rownames(b) <- make.unique(CCLE_sample_info$Primary.Site[iCL])
#   b$type <- "ATARiS"
#   b <- b[ord,]
#   b$ix <- 1:nrow(b)
#   
#   gg2 <- ggplot(data = b, aes(ix, Ataris)) + 
#     geom_line(aes(color = type), size = 1, linetype =2) + 
#     ggtitle(sprintf("ATARiS Score gene: %s", ataris_gn)) +
#     theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#   
#   return(list(gg1 = gg1, gg2 = gg2, gg3 = gg3))
# 
# }



plot_transcripts_essential2 <- function(tr, gn, getBM, CCLE_ataris, 
                                        CCLE_tpm_log2, CCLE_sample_info,
                                        iCL, method = "pearson", group_bmrk = "transcript",
                                        n=NULL, text = TRUE, subtype = FALSE, vline = -2, xnames = TRUE){
  hline <- vline
  CCLE_ataris <- as.matrix(CCLE_ataris)
  CCLE_tpm_log2 <- as.matrix(CCLE_tpm_log2)
  
  
  
  # Extract all transcripts of a specific gene
  
  
  trN <- getBM$Transcrip_name[which(getBM$Transcript_ID == tr)]
  gene <- getBM$Gene_name[which(getBM$Transcrip_name == trN)]
  
  if(group_bmrk == "gene"){
    trN <- gene
  }
  
  Trs <- which(getBM$Gene_name == gene)
  Essential_gn <- gn
  ord <- order(CCLE_ataris[Essential_gn, iCL], decreasing = T)
  a <- data.frame(CCLE_tpm_log2[Trs,iCL])
  
  if((ncol(a)==1)){
    gg1 <- plotAtaris(tr = tr, 
                      gn = gn, getBM = getBM, CCLE_tpm_log2 = CCLE_tpm_log2,
                      CCLE_sample_info = CCLE_sample_info,
                      CCLE_ataris = CCLE_ataris,
                      iCL = iCL, hline = hline,
                      method = method, n, text, subtype) +
      theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5)) +
      ggtitle(sprintf("Essentiality: %s  |  Biomarker: %s", gn, trN), subtitle = sprintf("(Biomarker type: %s)", group_bmrk))
    
    a <- data.frame(a[ord,])
    rownames(a) <- make.unique(CCLE_sample_info$CCLE.name[iCL])
    colnames(a) <- getBM$Transcrip_name[Trs]
    
    gg3 <- ggmatplot(a, highlight = c(gene, trN), xnames = xnames) + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      ylab("log2_tpm") + 
      xlab("Samples") + 
      ggtitle("log2 Expression", subtitle = sprintf("Gene %s  |  Biomarker: %s", gene, trN))
    
    gg4 <- NULL    
    
  }
  
  if((ncol(a)>1)){
    
    colnames(a) <- make.unique(CCLE_sample_info$CCLE.name[iCL])
    rownames(a) <- getBM$Transcrip_name[Trs]
    gnE <- log2(colSums(2^a-1)+1)
    auxE <- 2^a-1
    psi <- data.frame(t(t(auxE) / colSums(auxE)))
    a <- rbind(gnE,a)
    rownames(a)[1] <- gene
    
    # first plot: correlations ###########################################
    
    if(group_bmrk == "gene"){
      gg1 <- plotAtaris(tr = gene, 
                        gn = gn, getBM = getBM, CCLE_tpm_log2 = a,
                        CCLE_sample_info = CCLE_sample_info[iCL, ],
                        CCLE_ataris = CCLE_ataris[, iCL],
                        iCL = 1:ncol(a), hline = hline,
                        method = method, n, text, subtype) +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        ggtitle(sprintf("Essentiality: %s  |  Biomarker: %s", gn, trN), subtitle = sprintf("(Biomarker type: %s)", group_bmrk))
    }else{
      gg1 <- plotAtaris(tr = tr, 
                        gn = gn, getBM = getBM, CCLE_tpm_log2 = CCLE_tpm_log2,
                        CCLE_sample_info = CCLE_sample_info,
                        CCLE_ataris = CCLE_ataris,
                        iCL = iCL, hline = hline,
                        method = method, n, text, subtype) +
        theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        ggtitle(sprintf("Essentiality: %s  |  Biomarker: %s", gn, trN), subtitle = sprintf("(Biomarker type: %s)", group_bmrk))
      
    }
    
    gg3 <- ggmatplot(t(a[ord]), highlight = c(gene, trN), xnames = xnames) + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      ylab("log2_tpm") + 
      xlab("Samples") + 
      ggtitle("log2 Expression", subtitle = sprintf("Gene %s  |  Biomarker: %s", gene, trN))
    
    gg4 <- ggmatplot(t(psi[ord]), highlight = c(trN), xnames = xnames) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      ylab("PSI") + 
      xlab("Samples") + 
      ggtitle(sprintf("Gene %s: Marker: %s", gene, trN))
    
  }
  
  b <- data.frame(Essentiality = CCLE_ataris[Essential_gn, iCL])
  rownames(b) <- make.unique(CCLE_sample_info$CCLE.name[iCL])
  b$type <- "Essentiality"
  b <- b[ord,]
  b$ix <- 1:nrow(b)
  
  gg2 <- ggplot(data = b, aes(ix, Essentiality)) + 
    geom_line(aes(color = type), size = 1, linetype =2) + 
    ggtitle(sprintf("Essentiality of %s", Essential_gn)) +
    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_blank())
  
  
  if(vline != 0){
    xint <- which(sort(CCLE_ataris[gn, iCL], decreasing = T) <= vline)[1] - 0.5
    gg2 <- gg2 + geom_vline(xintercept = xint, col = 1, linetype=2)
    gg2 <- gg2 + geom_hline(yintercept = vline, col = 1, linetype=2)
    gg3 <- gg3 + geom_vline(xintercept = xint, col = 1, linetype=2)
    if(ncol(a)>1) gg4 <- gg4 + geom_vline(xintercept = xint, col = 1, linetype=2)
  }
  
  # ROC
  
  test <- data.frame(Predictor = as.numeric(a[trN, ]), Label = as.numeric((CCLE_ataris[gn, iCL] < vline)+0))
  
  if((ncol(a)==1)){test <- data.frame(Predictor = as.numeric(a[,1]), Label = as.numeric((CCLE_ataris[gn, iCL] < vline)+0))}
  
  pred <- prediction(predictions = test$Predictor, labels = test$Label)
  auc1 <- unlist(performance(pred, measure = "auc")@y.values)
  
  pred <- prediction(predictions = 1- test$Predictor, labels = test$Label)
  auc2 <- unlist(performance(pred, measure = "auc")@y.values)
  
  if(auc1 > auc2){
    auc <- auc1
    gg5 <- ggplot(test, aes(d = Label, m = Predictor)) + 
      theme_minimal() +
      ggtitle("ROC Curve", subtitle = sprintf("Essentiality: %s  |  Biomarker: %s", gn, trN)) +
      theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5)) +
      geom_roc(labels = FALSE)+
      annotate("text", label = sprintf("AUC = %s", round(auc, 2)), x=0.6875, y=0.3125) + xlab("FPR") + ylab("TPR")
    
  }else{
    auc <- auc2
    gg5 <- ggplot(test, aes(d = Label, m = 1-Predictor)) + 
      ggtitle("ROC Curve", subtitle = sprintf("Essentiality: %s  |  Biomarker: %s", gn, trN)) +
      theme_minimal() +
      theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5)) +
      geom_roc(labels = FALSE) + 
      annotate("text", label = sprintf("AUC = %s", round(auc, 2)), x=0.6875, y=0.3125) + xlab("FPR") + ylab("TPR")
    
    
  }
  
  
  
  return(list(gg1 = gg1, gg2 = gg2, gg3 = gg3, gg4 = gg4, gg5 = gg5))
  
}


resnik_ataris <- function(TR, GN, GeneGO){
  g1 <- as.character(GN); g1 <- getBM$Gene_ID[grep(g1, getBM$Gene_name)] [1]
  g2 <- getBM$Gene_ID[which(TR == getBM$Transcript_ID)] [1]
  
  isGO1 <- which(rownames(GeneGO)==g1)
  isGO2 <- which(rownames(GeneGO)==g2)
  
  if(is.na(g1) | length(isGO1) == 0 | length(isGO2) == 0){resnik <- NA
  }else{resnik <- quasiresnikDistance12(gen1 = g1, gen2 = g2, GeneGO = GxGO) }
  
  return(resnik)
}





plot_transcripts <- function(tr_expression, tr, getBM, n= "", text = TRUE, subtype = FALSE, vline = NULL, xnames = TRUE){
  
  # first plot: correlations ###########################################
  
  trN <- getBM$Transcrip_name[which(getBM$Transcript_ID == tr)]
  
  gene <- getBM$Gene_name[which(getBM$Transcript_ID == tr)]
  Trs <- which(getBM$Gene_name == gene)
  
  a <- data.frame(tr_expression[Trs,])
  
  if((ncol(a)==1)){
    a <- data.frame(a)
    colnames(a) <- getBM$Transcrip_name[Trs]
    gg3 <- ggmatplot(a, highlight = c(gene, trN), xnames = xnames) + 
      ylab("log2_tpm") + 
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Samples") + 
      ggtitle(sprintf("%s Gene %s: Marker: %s",n, gene, trN))
    
  }
  
  if((ncol(a)>1)){
    
    # colnames(a) <- make.unique(CCLE_sample_info$Primary.Site[iCL])
    rownames(a) <- getBM$Transcrip_name[Trs]
    gnE <- log2(colSums(2^a-1)+1)
    a <- rbind(gnE,a)
    rownames(a)[1] <- gene
    
    gg3 <- ggmatplot(t(a), highlight = c(gene, trN), xnames = xnames) + 
      ylab("log2_tpm") + 
      xlab("Samples") + 
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(sprintf("%s Gene %s: Marker: %s",n, gene, trN))
    
  }
  
  if(!is.null(vline)){
    gg3 <- gg3 + geom_vline(xintercept = vline, col = 1, linetype = 2)
  }
  
  # b <- data.frame(CCLE_ataris[ataris_gn, iCL])
  # colnames(b) <- "Ataris"
  # rownames(b) <- make.unique(CCLE_sample_info$Primary.Site[iCL])
  # b$type <- "ATARiS"
  # b <- b[ord,]
  # b$ix <- 1:nrow(b)
  # 
  # gg2 <- ggplot(data = b, aes(ix, Ataris)) + 
  #   geom_line(aes(color = type), size = 1, linetype =2) + 
  #   ggtitle(sprintf("ATARiS Score gene: %s", ataris_gn)) +
  #   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  # 
  # return(list(gg1 = gg1, gg2 = gg2, gg3 = gg3))
  return(gg3)
  
}
