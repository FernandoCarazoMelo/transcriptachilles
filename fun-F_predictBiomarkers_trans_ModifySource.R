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



#inputs
# iCL, iGn_ess, Essentiality,Exp_transcripts,getBM,
#aim cube where the dimensions are x=Exp_transcripts,y=EG,z=characteristics


# parameters
# co_demeter: essentiality threshold,under-2 a gene is labeled essential


F_predictBiomarkers_trans <- function(Essentiality, Exp_trancripts, Exp_genes, iCL, iGn_ess, iT_exp, iG_exp, getBM,
                                      co_demeter=-2, impute = FALSE, restCL = TRUE, string = FALSE, shiny = TRUE){
  require(limma)
  require(STRINGdb)
  require(scales)
  require(Matrix)
  require(impute)
  require(progress)
  require(qvalue)
  
  # dicotomise Essentiality
  Essentiality_d <- (Essentiality < co_demeter) * 1
  
  Exp_trancripts_log2 <- as.data.frame(log2(1 + Exp_trancripts))
  
  # comment in shiny!!!
  getBM_ntr <- as.data.frame(getBM %>% group_by(Gene_ID) %>% dplyr::mutate(tr = n()))
  
  # Gene expression
  Exp_genes_log2 <- as.data.frame(log2(1 + Exp_genes))
  
  ##########################################################################################
  
  nmTot_Genes <- unique(c(getBM$Gene_name[iT_exp], rownames(Exp_genes)[iG_exp]))
  nTot_Genes <- length(nmTot_Genes)
  
  Tr_Bmrk_list <- vector(length = length(iGn_ess), mode = "list")
  names(Tr_Bmrk_list) <- rownames(Essentiality_d)[iGn_ess]
  
  pb <- progress_bar$new(
    format = "  calculating biomarkers [:bar] :percent progress: :eta",
    total = length(iGn_ess), clear = FALSE, width= 80)
  
  for(m in 1:length(iGn_ess)){ # length(iM)
    pb$tick()
    
    # Just SHINY
    # Increment the progress bar, and update the detail text.
    if(shiny) incProgress(1/length(iGn_ess), detail = sprintf("gene %s/%s", m, length(iGn_ess)))
    
    essentiality <- rownames(Essentiality_d)[iGn_ess[m]]
    
    Dmatrix <- cbind(1, Essentiality_d[essentiality, iCL])
    if(sum(Dmatrix[, 2]) == nrow(Dmatrix)) next()
    colnames(Dmatrix) <- c("Intercept", essentiality)
    
    Cmatrix <- t(t(c(0,1)))
    
    # limma Transcripts------------------------------------------------------
    
    Fit_trans<-lmFit(Exp_trancripts_log2[iT_exp, iCL], Dmatrix)
    Fit_trans<-contrasts.fit(Fit_trans, Cmatrix)
    Fit_trans <- eBayes(Fit_trans)
    topT <- topTable(Fit_trans, number = Inf)
    topT$lfdr <- lfdr(topT$P.Value)
    topT$Group_bmkr <- "transcript"
    
    # Merge data with transcript info
    iB_t <- match(rownames(topT), getBM_ntr$Transcript_ID)
    topT_c <- cbind(topT, getBM_ntr[iB_t,])
    rownames(topT_c) <- NULL
    
    # Initial correction of P-values
    topT_c$P.Value_corr_tr <- 1 - (1 - topT_c$P.Value) ^ topT_c$tr
    
    # Save Fit limma
    # Fit_trans_list[[m]] <- Fit_trans
    
    
    # Limma Genes----------------------------------------------------
    Fit_gn <- lmFit(Exp_genes_log2[iG_exp, iCL], Dmatrix)
    Fit_gn <- contrasts.fit(Fit_gn, Cmatrix)
    Fit_gn <- eBayes(Fit_gn)
    topG <- topTable(Fit_gn, number = Inf)
    topG$lfdr <- lfdr(topG$P.Value)
    topG$Group_bmkr <- "gene"
    
    # Merge data with transcript info
    iB_g <- match(rownames(topG), getBM_ntr$Gene_name)
    topG_c <- cbind(topG, getBM_ntr[iB_g,])
    # topG_c$Transcript_ID <- NA
    # topG_c$Transcrip_name <- NA
    rownames(topG_c) <- NULL
    topG_c$Group_bmkr[which(topG_c$tr == 1)] <- "gene_oneTr"
    
    
    topG_c$P.Value_corr_tr <- topG_c$P.Value
    
    # Combine GT---------------------------------------------------------
    
    # B_tr_c$P.Value_corr_tr <- 1 - (1 - B_tr_c$P.Value) ^ B_tr_c$tr
    B_GT <- rbind(topT_c, topG_c)
    B_GT <- data.frame(B_GT %>% group_by(Gene_name) %>% slice(which.min(lfdr)))
    # B_GT <- B_GT[order(B_GT$P.Value), ]
    rownames(B_GT) <- NULL
    # B_GT_print <- B_GT %>% dplyr::select(Gene_name, Transcrip_name, Transcript_ID, Biotype, tr, logFC, AveExpr, P.Value, P.Value_corr_tr, adj.P.Val, lfdr, Group_bmkr)
    B_GT_print <- B_GT %>% dplyr::select(Gene_name, Transcrip_name, Transcript_ID, Biotype, tr, logFC, AveExpr, P.Value, lfdr, Group_bmkr)
    
    iP <- match(nmTot_Genes, B_GT_print$Gene_name)
    B_GT_print <- B_GT_print[iP, ]
    colnames(B_GT_print)[1:2] <- c("Gene_bmkr", "Transcrip_bmkr")
    Tr_Bmrk_list[[m]] <- B_GT_print
  }
  return(list(Tr_Bmrk_list = Tr_Bmrk_list, qvalue_Gn = qvalue(topG$P.Value), qvalue_Tr = qvalue(topT$P.Value)))
  
}


# if(restCL){
#   Dmatrix_rest <- cbind(1,Essentiality_d[which(rownames(Essentiality_d) == essentiality),-iCL])
#   colnames(Dmatrix_rest) <- c("Intercept", essentiality)
#   if(sum(Dmatrix_rest[,2]) > 0){
#     # with the rest of CL
#     Fit_trans_rest<-lmFit(Exp_trancripts_log2[iT_exp, -iCL], Dmatrix_rest)
#     Fit_trans_rest<-contrasts.fit(Fit_trans_rest, Cmatrix)
#     Fit_trans_rest <- eBayes(Fit_trans_rest)
#     topT_rest <- toptable(Fit_trans_rest, number = Inf, sort.by = "none")
#     
#     # topT <- data.frame(Gene = rownames(topT), topT)
#     # View(topT)
#     MagicCube_Trans[ m, , "RestCL_t"] <- topT_rest$t
#     MagicCube_Trans[ m, , "RestCL_P.Value"] <- topT_rest$P.Value
#     MagicCube_Trans[ m, , "RestCL_log2FC"] <- topT_rest$logFC
#     MagicCube_Trans[ m, ,"Same_Dir"] <- ((topT$t * topT_rest$t) > 0) * 1
#   }else{
#     MagicCube_Trans[ m, , "RestCL_t"] <- NA
#     MagicCube_Trans[ m, , "RestCL_P.Value"] <- NA
#     MagicCube_Trans[ m, , "RestCL_log2FC"] <- NA
#     MagicCube_Trans[ m, ,"Same_Dir"] <- NA
#   }
# }

# Adjusting pvalues by transcript
# Pvs <- MagicCube_Trans[ , , "P.Value"]
# adj.P.Val_by_Trans <- apply(Pvs, 2, function(pv){pvA <- p.adjust(pv, method = "BH"); return(pvA)})
# MagicCube_Trans[ , ,  "adj.P.Val_by_Trans"] <- t(adj.P.Val_by_Trans)


## String interactions

# if(string){
#   #(9606 for Human, 10090 for mouse).
#   string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="")
#   string_proteins <- string_db$get_proteins()
#   
#   gns <- unique(c(rownames(Essentiality[iGn_ess,]), getBM$Gene_name[iT_exp]))
#   ids <- string_proteins$protein_external_id[string_proteins$preferred_name %in% gns]
#   
#   Genes_not_mapped <- gns[!(gns %in% string_proteins$preferred_name)]
#   cat(sprintf("\n%s genes mapped in String (%s out of %s) \n", percent(length(ids) / length(gns)), length(ids), length(gns)))
#   cat("\nGenes not mapped in String: ", Genes_not_mapped)
#   string_db_interactions <-  string_db$get_interactions(ids)
#   
#   interactions <- string_db_interactions
#   interactions$from <- string_proteins$preferred_name[match(interactions$from, string_proteins$protein_external_id)]
#   interactions$to <- string_proteins$preferred_name[match(interactions$to, string_proteins$protein_external_id)]
#   
#   interactions <- interactions[, c("from", "to", "coexpression_transferred", "experiments", "database")]
#   interactions <- interactions[which(rowSums(interactions[, -c(1,2)]) > 0), ]
#   interactions_rev <- interactions
#   colnames(interactions_rev) <- colnames(interactions)[c(2,1,3:ncol(interactions))]
#   interactions <- rbind(interactions, interactions_rev)
#   rm(interactions_rev)
#   
#   GeneA <- as.factor(interactions$from)
#   GeneB <- as.factor(interactions$to)
#   AxB <- sparseMatrix(as.numeric(GeneA), as.numeric(GeneB),x=1)
#   rownames(AxB) <- levels(GeneA)
#   colnames(AxB) <- levels(GeneB)
#   dim(AxB)
#   
#   Nmap <- setdiff(gns, colnames(AxB))
#   aux1 <- matrix(0, ncol = length(Nmap), nrow = nrow(AxB))
#   colnames(aux1) <- Nmap
#   rownames(aux1) <- rownames(AxB)
#   
#   AxB <- cbind(AxB, aux1)
#   
#   aux2 <- matrix(0, ncol = ncol(AxB), nrow = length(Nmap))
#   colnames(aux2) <- colnames(AxB)
#   rownames(aux2) <- Nmap
#   
#   AxB <- rbind(AxB, aux2)
#   
#   AxB <- as.matrix(AxB)
#   
#   diag(AxB) <- 1
#   
#   
#   dim(AxB)
#   dim(MagicCube_Trans)
#   iR <- match(rownames(MagicCube_Trans), rownames(AxB))
#   
#   gnNm_Tr <- getBM$Gene_name[match(colnames(MagicCube_Trans), getBM$Transcript_ID)]
#   iC <- match(gnNm_Tr, colnames(AxB))
#   
#   AxB_m <- AxB[iR, iC]
#   
#   MagicCube_Trans[, ,"String"] <- AxB_m
#   
# }else{ MagicCube_Trans[, ,"Found in String"] <- NA; MagicCube_Trans[, ,"String"] <- NA }


# return(list(MagicCube_Trans = MagicCube_Trans, Fit_trans_list = Fit_trans_list))
