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



filter_tr_gn_and_getGN <- function(Exp_trancripts, getBM, iCL, qt = 0.5){
  #parameters
  # qntl_exp min % of CL for which a gene has expression
  # co_exp min epxresssion level to consider a gene as expressed
  # co_pess_iCL min % of CL for which a gene has to be essential
  
  # Get gene expression---------------------------------------------
  Exp_genes <- get_exp_genes(Exp_trancripts, getBM)
  
  
  # Selection of transcripts-----------------------------------------
  Exp_trancripts_log2 <- log2(1 + Exp_trancripts)
## Genes with one transcript
  getBM_ntr <- getBM %>% group_by(Gene_ID) %>% dplyr::mutate(tr = n())
  is_mult_tr <- getBM_ntr$tr > 1
  
  ## Non-expressed transcripts
  is_Ex_tr <- rowSums(Exp_trancripts_log2[, iCL]) > 0
  iT_aux <- which(is_mult_tr & is_Ex_tr)
  
  ## Expression quantile
  quan <- quantile(rowSums(Exp_trancripts_log2[iT_aux, iCL]),qt)
  is_quant_tr <- rowSums(Exp_trancripts_log2[, iCL]) > quan
  iT_sel <- which(is_mult_tr & is_Ex_tr & is_quant_tr)
  T_filt <- rownames(Exp_trancripts_log2)[iT_sel]

  # quantile <- quantile(as.numeric(as.matrix(Exp_trancripts_log2[iT_aux, iCL])),qt)
  # iEss <- which(Dmatrix[,2] == 1)
  # co_exp = 1
  # iaux <- which(rowMedians(as.matrix(Exp_trancripts_log2[,iEss])) > co_exp |
  #                 rowMedians(as.matrix(Exp_trancripts_log2[,-iEss])) > co_exp)
  
  
  # Selection of genes-----------------------------------------
  ## Non-expressed transcripts
  
  Exp_genes_log2 <- log2(1 + Exp_genes)
  is_Ex_gn <- rowSums(Exp_genes_log2[, iCL]) > 0

  ## Expression quantile
  quan <- quantile(rowSums(Exp_genes_log2[which(is_Ex_gn), iCL]),qt)
  is_quant_gn <- rowSums(Exp_genes_log2[, iCL]) > quan
  iG_sel <- which(is_Ex_gn & is_quant_gn)
  
  G_filt <- rownames(Exp_genes_log2)[iG_sel]

  
  return(list(Exp_gn = Exp_genes, T_filt = T_filt, G_filt = G_filt))
}