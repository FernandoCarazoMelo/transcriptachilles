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



##### select genes that are expressed in at least one cell line, putative essential genes
# gen expression superior to threshold
# gen expression existent for essential genes in the line for which it is essential
# list of genes 

F_EG_Selection<- function(Essentiality,Exp_gn,iCL,co_enr=1.5,co_pess_iCL=0.2,co_exp=1,
                          qntl_exp=0.75,co_demeter=-2, impute = FALSE){                    
  
  
  prmtrs<-rep(0,5)
  names(prmtrs)<-c("co_demeter","qntl_exp","co_exp","co_pess_iCL","co_enr")
  prmtrs["co_demeter"]<-co_demeter   # -2 essentiality threshold,under-2 a gene is labeled essential
  prmtrs["qntl_exp"]<- qntl_exp  #0.75 min % of CL for which a gene has expression
  prmtrs["co_exp"]<-co_exp  #1min epxresssion level to consider a gene as expressed
  prmtrs["co_pess_iCL"]<- co_pess_iCL   #0.2 min % of CL for which a gene has to be essential
  prmtrs["co_enr"]<- co_enr   #1.5 min enrichment ratio                      
  
  if(impute){  
    require(impute)
    
    output <- impute.knn(as.matrix(Essentiality) ,k = 10, rowmax = 0.75, colmax = 0.8, maxp = 1500, rng.seed=362436069)
    Essentiality <-output$data
    rm(output)
  }
  
  Essentiality_d<-(Essentiality<prmtrs["co_demeter"])*1
  Essentiality_iCL<-Essentiality[,iCL]
  Essentiality_d_iCL<-Essentiality_d[,iCL]
  
  
  Exp_gn_log2<-log2(Exp_gn+1)
  
  iM <- match(rownames(Essentiality_d_iCL), rownames(Exp_gn_log2))
  Essentiality_d_iCL <- Essentiality_d_iCL[!is.na(iM), ]
  Essentiality_d <- Essentiality_d[!is.na(iM), ]
  
  iM <- iM[!is.na(iM)]
  Exp_gn_log2_dem_iCL<-Exp_gn_log2[iM,iCL]
  
  Exp_gn_log2_dem_iCL_d<-(Exp_gn_log2_dem_iCL > log2(prmtrs["co_exp"] + 1))*1
  dim(Essentiality_d_iCL)==dim(Exp_gn_log2_dem_iCL)
  identical(rownames(Essentiality_d_iCL), rownames(Exp_gn_log2_dem_iCL))
  
  check_1<-Exp_gn_log2_dem_iCL_d*Essentiality_d_iCL
  check<-matrix(0,ncol=3,nrow=nrow(Essentiality_d_iCL))
  rownames(check)<-rownames((Essentiality_d_iCL))
  check[,1]<-rowSums(check_1)
  check[,2]<-rowSums(Essentiality_d_iCL)
  check[,3]<-(check[,1]*prmtrs["qntl_exp"]<check[,2])*1
  
  i_Gn<-which(check[,3]>0)
  
  
  ###### Selection of essential genes
  ##### Enrichment and % expression
  
  Enrichment_table<-matrix(0,ncol=6,nrow = length(i_Gn))
  colnames(Enrichment_table)<-c("Expected_nCL","Observed_nCL","Enrichment","%_iCL","%Rest_CL","EG")
  rownames(Enrichment_table)<-rownames(check)[i_Gn]
  n_CL_dem<-ncol(Essentiality_d)
  n_iCL<-length(iCL)
  Essentiality_frequency<-rowSums(Essentiality_d[i_Gn,])/n_CL_dem
  Essentiality_frequency_rest<-rowSums(Essentiality_d[i_Gn,-iCL])/(n_CL_dem-n_iCL)
  Enrichment_table[,1]<-t(Essentiality_frequency*n_iCL)
  Enrichment_table[,2]<-check[i_Gn,2]
  Enrichment_table[,3]<-Enrichment_table[,2]/Enrichment_table[,1]
  Enrichment_table[,4]<-Enrichment_table[,2]/n_iCL
  Enrichment_table[,5]<-Essentiality_frequency_rest
  Enrichment_table[,6]<-(Enrichment_table[,3]>=prmtrs["co_enr"] & Enrichment_table[,4]>prmtrs["co_pess_iCL"])*1
  
  
  EG<-which(Enrichment_table[,6]>0)
  EG<-rownames(Enrichment_table)[EG]
  EG_table_iCL<-Enrichment_table[EG,]
  length(EG)
  
  
  colnames(EG_table_iCL) <- c("Expected_nCL", "Observed_nCL", "Enrichment", "Percent_Ess_Sel" , "Percent_Ess_Rest", "EG")   
  EG_table_iCL <- EG_table_iCL[, c("Percent_Ess_Sel", "Percent_Ess_Rest", "Expected_nCL", "Observed_nCL", "Enrichment")]
  
  EG_table_iCL <- data.frame(Gene_Ess = rownames(EG_table_iCL), EG_table_iCL)
  # rownames(EG_table_iCL) <- NULL
  
  return(EG_table_iCL)
  
}




