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



getPairs_tr_list <- function(Tr_Bmrk_list, getBM, co_PV = 0.001,  co_log2FC = 2){
  
  # co_FDR = 0.4, Dir = FALSE, String = FALSE, trMax = 20
  require(dplyr)
  
  Gene_Ess <- names(Tr_Bmrk_list)
  
  aux2 <- do.call(rbind, Tr_Bmrk_list)
  
  aux2 <- cbind(Gene_Ess = unlist(strsplit( split = "[.]", x = rownames(aux2)))[c(T,F)], aux2)
  
  aux2 <- aux2[which((aux2$P.Value < co_PV) & (abs(aux2$logFC) > co_log2FC)), ]
  
  Pairs_tr <- aux2[order(aux2$P.Value), ]
  
  rownames(Pairs_tr) <- NULL
  
  # ilessTr <- which(getBM$tr < 20)
  
  return(Pairs_tr)
}

