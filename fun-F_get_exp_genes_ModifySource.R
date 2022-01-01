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



get_exp_genes <- function(Exp_trancripts, getBM){
  Exp_genes <- cbind(getBM[,c("Gene_name","Transcript_ID")], Exp_trancripts)
  Exp_genes <- Exp_genes %>% group_by(Gene_name) %>% summarise_each(funs(sum),-Transcript_ID)
  rownames(Exp_genes) <- Exp_genes[,1]$Gene_name
  Exp_genes <- dplyr::select(Exp_genes, -Gene_name)
  Exp_genes <- as.data.frame(Exp_genes)
  return(Exp_genes)
}