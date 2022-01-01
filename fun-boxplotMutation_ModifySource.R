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



boxplotMutation <- function(Essentiality, CLxMut, iCL, Gene_E, Mutation, tag = ""){
  
  mt <- CLxMut[, which(colnames(CLxMut) == Mutation)]
  
  par(mfrow = c(1, 2))
  
  boxplot(Essentiality[Gene_E, -iCL] ~ mt[-iCL],
               main = sprintf("%s  -  %s   (Rest: n = %s) \n nMut = %s", Mutation, Gene_E, ncol(Essentiality)-length(iCL),sum(mt[-iCL])),
               xlab = paste("Mutation  ", Mutation), 
               ylab = paste("Essentiality  ", Gene_E))
  
  abline(h = -2, lty = 2, col = 2)
  
  boxplot(Essentiality[Gene_E, iCL] ~ mt[iCL],
               main = sprintf("%s  -  %s   (%s: n = %s) \n nMut = %s", Mutation, Gene_E, tag, length(iCL),sum(mt[iCL])),
               xlab = paste("Mutation  ", Mutation), 
               ylab = paste("Essentiality  ", Gene_E))
  abline(h = -2, lty = 2, col = 2)
  
  
  par(mfrow = c(1, 1))
  

}