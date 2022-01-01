## ==================================================================================== ##
# transcriptAchilles Shiny App for predicting transcript biomarkers.
# Copyright (C) 2018
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



boxplotEssentiality <- function(Essentiality, iCL, Gene_E, tag = "", co_Ess = -2){
  
  cl <- rep(sprintf("Rest CLs (n = %s)", ncol(Essentiality) - length(iCL)), ncol(Essentiality))
  cl[iCL] <- sprintf("Selected CLs (n = %s)", length(iCL))
  cl <- as.factor(cl)
  cl <- relevel(cl, ref = sprintf("Selected CLs (n = %s)", length(iCL)))
  boxplot(Essentiality[Gene_E, ] ~ cl,
          main = sprintf("Essentiality score of %s",Gene_E),
          xlab = paste("Cell Lines"), 
          ylab = paste("Essentiality score of ", Gene_E, " (DEMETER)"))
  
  abline(h = co_Ess, lty = 2, col = 2)
  par(mfrow = c(1, 1))
}