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



ggmatplot <- function(dat, xnames=TRUE, highlight = NULL){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(tibble)
  
  if(!is.matrix(dat)){
    dat <- as.data.frame(dat)
    df <- dat %>% 
      rownames_to_column() %>% 
      gather(reading, value, -rowname) %>% 
      group_by(rowname) %>% 
      dplyr::mutate(x=1:n()) 
    colnames(df)[1] <- "reading"
    colnames(df)[2] <- "Legend"
    
    if(xnames){
      a <- ggplot(data = df, aes(x=1:nrow(df), y=value)) +
        geom_line(aes(color = Legend), size = 1.5) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(limits = df$reading)
    }else{
      a <- ggplot(data = df, aes(x=1:nrow(df), y=value)) +
        geom_line(aes(color = Legend), size = 1.5)
    }
    
  }else{
    dat <- as.data.frame(t(dat))
    df <- dat %>% 
      rownames_to_column() %>% 
      gather(reading, value, -rowname) %>% 
      group_by(rowname) %>% 
      dplyr::mutate(x=1:n()) 
    colnames(df)[1] <- "Legend"
    
    
    if(xnames){
      df$reading <- factor(df$reading,levels = unique(df$reading))
      a <- ggplot(data = df, aes(x=reading, y=value, group=Legend)) +
        geom_line(aes(color=Legend)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{
      a <- ggplot(df, aes(x=x, y=value, group=Legend)) +
        geom_line(aes(color=Legend))
    }
    
    if(length(highlight > 0)){a <- a + geom_line(data = df[which(df$Legend %in% highlight), ],aes(color=Legend), size = 1.5)}
  }
  return(a)
}
