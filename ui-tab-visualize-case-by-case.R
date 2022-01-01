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



tabPanel(title = "5. Visualize Case-by-Case",
         
         fluidRow(
           column(width = 4, selectizeInput(inputId = "I_GeneE", label = "Essential Gene (type and wait)", choices = NULL, options = list(maxOptions = 20))),
           column(width = 4, selectizeInput(inputId = "I_Marker", label = "Transcript Biomarker (type and wait)", choices = NULL, options = list(maxOptions = 20))),
           column(width = 4, sliderInput(inputId = "I_sld_demeter2", label = "Essentiality cutoff (Blue lines. If 0: not plotted)", min = -5, max = 0, value = 0, step = 0.1))),
         br(),
         downloadButton("plotA_p","Download plot"),
         br(),
         fluidRow(
           column(width = 6, withSpinner(plotOutput("O_gg_1"))),
           column(width = 6, withSpinner(plotOutput("O_gg_2")))),
         fluidRow(
           column(width = 6, withSpinner(plotOutput("O_gg_5"))),
           column(width = 6, withSpinner(plotOutput("O_gg_3")))),
         fluidRow(column(width = 6, withSpinner(plotOutput("O_gg_4")), offset = 6))
)
