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



tabPanel(title = "2. Find Essential Genes",
         sidebarPanel(
           h3("Your Analysis"),
           h5("Primary sites selected:"),
           h6(textOutput("O_txtout2")),
           h5("Number of samples:"),
           h6(textOutput("O_nCL_sel2")),
           h5("Number of essential genes:"),
           h6(textOutput("O_txtout_nEsGn")),
           br(),
           br(),
           h3("Essentiality filters"),
           sliderInput(inputId = "I_co_demeter", "Essentiality cutoff (DEMETER score)", min = -5, max = 0, value = -2, step = 0.1),
           sliderInput(inputId = "I_co_pess_iCL", "Minimal percentage of essentiality (Selected CLs)", min = 0, max = 1, value = 0.25),
           sliderInput(inputId = "I_co_enr", "Enrichment of essentiality in CLs (#Observed / #Expected)", min = 0, max = 10, value = 2,step = 0.1),
           br(),
           br(),
           h3("Expression filters"),
           sliderInput(inputId = "I_co_exp", "Minimal required expression of essential genes (TPM)...", min = 0, max = 20, value = 1),
           sliderInput(inputId = "I_qntl_exp", "Percentage of CLs that pass the expression filter", min = 0, max = 1, value = 0.75)
         ),
         
         mainPanel(
           h3("Number of Essential Genes:"),
           h4(textOutput("O_txtout_nEsGn2")),
           h3("a) Ranking of Essential Genes"),
           downloadButton("O_esGenes_download","Download table"),
           br(),
           br(),
           withSpinner(DT::dataTableOutput(outputId = "O_renderTable_EsGenes")),
           br(),
           h3("b) Plot a Gene (select a row)"),
           downloadButton("o_Plot_Ess_p","Download plot"),
           
           withSpinner(plotOutput("o_Plot_Ess"))
         )
)
