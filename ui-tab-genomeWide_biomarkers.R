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



tabPanel(title = "4. Predict Genome-Wide Biomarkers",
         sidebarPanel(width = 3,
                      h3("Your Analysis"),
                      h5("Primary sites selected:"),
                      h6(textOutput("O_txtout4")),
                      h5("Number of samples:"),
                      h6(textOutput("O_nCL_sel4")),
                      h5("Number of essential genes:"),
                      h6(textOutput("O_txtout_nEsGn4")),
                      br(),
                      br(),
                      h3("Parameters"),
                      sliderInput(inputId = "I_co_demeter_4", "Essentiality cutoff (DEMETER score)", min = -5, max = 0, value = -2, step = 0.1),
                      sliderInput(inputId = "I_qntl_exp_tr_4", "Quantile filter of transcripts", min = 0, max = 1, value = 0.5),
                      br(),
                      br(),
                      h3(textOutput("O_txtout_nEsGn_txt_4")),
                      actionButton("I_aB_predictBmrk_genomeWide", "Go! (takes a few minutes)")),
         
         mainPanel(# co_FDR = 0.4, co_log2FC = 2
           tabsetPanel(
             tabPanel("Ranking of biomarkers",
                      h3("Number of Pairs:"),
                      h4(textOutput("O_n_bmrks_all")),
                      h3("a) Ranking of genome-wide biomarkers:"),
                      br(),
                      h4("Filters:"),
                      fluidRow(column(width = 2, numericInput("I_co_PV_tr_4", "P-value < ...", min = 0, value =  0.001, max = 0.05, step = 0.00001)),
                               column(width = 2, numericInput("I_co_logFC_tr_4", "|Log2FC| > ...", min = .5, max = 10, value =  1.5, step = .5)),
                               column(width = 2, numericInput("I_co_lfdr_tr_4", "Local FDR < ...", min = 0, max = 0.95, value =  0.6,  step = 0.1))),
                      br(),
                      downloadButton("O_selected_bmrk_download_all","Download table"),
                      br(),
                      br(),
                      div(withSpinner(DT::dataTableOutput(outputId = "O_rT_Pairs_tr_all")), style = "font-size: 95%; width: 40%"),
                      h3("b) Plot biomarker (select a row):"),
                      downloadButton("plotA2_all_p","Download plot"),
                      fluidRow(
                        column(width = 6, withSpinner(plotOutput("O_gg_bmrk_1_all"))),
                        column(width = 6, withSpinner(plotOutput("O_gg_bmrk_2_all")))),
                      fluidRow(
                        column(width = 6, withSpinner(plotOutput("O_gg_bmrk_5_all"))),
                        column(width = 6, withSpinner(plotOutput("O_gg_bmrk_3_all"))))
                      
             )
             
             
           )
         )
)

