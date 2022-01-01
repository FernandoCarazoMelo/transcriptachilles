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



tabPanel("1. Select CLs",
         sidebarPanel(width = 3,
           h4("Your Analysis"),
           h5("Primary sites selected:"),
           h6(textOutput("O_txtout")),
           h5("Number of samples:"),
           h6(textOutput("O_nCL_sel")),
           br(),
           checkboxGroupInput(inputId = "I_primarySite", 
                              label = "Primary Site", 
                              choices = sort(unique(CCLE_sample_info_match$Primary.Site)), 
                              selected = c("kidney")),
           uiOutput("O_renderUI_subt")
         ),
         
         mainPanel(
           tabsetPanel(
             tabPanel(title = "Selected CLs",
                      h2("All Cell lines"),
                      withSpinner(plotOutput("O_gg_CL", height = 300, width = 1000)),
                      h2("Selected Cell lines"),
                      withSpinner(plotOutput("O_gg_CL_sbt", height = 300, width = 1000))),
             tabPanel(title = "List of Selected CLs",
                      br(),
                      uiOutput("O_renderUI_Exclude_CLs"),
                      downloadButton("O_selected_CLs_download","Download the data"),
                      br(),
                      br(),
                      DT::dataTableOutput("O_selected_CLs")))))