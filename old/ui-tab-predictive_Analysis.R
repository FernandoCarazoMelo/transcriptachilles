tabPanel(title = "2. Predictive analysis",
         tabsetPanel(
           tabPanel("Selection of essential (ES) genes",
                    sidebarPanel(
                      h4("Summary"),
                      h5("Primary sites:"),
                      h6(textOutput("O_txtout2")),
                      h5("Number of samples:"),
                      h6(textOutput("O_nCL_sel2")),
                      h5("Number of essential genes:"),
                      h6(textOutput("O_txtout_nEsGn")),
                      br(),
                      h4("Essentiality filters"),
                      sliderInput(inputId = "I_co_demeter", "Essentiality cutoff (DEMETER score)", min = -5, max = 0, value = -2, step = 0.1),
                      sliderInput(inputId = "I_co_pess_iCL", "Minimal percetage of essentiality (Selected CLs)", min = 0, max = 1, value = 0.2),
                      sliderInput(inputId = "I_co_enr", "Enrichment of essentiality (Expected vs Observed essential CLs)", min = 0, max = 5, value = 1,step = 0.1),
                      br(),
                      h4("Expression filters"),
                      sliderInput(inputId = "I_co_exp", "Minimal required expression (TPM)...", min = 0, max = 20, value = 1),
                      sliderInput(inputId = "I_qntl_exp", "... at least in XX percentage of essential CLs", min = 0, max = 1, value = 0.75)
                      ),
                    
                    mainPanel(
                      h3("Number of essential genes:"),
                      h4(textOutput("O_txtout_nEsGn2")),
                      br(),
                      h3("Table of essential genes"),
                      downloadButton("O_esGenes_download","Download the data"),
                      br(),
                      DT::dataTableOutput(outputId = "O_renderTable_EsGenes"))),
           
           tabPanel(title = "Prediction of transcript biomarkers of selected ES genes"),
           
           tabPanel("Visualize data",
                    fluidRow(
                      column(width = 4, selectizeInput(inputId = "I_GeneE", label = "Select essential gene (type and wait)", choices = NULL, options = list(maxOptions = 10))),
                      column(width = 4, selectizeInput(inputId = "I_Marker", label = "Select transcript biomarker (type and wait)", choices = NULL, options = list(maxOptions = 10))),
                      column(width = 4, sliderInput(inputId = "I_sld_demeter2", label = "Essentiality cutoff (Blue lines. If 0: not plotted)", min = -5, max = 0, value = -2, step = 0.1))),
                    br(),
                    fluidRow(
                      column(width = 6, withSpinner(plotOutput("O_gg_1"))),
                      column(width = 6, withSpinner(plotOutput("O_gg_2")))),
                    fluidRow(column(width = 6, withSpinner(plotOutput("O_gg_3")), offset = 6)),
                    fluidRow(column(width = 6, withSpinner(plotOutput("O_gg_4")), offset = 6))
           )
         )
)
