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



########################################################################################################
#  GLOBAL
########################################################################################################

source("helpers.R")

createLink <- function(tr_id) {
  sprintf('<a href="http://mar2016.archive.ensembl.org/id/%s" target="_blank" class="btn btn-primary btn-xs">Ensembl</a>',tr_id)
}

createLink2 <- function(tr_id) {
  sprintf('<a href="http://mar2016.archive.ensembl.org/id/%s" TARGET="_blank">%s</a>',tr_id, tr_id)
}


# Load data ---------------------------------------------------------------

load("dataShiny/DataForShiny2_opt.RData")

# load("./shinyApp/transcriptachilles/dataShiny/DataForShiny2_opt_small.RData")



# load("./shinyApp/transcriptachilles/dataShiny/DataForShiny2_opt.RData")
# ## save(CCLE_demeter_match, CCLE_sample_info_match, CCLE_tpm_gn_match, CCLE_tpm_match, getBM, file = "./code/transcriptachilles/dataShiny/DataForShiny2.RData")
# library(shiny)
# library(shinythemes)
# library(ggplot2)
# library(scales)
# library(dplyr)
# library(Matrix)
# library(matrixStats)
# library(tidyr)
# library(readr)
# library(psych)
# library(pheatmap)
# library(RColorBrewer)
# library(progress)
# library(grid)
# library(impute)
# library(shinycssloaders)
# library(limma)
# library(STRINGdb)
# library(qvalue)
# library(ROCR)
# library(plotROC)
# 
# 
# # ------------------------------------------------------------
# source("./shinyApp/transcriptachilles/fun-auxFunctions_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-boxplotMutation_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-F_EG_Selection_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-F_filter_tr_gn_and_getGN_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-F_get_exp_genes_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-F_predictBiomarkers_trans_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-getPairs_tr2_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-ggmatplot_ModifySource.R")
# source("./shinyApp/transcriptachilles/fun-multiplot_ModifySource.R")
# iCL <- which(CCLE_sample_info_match$Primary.Site == "kidney")



########################################################################################################
#  USER INTERFACE
########################################################################################################


# Define UI for application that draws a histogram
ui <- tagList(
  tags$head(includeScript("google-analytics.js")),
  
  
  # sandstone, spacelab
  navbarPage(theme = shinytheme("spacelab"),  title = "transcriptAchilles",
             
             # Overview
             source("ui-tab-overview.R",local=TRUE)$value,
             
             # 0: SELECTION OF CELL LINES
             source("ui-tab-selection_Of_CLs.R",local=TRUE)$value,
             
             # 1: FIND ESSENTIAL GENES
             source("ui-tab-find_Ess_genes.R",local=TRUE)$value,
             
             # 2: PREDICT BIOMARKERS OF A TARGET GENE/DRUG
             source("ui-tab-predict_biomarkers_of_a_target_gene.R",local=TRUE)$value,
             
             # 3: GENOME-WIDE PREDICTIVE ANALYSIS
             source("ui-tab-genomeWide_biomarkers.R",local=TRUE)$value,
             
             # 4: VISUALIZE CASE-BY-CASE
             source("ui-tab-visualize-case-by-case.R",local=TRUE)$value
             # ,
             # # 5: help
             # source("ui-tab-help.R",local=TRUE)$value
  ),
  
  # Footer
  source("ui-footer.R",local=TRUE)$value
  
)



########################################################################################################
#  SERVER
########################################################################################################


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  updateSelectizeInput(session = session, inputId = "I_GeneE", selected = "UBR4", choices = sort(rownames(CCLE_demeter_match)), server = TRUE)
  updateSelectizeInput(session = session, inputId = "I_Marker", selected = "BCR-001 | ENST00000305877", choices = paste(getBM$Transcrip_name, getBM$Transcript_ID, sep = " | "), server = TRUE)
  
  
  
  # 0: SELECTION OF CELL LINES---------------------------------------------------------------------------
  
  iCL <- reactive({ which(CCLE_sample_info_match$Primary.Site %in% input$I_primarySite) })
  
  output$O_gg_CL <- renderPlot({
    tb <- within(CCLE_sample_info_match,
                 Primary.Site <- factor(Primary.Site,
                                        levels=names(sort(table(Primary.Site),
                                                          decreasing=TRUE))))
    ggplot(tb, aes(Primary.Site)) + geom_bar() + ggtitle("All Cell Lines (n = 412)") + scale_fill_discrete(guide=FALSE) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  
  output$O_gg_CL_sbt <- renderPlot({
    tb2 <- within(CCLE_sample_info_match[iCL2(), ],
                  Primary.Site <- factor(Primary.Site,
                                         levels=names(sort(table(Primary.Site),
                                                           decreasing=TRUE))),
                  Subtype <- factor(Subtype,
                                    levels=names(sort(table(Subtype), decreasing=TRUE))))
    ggplot(tb2, aes(Primary.Site)) + geom_bar(aes(fill = Subtype)) + ggtitle(sprintf("Selected Cell Lines (n = %s)", length(iCL2()))) + theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$O_txtout <- renderText({paste(input$I_primarySite)})
  output$O_txtout2 <- renderText({paste(input$I_primarySite)})
  output$O_txtout3 <- renderText({paste(input$I_primarySite)})
  output$O_txtout4 <- renderText({paste(input$I_primarySite)})
  
  output$O_nCL_sel <- renderText({length(iCL2())})
  output$O_nCL_sel2 <- renderText({length(iCL2())})
  output$O_nCL_sel3 <- renderText({length(iCL2())})
  output$O_nCL_sel4 <- renderText({length(iCL2())})
  
  output$O_renderUI_subt <- renderUI({
    checkboxGroupInput("I_subtype",  label = "Subtype", 
                       choices = sort(unique(CCLE_sample_info_match$Subtype[iCL()])),
                       selected = unique(CCLE_sample_info_match$Subtype[iCL()]))
    
  })
  
  output$O_renderUI_Exclude_CLs <- renderUI({ 
    selectInput(inputId = "I_Exclude_CLs",
                label = "Exclude CLs from the selection", 
                choices = unique(CCLE_sample_info_match$Cell_Line[iCL()]), 
                multiple = T, 
                selectize = T)
  })
  
  
  iCL2 <- reactive({ which((CCLE_sample_info_match$Primary.Site %in% input$I_primarySite) & (CCLE_sample_info_match$Subtype %in% input$I_subtype) & !(CCLE_sample_info_match$Cell_Line %in% input$I_Exclude_CLs)) })
  
  # Data Table of Selected Cell lines
  
  selected_CLs <- reactive({ CCLE_sample_info_match[iCL2(), c("Cell_Line", "CCLE.name", "Gender", "Primary.Site", "Subtype")] })
  output$O_selected_CLs <- DT::renderDataTable(filter = "top", caption = "Table 1: Summary of selected cell lines.",
                                               options = list(pageLength = 5, autoWidth = TRUE),
                                               { selected_CLs() })
  output$O_selected_CLs_download <- downloadHandler(filename = function(){"TranscriptAchilles_SelectedCLs.csv"}, 
                                                    content = function(fname){write.csv2(selected_CLs(), fname, quote = F, row.names = F)})
  
  
  # 1. FIND ESSENTIAL GENES--------------------------------------------------------------------------------------------
  # input <- NULL; input$I_co_demeter <- -2; input$I_co_pess_iCL <- 0.25; input$I_co_enr <- 2; input$I_co_exp <- 1; input$I_qntl_exp <- 0.75
  
  # Build essential gene table
  esGene <- reactive({
    
    # if(!exists("iCL2()")){showNotification("ERROR. \n Please, select the Cell Lines.", type = "error")
    # }else{}
    aux <- as.data.frame(F_EG_Selection(Essentiality = CCLE_demeter_match,
                                        Exp_gn = CCLE_tpm_gn_match,
                                        iCL = iCL2(),
                                        co_demeter= input$I_co_demeter,
                                        co_pess_iCL = input$I_co_pess_iCL,
                                        co_enr= input$I_co_enr,
                                        co_exp = input$I_co_exp,
                                        qntl_exp= input$I_qntl_exp,
                                        impute = FALSE))
    aux <- aux[order(aux$Percent_Ess_Rest, decreasing = F), ]
    aux <- aux[order(aux$Percent_Ess_Sel, decreasing = T), ]
    aux
    
  })
  
  # esGene <- aux
  
  output$O_renderTable_EsGenes <- DT::renderDataTable(selection = list(mode = "single", selected = 1) , caption = "Table 2: Ranking of essential genes for selected cell lines.",
                                                      options = list(pageLength = 5, autoWidth = TRUE), 
                                                      {
                                                        aux <- esGene()
                                                        aux[,-1] <- round(aux[, -1],2)
                                                        aux[, 2] <- percent(aux[, 2])
                                                        aux[, 3] <- percent(aux[, 3])
                                                        aux$Enrichment[is.infinite(aux$Enrichment)] <- 999
                                                        rownames(aux) <- NULL
                                                        aux
                                                      })
  
  
  output$O_esGenes_download <- downloadHandler(filename = function(){"TranscriptAchilles_EssentialGenes.csv"},
                                               content = function(fname){
                                                 write.csv2(esGene(), fname, quote = F)})
  
  output$O_txtout_nEsGn <- renderText({nrow(esGene())})
  output$O_txtout_nEsGn2 <- renderText({nrow(esGene())})
  output$O_txtout_nEsGn3 <- renderText({nrow(esGene())})
  output$O_txtout_nEsGn4 <- renderText({nrow(esGene())})
  output$O_txtout_nEsGn_txt_4 <- renderText({sprintf("Predict Best Biomarkers of %s Essential Genes", nrow(esGene()))})
  
  iRow <- reactive({input$O_renderTable_EsGenes_rows_selected})
  
  output$o_Plot_Ess <- renderPlot({
    Gene_E <- as.character(esGene()$Gene_Ess[iRow()])
    boxplotEssentiality(Essentiality = CCLE_demeter_match, iCL = iCL2(), Gene_E = Gene_E, co_Ess = input$I_co_demeter)
  })
  
  output$o_Plot_Ess_p <- downloadHandler(filename = function(){paste("TranscriptAchilles_BoxplotEssentiality.pdf")}, 
                                         content = function(fname){
                                           pdf(fname, width = 12, height = 4.3)
                                           Gene_E <- as.character(esGene()$Gene_Ess[iRow()])
                                           boxplotEssentiality(Essentiality = CCLE_demeter_match, iCL = iCL2(), Gene_E = Gene_E, co_Ess = input$I_co_demeter)
                                           dev.off()})
  

  
  
  # 2. PREDICT BIOMARKERS OF A TARGET GENE------------------------------------------------------------------------------------
  # input$I_qntl_exp_tr <- 0.5; input$I_co_demeter2 <- -2; input$I_GeneE_to_predict_bmkr <- "IRAK1"; iGn_sel <- match(input$I_GeneE_to_predict_bmkr, rownames(CCLE_demeter_match)); input$O_rT_Pairs_tr_rows_selected <- 1
  
  
  output$O_GeneE_to_predict_bmkr <- renderUI({
    selectInput(inputId = "I_GeneE_to_predict_bmkr",
                label = "Select One or Multiple Genes:", 
                selected = rownames(esGene())[1],
                choices = rownames(esGene()),
                multiple = T, selectize = T)
  })
  # updateSelectizeInput(session = session, inputId = "I_GeneE_to_predict_bmkr", selected = rownames(esGene())[1], choices = rownames(esGene()), server = TRUE)
  
  
  iGn_sel <- reactive({
    match(input$I_GeneE_to_predict_bmkr, rownames(CCLE_demeter_match))
    
  })
  
  predictBiomarkers_tr <- eventReactive(input$I_aB_predictBmrk, {
    withProgress(message = "Predicting Biomarkers: ", value = 0, {
      incProgress(0, detail = "preparing data")
      Filt <- filter_tr_gn_and_getGN(Exp_trancripts = CCLE_tpm_match, getBM = getBM, 
                                     iCL = iCL2(), qt = input$I_qntl_exp_tr)
      
      iT_exp <- match(Filt$T_filt, rownames(CCLE_tpm_match))
      iG_exp <- match(Filt$G_filt, rownames(Filt$Exp_gn))
      Exp_genes <- Filt$Exp_gn
      predictBiomarkers_tr <- F_predictBiomarkers_trans(Essentiality = CCLE_demeter_match,
                                                        Exp_trancripts = CCLE_tpm_match,
                                                        Exp_genes = Exp_genes,
                                                        iGn_ess = iGn_sel(),
                                                        iT_exp = iT_exp,
                                                        iG_exp = iG_exp,
                                                        iCL = iCL2(),
                                                        getBM = getBM,
                                                        co_demeter = input$I_co_demeter2,
                                                        impute = F, restCL = F, string = F, shiny = T)
      predictBiomarkers_tr
    })
  })
  
  # Pairs_tr <- getPairs_tr_list(predictBiomarkers_tr$Tr_Bmrk_list, getBM = getBM, co_log2FC = .5, co_PV = 0.05)
  Pairs_tr <- reactive({
    getPairs_tr_list(predictBiomarkers_tr()$Tr_Bmrk_list, getBM = getBM, co_log2FC = .5, co_PV = 0.05)
    
  })
  
  nPairs <- reactive({nrow(Pairs_tr())})
  
  
  iSel <- reactive({
    Pairs_print <- Pairs_tr()
    which((Pairs_print$P.Value < input$I_co_PV_tr) &
            (abs(Pairs_print$logFC) > input$I_co_logFC_tr) &
            (Pairs_print$lfdr < input$I_co_lfdr_tr))
  })
  
  output$O_rT_Pairs_tr <- DT::renderDataTable(filter = 'top',escape = FALSE, selection = list(mode = "single", selected = 1), caption = "Table 3: Ranking of biomarkers for selected gene. The column Group_bmkr (Group Biomarker) represents whether the biomarker is a Transcript, a gene or a gene with only one transcript.",
                                              options = list(pageLength = 5, autoWidth = TRUE),
                                              {
                                                Pairs_print <- Pairs_tr()
                                                Pairs_print[,c(7,8,10)] <- round(Pairs_print[,c(7,8,10)], 2)
                                                Pairs_print$P.Value <- sprintf("%.2e",Pairs_print$P.Value)
                                                Pairs_print <- Pairs_print[iSel(), ]
                                                Pairs_print$Transcript_ID <- createLink2(Pairs_print$Transcript_ID)
                                                iOut <- which(Pairs_print$Group_bmkr =="gene")
                                                iOut2 <- which(Pairs_print$Group_bmkr %in% c("transcript", "gene_one_Tr"))
                                                Pairs_print$Transcript_ID[iOut] <- "N/A"
                                                Pairs_print$Transcrip_bmkr[iOut] <- "N/A"
                                                Pairs_print$Biotype[iOut] <- "N/A"
                                                Pairs_print$Gene_bmkr[iOut2] <- "N/A"
                                                Pairs_print <- Pairs_print %>% dplyr::select(-AveExpr)
                                              })
  
  output$O_selected_bmrk_download <- downloadHandler(filename = function(){"TranscriptAchilles_PredictedBmrkrs.csv"}, 
                                                     content = function(fname){write.csv2(Pairs_tr()[iSel(), ], fname, quote = F, row.names = F)})
  
  output$O_n_bmrks <- renderText({paste(length(iSel()))})
  
  output$O_gg_qvalue_gn <- renderPlot({
    hist(predictBiomarkers_tr()$qvalue_Gn)
  })
  
  output$O_gg_qvalue_tr <- renderPlot({
    hist(predictBiomarkers_tr()$qvalue_Tr)
  })
  
  
  output$plot_qvalue_p <- downloadHandler(filename = function(){paste("TranscriptAchilles_localFDR.pdf")}, 
                                     content = function(fname){
                                       pdf(fname, width = 10, height = 4.5); 
                                       multiplot(plotlist = list(hist(predictBiomarkers_tr()$qvalue_Gn) + ggtitle("Genes", subtitle = "P-value density histogram"), hist(predictBiomarkers_tr()$qvalue_Tr) + ggtitle("Transcripts", subtitle = "P-value density histogram")),cols=2);
                                       dev.off()})
  
  
  
  plotA2 <- reactive({
    Pairs <- Pairs_tr()[iSel(), ]
    rBmrk <- Pairs[input$O_rT_Pairs_tr_rows_selected, ]
    
    tr <- as.character(rBmrk$Transcript_ID)
    gn <- as.character(rBmrk$Gene_Ess)
    dm <- input$I_co_demeter2
    group_bmrk <- as.character(rBmrk$Group_bmkr)
    plot_transcripts_essential2(tr = tr, group_bmrk = group_bmrk,
                                gn = gn, getBM = getBM,
                                CCLE_ataris = CCLE_demeter_match,
                                CCLE_tpm_log2 = log2(1+CCLE_tpm_match),
                                CCLE_sample_info = CCLE_sample_info_match,
                                iCL = iCL2(), text = F, subtype = F, vline = dm)
  })
  
  output$plotA2_p <- downloadHandler(filename = function(){paste("TranscriptAchilles_PlotPairs.pdf")}, 
                                                    content = function(fname){
                                                      pdf(fname, width = 12, height = 10); 
                                                      multiplot(plotlist = list(plotA2()$gg1,plotA2()$gg2, plotA2()$gg5, plotA2()$gg3),cols=2);
                                                      dev.off()})
  
  
  ## Plot isoforms
  
  output$O_gg_bmrk_1 <- renderPlot({
    plotA2()$gg1
  })
  output$O_gg_bmrk_2 <- renderPlot({
    plotA2()$gg2
  })
  output$O_gg_bmrk_3 <- renderPlot({
    plotA2()$gg3
  })
  output$O_gg_bmrk_4 <- renderPlot({
    plotA2()$gg4
  })
  output$O_gg_bmrk_5 <- renderPlot({
    plotA2()$gg5
  })
  
  
  # 3: PREDICT GENOME-WIDE BIOMARKERS-------------------------------------------------------------------------
  # input$I_co_demeter_4 <- -2; input$O_rT_Pairs_tr_all_rows_selected <- 1
  
  predictBiomarkers_tr_all <- eventReactive(input$I_aB_predictBmrk_genomeWide, {
    limit <- 125 # 125
    if (nrow(esGene()) > limit){
      showNotification(sprintf("%s genes selected. The limit is %s. Running for top %s genes.", nrow(esGene()), limit, limit), type = "warning", duration = 15)
      esGene_red <- esGene()[1:limit, ]
    }else{
      esGene_red <- esGene()
    }
    withProgress(message = "Predicting Biomarkers: ", value = 0, {
      incProgress(0, detail = "preparing data")
      Filt <- filter_tr_gn_and_getGN(Exp_trancripts = CCLE_tpm_match, getBM = getBM, iCL = iCL2(), qt = input$I_qntl_exp_tr_4)
      
      iT_exp <- match(Filt$T_filt, rownames(CCLE_tpm_match))
      iG_exp <- match(Filt$G_filt, rownames(Filt$Exp_gn))
      Exp_genes <- Filt$Exp_gn
      
      iGn <- match(rownames(esGene_red), rownames(CCLE_demeter_match))
      
      predictBiomarkers_tr_all <- F_predictBiomarkers_trans(Essentiality = CCLE_demeter_match,
                                                            Exp_trancripts = CCLE_tpm_match,
                                                            Exp_genes = Exp_genes,
                                                            iGn_ess = iGn,
                                                            iT_exp = iT_exp,
                                                            iG_exp = iG_exp,
                                                            iCL = iCL2(),
                                                            getBM = getBM,
                                                            co_demeter = input$I_co_demeter_4,
                                                            impute = F, restCL = F, string = F, shiny = T)
      predictBiomarkers_tr_all
    })
    
    
  })
  
  
  Pairs_tr_all <- reactive({
    getPairs_tr_list(predictBiomarkers_tr_all()$Tr_Bmrk_list, getBM = getBM, co_log2FC = .5, co_PV = 0.05)
    
  })
  
  nPairs_all <- reactive({nrow(Pairs_tr_all())})
  
  
  iSel_all <- reactive({
    Pairs_print <- Pairs_tr_all()
    which((Pairs_print$P.Value < input$I_co_PV_tr_4) &
            (abs(Pairs_print$logFC) > input$I_co_logFC_tr_4) &
            (Pairs_print$lfdr < input$I_co_lfdr_tr_4))
  })
  
  output$O_rT_Pairs_tr_all <- DT::renderDataTable(filter = 'top', escape = FALSE, selection = list(mode = "single", selected = 1), caption = "Table 4: Genome-Wide ranking of pairs (Essential Gene & Biomarkers). The column Group_bmkr (Group Biomarker) represents whether the biomarker is a Transcript, a gene or a gene with only one transcript.",
                                                  options = list(pageLength = 5, autoWidth = TRUE),
                                                  {
                                                    Pairs_print <- Pairs_tr_all()
                                                    Pairs_print[,c(7,8,10)] <- round(Pairs_print[,c(7,8,10)], 2)
                                                    Pairs_print$P.Value <- sprintf("%.2e",Pairs_print$P.Value)
                                                    Pairs_print <- Pairs_print[iSel_all(), ]
                                                    Pairs_print$Transcript_ID <- createLink2(Pairs_print$Transcript_ID)
                                                    iOut <- which(Pairs_print$Group_bmkr =="gene")
                                                    iOut2 <- which(Pairs_print$Group_bmkr %in% c("transcript", "gene_one_Tr"))
                                                    Pairs_print$Transcript_ID[iOut] <- "N/A"
                                                    Pairs_print$Transcrip_bmkr[iOut] <- "N/A"
                                                    Pairs_print$Biotype[iOut] <- "N/A"
                                                    Pairs_print$Gene_bmkr[iOut2] <- "N/A"
                                                    Pairs_print <- Pairs_print %>% dplyr::select(-AveExpr)
                                                  })
  
  output$O_n_bmrks_all <- renderText({paste(length(iSel_all()))})
  
  
  # Dowload  
  output$O_selected_bmrk_download_all <- downloadHandler(filename = function(){"TranscriptAchilles_GenomeWideBmrkrs.csv"}, 
                                                         content = function(fname){write.csv2(Pairs_tr_all()[iSel_all(), ], fname, quote = F, row.names = F)})
  
  
  plotA2_all <- reactive({
    Pairs <- Pairs_tr_all()[iSel_all(), ]
    rBmrk <- Pairs[input$O_rT_Pairs_tr_all_rows_selected, ]
    
    tr <- as.character(rBmrk$Transcript_ID)
    gn <- as.character(rBmrk$Gene_Ess)
    dm <- input$I_co_demeter_4
    group_bmrk <- as.character(rBmrk$Group_bmkr)
    
    plot_transcripts_essential2(tr = tr,
                                gn = gn, getBM = getBM, group_bmrk = group_bmrk,
                                CCLE_ataris = CCLE_demeter_match,
                                CCLE_tpm_log2 = log2(1+CCLE_tpm_match),
                                CCLE_sample_info = CCLE_sample_info_match,
                                iCL = iCL2(), text = F, subtype = F, vline = dm)
  })
  
  output$plotA2_all_p <- downloadHandler(filename = function(){paste("TranscriptAchilles_PlotPairs.pdf")}, 
                                     content = function(fname){
                                       pdf(fname, width = 12, height = 10); 
                                       multiplot(plotlist = list(plotA2_all()$gg1,plotA2_all()$gg2, plotA2_all()$gg5, plotA2_all()$gg3),cols=2);
                                       dev.off()})
  
  
  ## Plot isoforms
  
  output$O_gg_bmrk_1_all <- renderPlot({
    plotA2_all()$gg1
  })
  output$O_gg_bmrk_2_all <- renderPlot({
    plotA2_all()$gg2
  })
  output$O_gg_bmrk_3_all <- renderPlot({
    plotA2_all()$gg3
  })
  output$O_gg_bmrk_5_all <- renderPlot({
    plotA2_all()$gg5
  })
  
  
  # 4: CASE BY CASE-----------------------------------------------------------------------------------------------
  
  plotA <- reactive({
    trN <- unlist(strsplit(x = input$I_Marker, split = "[ ][|]"))[1]
    tr <-  getBM$Transcript_ID[which(getBM$Transcrip_name == trN)]
    gn <- input$I_GeneE

    plot_transcripts_essential2(tr = tr,
                                gn = gn, getBM = getBM,
                                CCLE_ataris = CCLE_demeter_match,
                                CCLE_tpm_log2 = log2(1+CCLE_tpm_match),
                                CCLE_sample_info = CCLE_sample_info_match,
                                iCL = iCL2(), text = F, subtype = F, vline = input$I_sld_demeter2)
  })
  
  
  output$plotA_p <- downloadHandler(filename = function(){paste("TranscriptAchilles_PlotPairs.pdf")}, 
                                    content = function(fname){
                                      pdf(fname, width = 12, height = 10); 
                                      multiplot(plotlist = list(plotA()$gg1,plotA()$gg2, plotA()$gg5, plotA()$gg3),cols=2);
                                      dev.off()})
  
  
  ## Plot isoforms
  
  output$O_gg_1 <- renderPlot({
    plotA()$gg1
  })
  output$O_gg_2 <- renderPlot({
    plotA()$gg2
  })
  output$O_gg_3 <- renderPlot({
    plotA()$gg3
  })
  output$O_gg_4 <- renderPlot({
    plotA()$gg4
  })
  output$O_gg_5 <- renderPlot({
    plotA()$gg5
  })
}



# Run the application 
shinyApp(ui = ui, server = server)