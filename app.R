library(shiny)
if (!require("fgsea", quietly = TRUE))
  BiocManager::install("fgsea")
library(fgsea)
if (!require("data.table", quietly = TRUE))
  install.packages("data.table")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!require("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!require("ggsignif", quietly = TRUE))
  install.packages("ggsignif")
if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")
if (!require("rstatix", quietly = TRUE))
  install.packages("rstatix")
if (!require("ggpattern", quietly = TRUE))
  remotes::install_github("coolbutuseless/ggpattern")



load(".RData")
pathways <- gmtPathways("m2.all.v2023.1.Mm.symbols.gmt.txt")
pathways2 <- gmtPathways("m5.all.v2023.1.Mm.symbols.gmt.txt")
pathways <- c(pathways, pathways2)
set.seed(123)
DefaultAssay(whole.integrated) <- "RNA"
whole.integrated <- NormalizeData(object = whole.integrated, normalization.method = "LogNormalize",
                      scale.factor = 10000)
whole.integrated <- ScaleData(whole.integrated, verbose = FALSE)
# counts <- GetAssay(object = whole.integrated, assay=DefaultAssay(whole.integrated))
# count_matrix <- as.matrix(GetAssayData(object = whole.integrated, assay="RNA"))

#USING THIS TO GET RNA
#GetAssayData(object = count_matrix, assay="RNA")

# umap_pos <- as.data.frame(Embeddings(whole.integrated, reduction = "umap"))
# tsne_pos <- as.data.frame(Embeddings(whole.integrated, reduction = "tsne"))
# data_merg <- whole.integrated@meta.data
# data_merg$cells <- rownames(data_merg)
# umap_pos$cells <- rownames(umap_pos)
# addGesecaScores <- function(pathways,
#                             matrix,
#                             prefix="",
#                             scale=FALSE) {
#   x <- matrix
#   E <- x@scale.data
# 
# 
# 
#   for (i in seq_along(pathways)) {
#     pathway <- pathways[[i]]
#     pathway <- intersect(unique(pathway), rownames(E))
# 
#     score <- colSums(E[pathway, , drop=FALSE])/sqrt(length(pathway))
#     score <- scale(score, center=TRUE, scale=scale)
#   }
# 
# 
#   return(score)
# }


# 
res <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9")
# Define UI for application that draws a histogram
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      # Application title
      titlePanel("Pathway highlight"),
      
      textAreaInput("genes", "Add genes of interest", value = ""),
      actionButton(inputId = "addAGIs", label = "Submit genes"),
      radioButtons("path", "Choose a pathway from the list instead?",
                   choices = c("No" = "no",
                               "Yes" = "yes"),
                   selected = "no"),
      selectizeInput(
        inputId = "searchInput",
        label = "Input pathway",
        choices = NULL,
        multiple = FALSE,
        options = list(
          placeholder = 'Type to search...',
          openOnFocus = TRUE,
          autocomplete = TRUE,
          create = FALSE
        )
      ),
      radioButtons("split", "Split",
                   choices = c("No" = "no",
                               "By clusters" = "yes",
                               "By sample" = "yes_sam"),
                   selected = "no"),
      radioButtons("umap", "Dimentional reduction",
                   choices = c("UMAP" = "umap",
                               "t-SNE" = "tsne"),
                   selected = "umap"),
      radioButtons("norm", "Normalization", choices = c("Log-normalize"="log", 
                                                        "Integrated"="scale", 
                                                        "Z-score"="zsc",
                                                        "Raw counts" = "raw"),
                   selected =  "raw"),
      selectInput("residual", "Select residual for UMAP", res,
                  selected = "0.1"),
      radioButtons("add_pval", "Add P-value to violin plot",
                   choices = c("No" = "no",
                               "Yes" = "yes"),
                   selected = "no"),
      radioButtons("chg_shov", "Combine samples on violin plot",
                   choices = c("No" = "no",
                               "Yes" = "yes"),
                   selected = "no"),
      
      downloadButton('downloadmark', 'Download Plot'),
      downloadButton('downloadhist', 'Download Histogram'),
      downloadButton('downloadviol', 'Download Violin')
    ),
    
    
    
    mainPanel(
      
      mainPanel(
        plotly::plotlyOutput("myPlot", width = "140%", height = "600px"),
        plotly::plotlyOutput("myHist", width = "140%", height = "400px"),
        plotOutput("myViolin", width = "140%", height = "400px")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, server) {
  my_res <- reactive({
    res_us <- input$residual
    return(res_us)})
  
  my_dim <- reactive({
    res_us <- input$umap
    return(dim_us)})
  
  
  my_genes <- reactiveVal(character(0))
  
  # Observe the click on "Submit genes" button
  observeEvent(input$addAGIs, {
    data_genes <- strsplit(input$genes, "\n")[[1]]
    data_genes <- toupper(data_genes)
    data_genes <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", data_genes, perl=TRUE)
    my_genes(data_genes)  # Store the submitted genes in the reactive value
  })
  
  observeEvent(pathways, {
    updateSelectizeInput(
      inputId = "searchInput",
      choices = names(pathways),
      selected = NULL
    )
  })
  
  filteredData <- reactive({
    keyword <- input$searchInput
    if (!is.null(keyword) && keyword != "") {
      vec <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", as.vector(pathways[[keyword]]), perl=TRUE)  
    }
  })
  
  my_value2 <- reactive({
    if("log" %in% input$norm){
      DefaultAssay(whole.integrated) <- "RNA"
      count_matrix <- GetAssayData(whole.integrated, slot = "data")
    }
    if("zsc" %in% input$norm){
      DefaultAssay(whole.integrated) <- "RNA"
      count_matrix <- GetAssayData(whole.integrated, slot = "scale.data")
    }
    if("scale" %in% input$norm){
      DefaultAssay(whole.integrated) <- "integrated"
      count_matrix <- GetAssayData(whole.integrated, slot = 'data')
    }
    if("raw" %in% input$norm){
      DefaultAssay(whole.integrated) <- "RNA"
      count_matrix <- GetAssayData(whole.integrated, slot = 'counts')
    }
    
    # Subset the data frame to remove columns with any missing values
    return(count_matrix)
  })
  
  my_value <- function(){
    count_matrix <- as.data.frame(t(as.matrix(my_value2())))
    count_matrix$cells <- rownames(count_matrix)
    if (input$path=="no" & length(my_genes())<=1){
      count_matrix <- as.data.frame(count_matrix[,colnames(count_matrix)%in% my_genes()])
      whole_data <- cbind.data.frame(umap_pos, tsne_pos, count_matrix)
      colnames(whole_data)[ncol(whole_data)] <- "Expression"
      whole_data <- full_join(whole_data, data_merg, by = "cells")
    }
    else{
      count_matrix <- as.data.frame(count_matrix[,colnames(count_matrix)%in% "Acp5"])
      whole_data <- cbind.data.frame(umap_pos, tsne_pos, count_matrix)
      if (input$path=="no" & length(my_genes())>1){
        whole_data$Expression <- my_plot_path_list()[,1]
      }
      if (input$path=="yes"){
        whole_data$Expression <-my_plot_path()[,1]
      }
      colnames(whole_data)[ncol(whole_data)] <- "Expression"
      whole_data <- full_join(whole_data, data_merg, by = "cells")
    }
    rm(count_matrix)
    gc()
    return(whole_data)
  }
  
  
  my_plot_path<- function(){
    if("scale" %in% input$norm){
      DefaultAssay(whole.integrated) <- "integrated"
    }
    else{
      DefaultAssay(whole.integrated) <- "RNA"
    }
    count_matrix <- GetAssay(object = whole.integrated, assay=DefaultAssay(whole.integrated))
    plot_fgsea_path <- addGesecaScores(list(pathway=filteredData()), count_matrix,
                                       scale=TRUE)
    return(plot_fgsea_path)
  }
  
  
  my_plot_path_list <- function(){
    vec <- as.vector(my_genes())
    if("scale" %in% input$norm){
      DefaultAssay(whole.integrated) <- "integrated"
    }
    else{
      DefaultAssay(whole.integrated) <- "RNA"
    }
    count_matrix <- GetAssay(object = whole.integrated, assay=DefaultAssay(whole.integrated))
    plot_fgsea <- addGesecaScores(list(pathway=vec),count_matrix, scale=TRUE) 
    return(plot_fgsea)
  }
  
  reduction <- list("umap" = c("UMAP_1", "UMAP_2"), 
                    "tsne" = c("tSNE_1", "t_SNE2"))
  
  t <- list(
    size = 14)
  
  my_data <- function(){
    whole_data <- my_value()
    colnames(whole_data)[which(colnames(whole_data)==paste0("integrated_snn_res.", my_res()))] <- "clusters"
    if (input$umap == "umap"){
      whole_data$Dim1 <- whole_data$UMAP_1
      whole_data$Dim2 <- whole_data$UMAP_2
    }
    if (input$umap == "tsne"){
      whole_data$Dim1 <- whole_data$tSNE_1
      whole_data$Dim2 <- whole_data$tSNE_2
    }
    if (input$path=="no" & length(my_genes())<=1){
      whole_data <- whole_data
    }
    if (input$path=="no" & length(my_genes())>1){
      whole_data$Expression <- my_plot_path_list()[,1]
    }
    if (input$path=="yes"){
      whole_data$Expression <-my_plot_path()[,1]
    }
    
    centers <- whole_data %>%
      group_by(clusters) %>%
      summarize(center_x = mean(Dim1),
                center_y = mean(Dim2)) %>%
      ungroup()
    
    return(list(whole_data, centers))
  }
  
  
  my_plot <- function(){
    whole_data <- my_data()[[1]]
    centers <- my_data()[[2]]
    plot <- whole_data%>%ggplot(aes(x = Dim1, 
                                    y = Dim2))+
      geom_point(aes(color = Expression, text = paste("Expression:", Expression,
                                                      "\n", "UMAP 1:", UMAP_1,
                                                      "\n", "UMAP 2:", UMAP_2,
                                                      "\n", "t-SNE 1:", tSNE_1,
                                                      "\n", "t-SNE 2:", tSNE_2,
                                                      "\n", "Cluster:", clusters,
                                                      "\n", "Sample:", sample)), size = 1)+
      theme_classic()+geom_text(data = centers, aes(x = center_x, 
                                                    y = center_y, 
                                                    label = clusters), 
                                size = 5)
    if (input$split=="no") {
      plot <- plot
    }
    
    if (input$split=="yes"){
      plot <- plot+facet_wrap(.~clusters)+coord_fixed()
    }
    if (input$split=="yes_sam"){
      plot <- plot+facet_wrap(.~sample, ncol = 2)+coord_fixed()
    }
    if (input$path=="no" & length(my_genes())<=1){ 
      plot <- plot+
        scale_color_gradient2(low = "darkblue", mid = "lightgray", 
                              high = "#E41A1C", midpoint = 0)+
        theme(text = element_text(size = 15))+
        labs(color = "Expression")
    }
    else{
      plot <- plot+
        scale_color_gradient2(low = "darkblue", mid = "lightgray", 
                              high = "#E41A1C", midpoint = 0,
                              limits=c(-3, 3), breaks=c(-3, 0, 3))+
        theme(text = element_text(size = 15))+
        labs(color = "Expression")
      
    }
    if (input$umap == "umap"){
      plot <- plot+xlab("UMAP 1")+ylab("UMAP 2")
    }
    if (input$umap == "tsne"){
      plot <- plot+xlab("t-SNE 1")+ylab("t-SNE 2")
    }
    return(plotly::ggplotly(plot, 
                            tooltip = c("text")))
  }
  
  
  my_plot2 <- function(){
    whole_data <- my_data()[[1]]
    centers <- my_data()[[2]]
    
    plot1 <- whole_data%>%ggplot(aes(x = Dim1, y = Dim2))+
      geom_point(aes(color = clusters))+
      theme_classic()+
      theme(legend.position = "none", text = element_text(size = 20))+
      geom_text(data = centers, aes(x = center_x,
                                    y = center_y,
                                    label = clusters),
                size = 10)
    
    if(input$split=="yes_sam"){
      df_plot <- whole_data %>%
        group_by(sample, clusters) %>%
        summarise(sum_value = n()) %>%
        ungroup()%>%
        group_by(sample)%>%
        summarise(perc = sum_value/sum(sum_value)*100, clusters = clusters)%>%
        ungroup()
      
      plot2 <- df_plot%>%ggplot(aes(x=clusters, y = perc))+
        geom_bar(stat="identity") +
        ylab("Percentage of cells, %")+
        theme_classic()+
        theme(text = element_text(size = 20))+
        xlab("Clusters")+
        facet_wrap(.~sample)
      
      plot1 <- plot1+
        facet_wrap(.~sample)
    }
    
    else{
      plot2 <- whole_data%>%ggplot(aes(x=clusters))+
        geom_bar(aes(y = ..count.. / sum(..count..)), stat = "count") +
        scale_y_continuous(labels = scales::percent_format())+
        ylab("Percentage of cells")+
        theme_classic()+
        theme(text = element_text(size = 20))+
        xlab("Clusters")
      plot1 <- plot1
    }
    if(input$umap=="umap"){
      plot1 <- plot1+xlab("UMAP 1")+ylab("UMAP 2")
    }
    if(input$umap=="tsne"){
      plot1 <- plot1+xlab("t-SNE 1")+ylab("t-SNE 2")
    }
    plot2 <- plotly::ggplotly(plot2, tooltip = "all")
    plot1 <- plotly::ggplotly(plot1, tooltip = "all")
    return(plotly::subplot(plot2, plot1))
  }
  
  my_plot3 <- function(){
    whole_data <- my_data()[[1]]
    plot_vln <- whole_data%>%ggplot(aes(y = Expression, x = clusters))+
      geom_violin(aes(fill = clusters), width=0.5)+
      geom_boxplot(aes(fill = clusters), width=0.05, outlier.shape = NA)+
      theme_classic()+
      xlab("Clusters")+
      ylab("Expression")+
      theme(legend.position = "none", text = element_text(size = 20))
    
    if(input$split=="yes_sam"){
      if(input$chg_shov=="yes"){
        plot_vln <- whole_data%>%ggplot(aes(y = Expression, x = clusters))+
          geom_violin_pattern(aes(fill = clusters,
                                  pattern = sample),
                              pattern_angle = 45,
                              colour  = 'black',
                              pattern_fill = "black",
                              position = position_dodge(0.75))+
          geom_boxplot_pattern(aes(fill = clusters,
                                   pattern = sample),
                               pattern_angle = 45,
                               colour  = 'black',
                               pattern_fill = "black",
                               outlier.shape = NA,
                               width = 0.05,
                               position = position_dodge(0.75))+
          scale_pattern_manual(values = c("none","stripe"))+
          theme_classic()+
          xlab("Clusters")+
          ylab("Expression")+
          theme(text = element_text(size = 20), legend.title=element_text(size=20))+
          guides(fill_alpha = guide_legend(title="Sample"), fill = "none")+
          theme(legend.key.size = unit(3, "lines"),
                legend.spacing = unit(0.5, "lines"),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16))+
          ylim(c(min(whole_data$Expression), max(whole_data$Expression)+1))
        
        stat.test <- whole_data %>%
          group_by(clusters)%>%
          wilcox_test(Expression ~ sample) %>%
          adjust_pvalue(method = "bonferroni")%>%
          mutate(
            p.adj.signif = case_when(
              p.adj >= 0.05 ~ "ns",
              p.adj < 0.001 ~ "****",
              p.adj < 0.01 ~ "***",
              p.adj < 0.05 ~ "**",
              TRUE ~ "*"
            ))%>%
          add_xy_position(x = "supp")%>%
          dplyr::select(-c("xmin", "xmax"))
        
      }
      else{
        plot_vln <-  plot_vln+
          facet_wrap(.~sample)
        
        stat.test <- whole_data %>%
          group_by(sample)%>%
          wilcox_test(Expression ~ clusters, ref.group = "all") %>%
          adjust_pvalue(method = "bonferroni")%>%
          add_xy_position(x = "supp")%>%
          dplyr::select(-c("xmin", "xmax"))
      }
    }
    else{
      plot_vln <- plot_vln
      
      stat.test <- whole_data %>% wilcox_test(Expression ~ clusters, ref.group = "all") %>%
        adjust_pvalue(method = "bonferroni")%>%
        add_xy_position(x = "supp")%>%
        dplyr::select(-c("xmin", "xmax"))
      
    }
    if(input$add_pval=="yes"){
      if(input$split=="yes_sam"&input$chg_shov=="yes"){
        
        plot_vln <- plot_vln+
          geom_text(data = stat.test %>% filter(p.adj < 0.05),
                    aes(x = clusters, y = max(whole_data$Expression), label = p.adj.signif),
                    size = 6, vjust = -1)
      }
      else{
        plot_vln <- plot_vln+
          stat_pvalue_manual(stat.test, label = "p.adj.signif",
                             y.position = "y.position")
      }
      return(plot_vln)
    }
    else{
      return(plot_vln)
    }
    gc()
  }
  
  output$myPlot <- plotly::renderPlotly({
    my_plot()
  })
  
  output$myHist <- plotly::renderPlotly({
    my_plot2()
  })
  
  output$myViolin <- renderPlot({
    print(my_plot3())
  })
  
  output$downloadmark <- downloadHandler(
    filename = "Markers.png",
    content = function(file) {
      png(file, width = 600, height=500)
      print(my_plot())
      dev.off()
    })  
  
  output$downloadhist <- downloadHandler(
    filename = "Histogram.png",
    content = function(file) {
      png(file, width = 800, height=500)
      print(my_plot2())
      dev.off()
    })
  
  output$downloadviol <- downloadHandler(
    filename = "Violin plot.png",
    content = function(file) {
      png(file, width = 600, height=500)
      print(my_plot3())
      dev.off()
    })
}

# Run the application

#Add 
shinyApp(ui = ui, server = server)