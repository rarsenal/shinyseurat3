library(shinyjs)

SeuratObject <- RenameIdents(
  object = SeuratObject,
  "FIB-1" = "Fibroblasts-1",
  "FIB-2" = "Fibroblasts-2",
  "FIB-3" = "Fibroblasts-3",
  "FIB-4" = "Fibroblasts-4",
  "FIB-5" = "Fibroblasts-5",
  "ENDO" = "Endothelial cells",
  "MYL" = "Myeloid cells",
  "DEN" = "Dendritic cells",
  "TCELL" = "T cells",
  "BCELL" = "B cells",
  "LYME" = "Lymphatic cells",
  "SCH" = "Schwann cells",
  "RBC" = "Red blood cells"
)

levels(x = SeuratObject) <- c(
  "Fibroblasts-1",
  "Fibroblasts-2",
  "Fibroblasts-3",
  "Fibroblasts-4",
  "Fibroblasts-5",
  "Endothelial cells",
  "Lymphatic cells",
  "Myeloid cells",
  "Dendritic cells",
  "T cells",
  "B cells",
  "Schwann cells",
  "Red blood cells"
)


function(input, output, session) {
  output$Draw <- renderPlot({staticImage})
  output$Logo <- renderPlot({logoImage})
  
  # Display Main Reduction Plots as overview toggled options for display legends/labels and setting size. Default displays demo dataset

  output$Main <- renderPlot({
    
    pt1 <- DimPlot(get("SeuratObject"), reduction = "tsne", label = FALSE, pt.size = 0.1, label.size = 6) + theme(aspect.ratio = 1, legend.position = "none") + NoAxes()

    #pt1 <- DimPlot(get("SeuratObject"), reduction = input$pt1Type, label = input$labelBoolean, pt.size = input$dotSize, label.size = input$labelSize) + theme(aspect.ratio = 1, legend.position = "none") + NoAxes()
    #pt2 <- staticImage
    ggsave(file="wound_tsne.png", plot=pt1, dpi=600)
    #pt2 <- DimPlot(get("SeuratObject"), reduction = input$pt2Type, label = input$labelBoolean, pt.size = input$dotSize, label.size = input$labelSize)
    #CombinePlots(plots = list(pt1, pt2), legend = input$legendBoolean)
    #cowplot::plot_grid(pt1, pt2, rel_widths = c(1, 2))
    #hoverdata <-FetchData(get("SeuratObject"), vars = c("ident"))

    pt1
    #pt2 <- HoverLocatorP(plot = pt1, information = hoverdata)
    #pt2 %>%
    #  config(displayModeBar = F) %>% layout(xaxis=list(fixedrange=TRUE))
    #pt2 %>%
    #  config(displayModeBar = F) %>% layout(height = "auto")
  })

  updateSelectizeInput(session, 'geneNames', choices = choices3, selected=NULL, server = TRUE)

  select_genes <- reactiveValues(
      genes = NULL
  )

  observeEvent(input$submit_loc, {
    if (is.null(input$geneNames)) {
	selected=sample(top100,3)
	updateSelectizeInput(session, 'geneNames', choices = choices3, selected=selected, server = TRUE)
        select_genes$genes <- input$geneNames
    }
    
    select_genes$genes <- input$geneNames

  })


  fPlot <- eventReactive(select_genes$genes, {
    
    plots <- FeaturePlot(get("SeuratObject"), reduction = "tsne", select_genes$genes, pt.size = 0.1, coord.fixed=TRUE, combine = FALSE)
    #plots <- FeaturePlot(get("SeuratObject"), reduction = input$pt1Type, input$geneNames, pt.size = input$dotSize, coord.fixed=FALSE, combine = FALSE)
    #plots <- lapply(X = plots, FUN = function(x) x + theme(aspect.ratio = 0.65, plot.title = element_text(size = 20, face="italic")) + NoAxes())
    #hoverdata <-FetchData(get("SeuratObject"), vars = c("ident"))
    plots <- lapply(X = plots, FUN = function(x) x + theme(plot.title = element_text(size = 20, face="italic")) + NoAxes())
    cowplot::plot_grid(plotlist = plots, ncol=3)
    #plots <- lapply(X = plots, FUN = function(x) HoverLocatorP(plot = x, information = hoverdata) 
    #p1 <- HoverLocatorP(plot = plots[[1]], information = hoverdata) %>% config(displayModeBar = F)
    #p2 <- HoverLocatorP(plot = plots[[2]], information = hoverdata) %>% config(displayModeBar = F)
    #p3 <- HoverLocatorP(plot = plots[[3]], information = hoverdata) %>% config(displayModeBar = F)
    #subplot(p1, p2, p3)
    #FeaturePlot(get("SeuratObject"), reduction = input$pt1Type, input$geneNames, pt.size = input$dotSize, coord.fixed=TRUE, ncol=3)
  })
  
  vPlot <- eventReactive(select_genes$genes, {

    Vplots <- VlnPlot(get("SeuratObject"), features = select_genes$genes, combine = FALSE, pt.size = 0)
    Vplots <- lapply(X = Vplots, FUN = function(x) x + theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(size = 20, face="italic")) + geom_jitter(colour = "#404040", alpha=0.4, size=1))
    cowplot::plot_grid(plotlist = Vplots, ncol=3)
  })
  
  rPlot <- eventReactive(select_genes$genes, {

    Rplots <- RidgePlot(get("SeuratObject"), features = select_genes$genes, combine = FALSE)
    Rplots <- lapply(X = Rplots, FUN = function(x) x + theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 20, face="italic")))
    cowplot::plot_grid(plotlist = Rplots, ncol=3)
  })
  
  hPlot <- eventReactive(select_genes$genes, { DoHeatmap(subset(get("SeuratObject"), downsample = 100), features = select_genes$genes, size = 5) + theme(legend.position = "bottom", axis.text.y = element_text(size=12, face="bold.italic"))+ 
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  })
  

  output$featureText <- renderText({'Start by typing gene names above. Note that only significantly expressed genes are searchable.'})

  output$featurePlot <- renderPlot(fPlot(), height = function() {
                   session$clientData$output_featurePlot_width/3 })

  output$feature1 <- renderUI(
  {
     if (is.null(select_genes$genes)) {
         verbatimTextOutput("featureText")
     } else {
         plotOutput("featurePlot", height="auto") %>% withSpinner(color="#3399cc")
     }
  })

  output$vlnPlot <- renderPlot(vPlot(),height = function() {
    session$clientData$output_vlnPlot_width/3})
  output$ridgePlot <- renderPlot(rPlot(),height = function() {
    session$clientData$output_ridgePlot_width/3})
  output$heatPlot <- renderPlot(hPlot())
}
