function(input, output, session) {
  
  # Update the Input field with 3 randomly selected genes from a list of 100 top expressing genes (see global.R for top100 and choices2 definitions)
  updateSelectizeInput(session, 'geneNames', choices = choices2, selected=sample(top100,3), server = TRUE)


  # Display Main Reduction Plots as overview toggled options for display legends/labels and setting size. Default displays demo dataset
  output$Main <- renderPlot({
    
    pt1 <- DimPlot(get("SeuratObject"), reduction = input$pt1Type, label = input$labelBoolean, pt.size = input$dotSize, label.size = input$labelSize)
    pt2 <- DimPlot(get("SeuratObject"), reduction = input$pt2Type, label = input$labelBoolean, pt.size = input$dotSize, label.size = input$labelSize)
    CombinePlots(plots = list(pt1, pt2), legend = input$legendBoolean)
    
  })
  
  # FeaturePlot #1     
  fPlot <- eventReactive(input$geneNames, {
    FeaturePlot(get("SeuratObject"), reduction = input$pt1Type, input$geneNames,pt.size = input$dotSize,coord.fixed=TRUE,ncol=3)
  })
  
  # FeaturePlot #2
  fPlot2 <- eventReactive(input$geneNames, {
    FeaturePlot(get("SeuratObject"), reduction = input$pt2Type, input$geneNames,pt.size = input$dotSize,coord.fixed=TRUE,ncol=3)
  })
  
  # ViolnPlot
  vPlot <- eventReactive(input$geneNames, {
    VlnPlot(get("SeuratObject"), features = input$geneNames)
  })
  
  # RidgePlot
  rPlot <- eventReactive(input$geneNames, {
    RidgePlot(get("SeuratObject"), features = input$geneNames)
  })
  
  # DotPlot
  dPlot <- eventReactive(input$geneNames, {DotPlot(get("SeuratObject"), features = input$geneNames)})
  
  # HeatMapPlot
  hPlot <- eventReactive(input$geneNames, {DoHeatmap(subset(get("SeuratObject"), downsample = 100), features = input$geneNames, size = 3)})
  
  output$feature1 <- renderPlot(fPlot(),height = function() {
    session$clientData$output_feature1_width/3})
  output$feature2 <- renderPlot(fPlot2(),height = function() {
    session$clientData$output_feature2_width/3})
  output$vlnPlot <- renderPlot(vPlot())
  output$ridgePlot <- renderPlot(rPlot())
  output$dotPlot <- renderPlot(dPlot())
  output$heatPlot <- renderPlot(hPlot())
  
}
