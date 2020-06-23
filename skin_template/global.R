library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)
library(shinyWidgets)
library(shinythemes)
library(ggplot2)
library(plotly)

SeuratObject <-	readRDS(file="seurat3.rds")
choices2 = rownames(SeuratObject)
reductionChoices = names(SeuratObject@reductions)
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst")
choices3 <- VariableFeatures(SeuratObject)
top100 = head(choices3, 100)
con = file("title.txt", "r")
titlename = readLines(con, n = 1)
close(con)
#staticImage = cowplot::ggdraw() + cowplot::draw_image("skin_demo.svg", scale = 0.8)
#logoImage = cowplot::ggdraw() + cowplot::draw_image("logo.png", scale = 1)



HoverLocatorP <- function(
  plot,
  information = NULL,
  dark.theme = FALSE,
  ...
) {
  #   Use GGpointToBase because we already have ggplot objects
  #   with colors (which are annoying in plotly)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE)
  data <- ggplot_build(plot = plot)$plot$data
  rownames(x = plot.build) <- rownames(x = data)
  #   Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  #   Add the names we're looking for (eg. cell name, gene name)
  if (is.null(x = information)) {
    plot.build$feature <- rownames(x = data)
  } else {
    info <- apply(
      X = information,
      MARGIN = 1,
      FUN = function(x, names) {
        return(paste0("Cell type", ': ', x, collapse = '<br>'))
      },
      names = colnames(x = information)
    )
    data.info <- data.frame(
      #feature = paste(rownames(x = information), info, sep = '<br>'),
      feature = info,
      row.names = rownames(x = information)
    )
    plot.build <- merge(x = plot.build, y = data.info, by = 0)
  }
  #   Set up axis labels here
  #   Also, a bunch of stuff to get axis lines done properly
  xaxis <- list(
    title = "",
    showgrid = FALSE,
    zeroline = FALSE,
    showticklabels = FALSE,
    showline = FALSE
  )
  yaxis <- list(
    title = "",
    showgrid = FALSE,
    zeroline = FALSE,
    showticklabels = FALSE,
    showline = FALSE
  )
  #   Check for dark theme
  if (dark.theme) {
    title <- list(color = 'white')
    xaxis <- c(xaxis, color = 'white')
    yaxis <- c(yaxis, color = 'white')
    plotbg <- 'black'
  } else {
    title = list(color = 'black')
    plotbg = 'white'
  }
  #   The `~' means pull from the data passed (this is why we reset the names)
  #   Use I() to get plotly to accept the colors from the data as is
  #   Set hoverinfo to 'text' to override the default hover information
  #   rather than append to it
  p <- plotly::layout(
    p = plot_ly(
      data = plot.build,
      x = ~x,
      y = ~y,
      type = 'scatter',
      mode = 'markers',
      color = ~I(color),
      hoverinfo = 'text',
      text = ~feature,
      marker = list(size = 2)
    ),
    xaxis = xaxis,
    yaxis = yaxis,
    title = plot$labels$title,
    titlefont = title,
    paper_bgcolor = plotbg,
    plot_bgcolor = plotbg,
    ...
  )
  # add labels
  label.layer <- which(x = sapply(
    X = plot$layers,
    FUN = function(x) class(x$geom)[1] == "GeomText")
  )
  if (length(x = label.layer) == 1) {
    p <- plotly::add_annotations(
      p = p,
      x = plot$layers[[label.layer]]$data[, 1],
      y = plot$layers[[label.layer]]$data[, 2],
      xref = "x",
      yref = "y",
      text = plot$layers[[label.layer]]$data[, 3],
      xanchor = 'right',
      showarrow = FALSE,
      font = list(size = plot$layers[[label.layer]]$aes_params$size * 4)
    )
  }
  return(p)
}


GGpointToBase <- function(plot, do.plot = TRUE, ...) {
  plot.build <- ggplot_build(plot = plot)
  cols <- c('x', 'y', 'colour', 'shape', 'size')
  build.use <- which(x = vapply(
    X = plot.build$data,
    FUN = function(dat) {
      return(all(cols %in% colnames(x = dat)))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  if (length(x = build.use) == 0) {
    stop("GGpointToBase only works on geom_point ggplot objects")
  }
  build.data <- plot.build$data[[min(build.use)]]
  plot.data <- build.data[, cols]
  names(x = plot.data) <- c(
    plot.build$plot$labels$x,
    plot.build$plot$labels$y,
    'color',
    'pch',
    'cex'
  )
  if (do.plot) {
    PlotBuild(data = plot.data, ...)
  }
  return(plot.data)
}
