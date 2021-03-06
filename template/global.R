library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)
library(shinyWidgets)
library(shinythemes)

SeuratObject <-	readRDS(file="seurat3.rds")
choices2 = rownames(SeuratObject)
reductionChoices = names(SeuratObject@reductions)
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst")
top100 = head(VariableFeatures(SeuratObject), 100)
con = file("title.txt", "r")
titlename = readLines(con, n = 1)
close(con)
