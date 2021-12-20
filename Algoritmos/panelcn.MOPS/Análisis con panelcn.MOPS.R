
### Script para detectar CNVs con la herramienta panelcn.MOPS


# Se cargan los paquetes necesarios

library(panelcn.mops)

library(dplyr)

library(writexl)


# Se lee el archivo BED con las coordenadas de los exones

BED_file <- "C:/Users/Usuario/Desktop/TFM/archivo_bed.bed"

countWindows <- getWindows(BED_file) 


# Se hace un recuento de lecturas a partir de los archivos BAM

BAMFiles <- list.files("C:/Users/Usuario/Desktop/TFM/BAMfiles/", pattern= ".bam$", full.names =TRUE)

BAMFiles_cntrl <- list.files("C:/Users/Usuario/Desktop/TFM/BAMfilesCntrl", pattern= ".bam$", full.names =TRUE)

test <- countBamListInGRanges(countWindows = countWindows,
                              bam.files = BAMFiles, read.width = 101)

control <- countBamListInGRanges(countWindows = countWindows,
                              bam.files = BAMFiles_cntrl, read.width = 101)


# Se corre el algoritmo y se obtiene una lista de tablas con los resultados

XandCB <- test

elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
                                 elementMetadata(control))


resultlist <- runPanelcnMops(XandCB,
                             1:ncol(elementMetadata(test)),
                             countWindows = countWindows,
                             I= c(0.025, 0.57, 1, 1.46, 2), normType = "quant",
                             sizeFactor = "quant", qu = 0.25, quSizeFactor = 0.75, norm = 1,
                             priorImpact = 1, minMedianRC = 30, maxControls = 25)

sampleNames <- colnames(elementMetadata(test))

resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                                 countWindows = countWindows,
                                 sampleNames = sampleNames)


# Se crean dos Excel
# El Excel resultados.panelcn.originales contiene todos los registros que se generan en el análisis
# El Excel resultados.filtrados contiene solamente aquellos registros que se corresponden con delecciones y duplicaciones (es decir, las CNVs detectadas)

tablas.filtradas <- lapply(resulttable,FUN=function(df){
    df <- 
      df %>% 
      filter(lowQual!="lowQual"&CN!="CN2")
    df
  })

nombres <- sapply(resulttable,FUN=function(x) unique(x[["Sample"]]))

names(resulttable) <- nombres

names(tablas.filtradas) <- nombres

writexl::write_xlsx(resulttable,path = "resultados.panelcn.originales.xlsx")

writexl::write_xlsx(tablas.filtradas,path="resultados.filtrados.xlsx")

