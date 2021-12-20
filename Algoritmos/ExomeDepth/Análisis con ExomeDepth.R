
### Script para detectar CNVs con la herramienta ExomeDepth


# Se cargan los paquetes necesarios

library(ExomeDepth)

library(dplyr)


# Se cargan los archivos que se requieren como input:
# el archivo BED, los ficheros BAM y el genoma de referencia en formato FASTA

BED_file <- "C:/Users/Usuario/Desktop/TFM/archivo_bed.bed"

archivo_bed <- read.table(BED_file, header = F)

BAMFiles <- list.files("C:/Users/Usuario/Desktop/TFM/BAMfiles/", pattern= ".bam$", full.names =TRUE)

BAMFiles_cntrl <- list.files("C:/Users/Usuario/Desktop/TFM/BAMfilesCntrl", pattern= ".bam$", full.names =TRUE)

hg19 <- "C:/Users/Usuario/Desktop/TFM/Homo_sapiens.GRCh37.dna.primary_assembly.fa"


# Se hace un recuento de lecturas a partir de los archivos BAM de las muestras control

# En la ejecución del código se detecta un error que impide correr el código si se tienen
# en cuenta los cromosomas 5-9, X e Y. Finalmente, se decide no utilizar esos cromosomas

ExomeCount_Cntrl <- getBamCounts(bed.frame = archivo_bed %>% filter(!V1%in%c("5","6","7","8","9","X","Y")),
                                 bam.files = BAMFiles_cntrl,
                                 include.chr = FALSE,
                                 referenceFasta = hg19)


# Se hace un recuento de lecturas a partir de los archivos BAM del resto de muestras

# En la ejecución del código se detecta el mismo error que en las muestras control, este error impide 
# correr el código si se tienen en cuenta los cromosomas 5-9, X e Y. Finalmente, se decide no utilizar esos cromosomas

ExomeCount <- getBamCounts(bed.frame = archivo_bed %>% filter(!V1%in%c("5","6","7","8","9","X","Y")),
                           bam.files = BAMFiles,
                           include.chr = FALSE,
                           referenceFasta = hg19)


# Se convierten ambos recuentos de lectura a data frames

ExomeCount.dafr_Cntrl <- as(ExomeCount_Cntrl, 'data.frame')

print(head(ExomeCount.dafr_Cntrl))

ExomeCount.dafr <- as(ExomeCount, 'data.frame')

print(head(ExomeCount.dafr))


# Se crea el set de referencia con las muestras control

my.ref.samples <-c ("X17299_alin.bam",
                    "X17300_alin.bam",
                    "X17312_alin.bam",
                    "X17327_alin.bam",
                    "X17356_alin.bam",
                    "X17373_alin.bam",
                    "X17399_alin.bam",
                    "X17401_alin.bam")

my.reference.set <- as.matrix(ExomeCount.dafr_Cntrl[, my.ref.samples])

ExomeCount_mat <- as.matrix(ExomeCount[,-c(1:5)])


# Se crea un bucle para analizar las muestras no-control y detectar las CNVs
# Los resultados de las CNVs detectadas se guardan en archivos .txt, se crea un archivo de texto
# para cada muestra

estudio_nsamples <- ncol(ExomeCount)-5

estudio_samples <-  colnames(ExomeCount)[-c(1:5)]

i<-1

for (i in 1:estudio_nsamples) {
  
  samplename <- estudio_samples[i]
  
  print(samplename)
  
  my.choice_1 <- select.reference.set (test.counts = ExomeCount[,i+5],
                                           reference.counts = my.reference.set,
                                           bin.length = (ExomeCount_Cntrl$end - ExomeCount_Cntrl$start)/1000, 
                                           n.bins.reduced = 10000)
  
  my.matrix <- as.matrix( ExomeCount.dafr_Cntrl[, my.choice_1$reference.choice]) 
  
  my.reference.selected <- apply(X = my.matrix,MAR = 1, FUN = sum)
  
  all.exons_1 <- new('ExomeDepth',
                         test = ExomeCount_mat[,i],
                         reference = my.reference.selected, 
                         formula = 'cbind(test, reference) ~ 1')
  
  show(all.exons_1)
  
  all_exons_CNVs <- CallCNVs(x = all.exons_1,
                             transition.probability = 10^-3,
                             chromosome = ExomeCount.dafr_Cntrl$chromosome,
                             start = ExomeCount.dafr_Cntrl$start,
                             end = ExomeCount.dafr_Cntrl$end,
                             name = ExomeCount.dafr_Cntrl$exon)
  
  head(all_exons_CNVs@CNV.calls)
  
  if(nrow(all_exons_CNVs@CNV.calls)>0){
    
    output.file <- paste0(i ,".cnv.txt")
    
    write.table(file = file.path(output.file),x = cbind(i,all_exons_CNVs@CNV.calls),row.names = FALSE,quote=F,sep="\t")
  }
}





