
### Script para detectar CNVs con la herramienta CODEX2


## Se cargan los paquetes necesarios

library(CODEX2)

library(CODEX)

library(dplyr)


## Se cargan los archivos BAM y BED y se obtiene el nombre de las muestras

bamdir  <- list.files("C:/Users/Usuario/Desktop/TFM/BAMfiles/", pattern= ".bam$", full.names =TRUE)

sampname <- as.matrix(unlist(strsplit(bamdir,"\\.bam")))

bedFile  <- "C:/Users/Usuario/Desktop/TFM/archivo_bed.bed"


## Se crea un bucle para ejecutar el código con todos los cromosomas presentes en el archivo BED

for(chr in c(1:19,21,22)){

bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "TFM", chr)

bamdir <- bambedObj$bamdir

sampname <- bambedObj$sampname

ref <- bambedObj$ref

projectname <- bambedObj$projectname

chr <- bambedObj$chr


# Se analiza la cobertura y el contenido GC  

coverageObj <- getcoverage(bambedObj, mapqthres = 20)

Y <- coverageObj$Y

readlength <- coverageObj$readlength

gc <- getgc(chr, ref)

mapp <- getmapp(chr, ref)


# Se realiza un control de calidad 

qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000), length_thresh = c(20, 2000), mapp_thresh = 0.9, gc_thresh = c(20, 80))

Y_qc <- qcObj$Y_qc

sampname_qc <- qcObj$sampname_qc

gc_qc <- qcObj$gc_qc

mapp_qc <- qcObj$mapp_qc

ref_qc <- qcObj$ref_qc

qcmat <- qcObj$qcmat


# Se normalizan los datos

normObj <- normalize(Y_qc, gc_qc, K = 1:8)

Yhat <- normObj$Yhat

AIC <- normObj$AIC

BIC <- normObj$BIC

RSS <- normObj$RSS

K <- normObj$K


# Se determinan las CNVs presentes en las muestras

optK = K[which.max(BIC)]

finalcall <- segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc, ref_qc, chr, lmax = 200, mode = "fraction")


# Los resultados de las CNVs detectadas se guardan en archivos .txt, se crea un archivo de texto
# para cada cromosoma

write.table(finalcall, file = paste(projectname, '_', chr, '_', optK, '_CODEX_frac.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)

}

