
### Script de la modificación del archivo BED


# Se cargan los paquetes necesarios

library(dplyr)

library(stringr)


# Se lee el archivo BED con las coordenadas de los exones de los genes del panel y se muestran los primeros registros

file.bed <- read.delim("BED.txt",header = F)

head(file.bed)


# Se crea un objeto lista para separar la cadena de texto de la variable V4 en subcadenas que contienen la información que está delimitada por puntos

obj <- stringr::str_split(string = file.bed$V4,pattern = "\\.")

  
# Se unen las cuatro subcadenas en un data frame
  
obj2 <- do.call("rbind.data.frame",obj)

colnames(obj2) <- c("gen","chr","inicio","final")

  
# Se une el data frame con el fichero bed original

obj3 <- cbind.data.frame(file.bed,obj2)


# Se eliminan las letras "chr", se ordenan los registros por cromosoma y coordenada de inicio y fin y se eliminan los registros "rs" correspondientes a datos de SNPs

obj3 <- obj3 %>% 
        mutate(chr.number=gsub(pattern = "chr",
                           replacement = "",
                           chr)) %>% 
        mutate(chr.number2=ifelse(chr.number=="X","23",
                              ifelse(chr.number=="Y","24",chr.number)),
           chr.number2=as.numeric(paste0(chr.number2)),
           inicio=as.numeric(paste0(inicio)),
           final=as.numeric(paste0(final))) %>% 
        arrange(chr.number2,inicio) %>% 
        select(-chr.number2) %>% 
        filter(!grepl("rs",gen))
  

# Se guarda en un nuevo fichero bed la información con las columnas y el formato adecuado para las herramientas de detección de CNVs

write.table(obj3 %>% select(chr.number,inicio,final,gen), file = "archivo_bed.bed", sep="\t", quote = F, row.names = F, col.names = F)
  
