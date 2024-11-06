# Cargar las librerías necesarias
library(SummarizedExperiment)
library(readxl)
library(ggplot2)
library(limma)

# Especificar la ruta del archivo en el repositorio
file_path <- "path/to/TIO2+PTYR-human-MSS+MSIvsPD.XLSX"

# Extraer datos desde el Excel
MSdata <- read_excel("TIO2+PTYR-human-MSS+MSIvsPD.XLSX", sheet = 1)

# Seleccionar data de las columns de interes, columna 1 para los fofolipidos y columnas 5 a 16 para los datos de abundancia5 to 16
text_data <- MSdata[, c(1, 5:16)]

# Escribe los datos en de la tabla en formato texto
write.table(text_data, file = "Text_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Seleccionar las columnas 5 a 16 como los datos de ensayo
assay_data <- as.matrix(MSdata[, 5:16])

# Poner nombres de las filas de los datos de ensayo con la columna de SequenceModifications
rownames(assay_data) <- MSdata$SequenceModifications

# Extraer el nombre muestra y el fenotipo a partir de los nombres de las columnas
# Extraer el ID de muestra
sample_ids <- gsub("_.*", "", colnames(assay_data)) 

# Asignar el fenotipo basado en MSS o PD en los nombres
phenotypes <- ifelse(grepl("MSS", colnames(assay_data)), "MSS", "PD")  

# Crear un data frame para colData
col_data <- data.frame(SampleID = sample_ids,Phenotype = phenotypes)

# Poner los nombre de fila a los datos de la s columnas
rownames(col_data) <- colnames(assay_data) 

# Anhadir metadatos
experiment_metadata <- list(conditions = c("MSS", "PD"), description = "El conjunto de datos adjunto se obtuvo de un experimento de fosfoproteómica realizado para analizar (3 + 3) modelos PDX de dos subtipos diferentes utilizando muestras enriquecidas de fosfopeptidos. Se realizó un análisis LC-MS de 2 duplicados técnicos en cada muestra. El conjunto de resultados consistió en abundancias normalizadas de señales de MS para aproximadamente 1400 fosfopeptidos. Objetivo del análisis: buscar fosfopeptidos que permitan diferenciar los dos grupos de tumores. Esto debe realizarse con análisis estadístico y visualización. Los datos se han proporcionado en un archivo de Excel: TIO2+PTYR-human-MSS+MSIvsPD.XLSX.")

# Crear el objeto SummarizedExperiment
se <- SummarizedExperiment( assays = list (counts = assay_data),colData = col_data, metadata = experiment_metadata)

#Guardar el objeto SummarizedExperiment en un contenedor de formato .rda
save(se,file = "summarized_experiment.rda")

# Ver el objeto SummarizedExperiment
print (se)

# Análisis 

# Extraemos datos de se
counts <- assay(se, "counts")               
phenotypes <- colData(se)$Phenotype          
# Dividir los datos en grupos MSS y PD, se consideran las replicas technicas como muestras
mss_counts <- counts[, phenotypes == "MSS"]
pd_counts <- counts[, phenotypes == "PD"]

# comparar con un test de student las medias de los fenotipos MDD y PD de cada fosfopéptido 
p_values <- apply(counts, 1, function(x) {t.test(x[phenotypes == "MSS"], x[phenotypes == "PD"])$p.value})

# Ajustar los valores p utilizando el método de la Tasa de Descubrimiento Falso (FDR)
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Calcular el cambio doblaje logarítmico entre los grupos
logCD <- rowMeans(mss_counts) - rowMeans(pd_counts)

# Compilar los resultados en un data frame
comparison_results <- data.frame(Phosphopeptide = rownames(counts),  logCD = logCD,  P.Value = p_values,  adj.P.Val = adjusted_p_values)

# Ver fosfopéptidos significativos valor p ajustado menor que 0.05
sig_fosfo <- comparison_results[comparison_results$adj.P.Val < 0.05,]

# Crear un gráfico de volcano para visualizar los fosfopéptidos significativos
ggplot(comparison_results, aes(x = logCD, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +   # Colorear los puntos según significancia
  scale_color_manual(values = c("black", "green")) +           # Negro para no significativos verde para significativos
  labs(title = "Gráfico de Volcano Diferencias de fosfopéptidos") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Línea del limite 0.05
  guides(color = guide_legend(title = "Significativo"))

# Filtrar los fosfopéptidos significativos (valor p ajustado < 0.05)
sig_fosfo <- comparison_results[comparison_results$adj.P.Val < 0.05, ]

# Extraer los datos de fosfopéptidos significativamente deferented para ambos grupos
significant_counts <- counts[rownames(counts) %in% sig_fosfo$Phosphopeptide, ]

# Calcular las medias (MSS y PD)
mss_means <- rowMeans(significant_counts[, phenotypes == "MSS"])
pd_means <- rowMeans(significant_counts[, phenotypes == "PD"])

# Crear un data frame 
plot_data <- data.frame(
  Phosphopeptide = rep(rownames(significant_counts), 2),
  Group = rep(c("MSS", "PD"), each = nrow(significant_counts)),
  Mean_Counts = c(mss_means, pd_means))

# Gráfico de barras 
ggplot(plot_data, aes(x = Phosphopeptide, y = Mean_Counts, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Media de counts de Fosfopeptidos MSS vs PD",
       x = "Fosfopeptidos",
       y = "Counts") +
  theme_minimal() +
  theme(axis.text.x = element_blank())  +
    scale_fill_manual(values = c("MSS" = "blue", "PD" = "orange"))

# Imprimir lista de fosfopeptidos significativos. 

print(sig_fosfo)

# Realizar PCA

# Reordenar y limpiar datos para que no de errores
counts_t <- t(counts)
counts_t <- counts_t[, apply(counts_t, 2, var) != 0]

# Realizar PCA en los datos 
pca_result <- prcomp(counts_t, scale. = TRUE)
print(pca_result)

# Crear un data frame con los resultados de la PCA 
pca_data <- as.data.frame(pca_result$x)
pca_data$Phenotype <- phenotypes  

# Grafica de los resultados de la PCA (PC1 vs PC2) por fenotipo
ggplot(pca_data, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size = 3) +
  labs(title = "PCA de Recuentos de fosfopéptidos",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("MSS" = "blue", "PD" = "orange")) +
  theme(legend.position = "top")

# Extraer los loadings de cada fosfopéptido en PC1
loadings <- pca_result$rotation

# Ordenar las cargas absolutas de PC1 en orden descendente para identificar los más importantes
important_features_pc1 <- sort(abs(loadings[, "PC1"]), decreasing = TRUE)
head(important_features_pc1)

