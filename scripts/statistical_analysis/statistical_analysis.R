# Cargar librerias
library(phyloseq)
library(vegan)
library(microbiome)
library(ComplexHeatmap)
library(circlize)

# Guardar variantes
otu_mat <- as.matrix(read.table("otu_table.tsv", sep = "\t"))
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
tax_mat <- as.matrix(read.table("tax_table_phyloseq.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = 1))
tax <- tax_table(tax_mat)
metadata <- (read.table("metadata.tsv", sep ="\t"))

# Manipulación metadatos
header.true <- function(df) {
  names(df) = as.character(unlist(df[1,]))
  df[-1,]
}
metadata <- header.true(metadata)

# Comprobación eliminación V1, V2, V3
head(metadata)

# Asignación nombre de las filas igual a la columna metadata sample ID
rownames(metadata) <- metadata$Seq_ID

# Comprobación
head(metadata)

# Creación objeto metadata en formato para phyloseq
metadata <- sample_data(metadata)

# Arbol
# tree = read_tree("tree.tree")
# tree

# Comprobación
sample_names(otu)
sample_names(tax)
sample_names(metadata)

# Creación objeto physeq
physeq = phyloseq(otu, tax, metadata)
physeq

# Comprobación primeras filas de los metadatos
head(meta(physeq))

# Cambio nombres otu
taxa_names(physeq)
#taxa_names(physeq) <- paste("ASV", seq(ntaxa(physeq)))

library(ggplot2)


# Gráfico para ver especies
plot_bar(physeq, fill = "Species") + geom_bar(aes(color = Species, fill = Species), stat = "identity", position = "stack")

# Gráfico para ver familias
colnames(tax)
plot_bar(physeq, fill = "Phylum")
# Histograma de total de lecturas

hist(sample_sums(physeq), main="Histogram: Read Counts", xlab= "Total Reads", border="blue", col="green", las=1, breaks=12)

# Filtro para quitar muestras que tengan menos de 1000 lecturas
physeq_reads_1k = prune_samples(sample_sums(physeq) > 1000, physeq)

# Calculo de la diversidad alfa
alpha <- microbiome::alpha
tab <- alpha(physeq, index= "all")
colnames(tab)
# Visualización
head(tab)

# Calculo de la diversity indicando solo el índice de Shannon
metadatos_diversidad <- meta(physeq)
metadatos_diversidad$Shannon <- microbiome::diversity(
  physeq,
  index = "Shannon"
)$shannon
# Visualización
metadatos_diversidad

# Visualización gráfica en función edad y diversidad
grafico_1 <- ggplot(metadatos_diversidad, aes(x= Age_Group, y = Shannon, fill = Age_Group, )) + geom_violin(trim=FALSE)
grafico_1

# Separación datos diversidad y grupo de edad
spl = split(metadatos_diversidad$Shannon, metadatos_diversidad$Age_Group)
spl

# Media de Advanced aging
mean(spl$`Advanced aging`)
# Media de Older adults
mean(spl$`Older adults`)
# Media de Young adults
mean(spl$`Young adults`)

# P-valor Young adults vs Older adults
p_value_1 = ks.test(spl$`Young adults`, spl$`Older adults`)$p.value
p_value_1 # Hay diferencias significativas

# P-valor Young adults vs Advanced aging
p_value_2 = ks.test(spl$`Young adults`, spl$`Advanced aging`)$p.value
p_value_2

# P-valor Older adults vs Advanced aging
p_value_3 = ks.test(spl$`Older adults`, spl$`Advanced aging`)$p.value
p_value_3



library(rstatix)

dunn_test(
  metadatos_diversidad,
  Shannon ~ Age_Group,
  p.adjust.method = "BH"
)

# Hay que filtrar otu debido a que tiene NA
otu_matrix <- as(otu, "matrix")
otu_clean <- otu_matrix[complete.cases(otu_matrix), ]

# Verificar dimensiones
dim(otu_clean)

# Verificar que ya no hay NAs
any(is.na(otu_clean))

# Calcular diversidad
diversity_results <- microbiome::diversity(otu_clean)
head(diversity_results)

microbiome::diversity(otu_clean)[order(as.numeric(as.character(microbiome::diversity(otu_clean)$shannon)), decreasing = TRUE),]

# Análisis diversidad beta con Indices de Bray-Curtis
ps.prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
otu_table(physeq)
# Abundancias relativizadas
otu_table(ps.prop)
# Iteración entre las distancias entre las muestras
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "Age_Group", title="Bray NMDS")

# Heatmap
plot_heatmap(physeq, "NMDS", "bray", "Age_Group")


# Creación heatmap
# Filtración de los Top 30 taxones más abundantes
top30_taxas <- names(sort(taxa_sums(physeq), decreasing = TRUE) [1:30])
top30_taxas
phylo_top30 <- prune_taxa(top30_taxas, physeq)
phylo_top30

# Extracción y escala de la matriz (Z-score por fila)
matriz_heatmap <- as(otu_table(phylo_top30), "matrix")
if(!taxa_are_rows(phylo_top30)){matriz_heatmap <- t(matriz_heatmap) }
matriz_escalada <- t(scale(t(matriz_heatmap)))

# Colores para grupos de edad
anotacion_edad <- HeatmapAnnotation(
  Grupo_edad = metadata$Age_Group,
  col = list(Age_group = c(
    "Young adults" = "#8dd3c7",
    "Older adult" = "#fdb462",
    "Advanced aging"= "#fb8072"
  )),
  show_annotation_name = TRUE
)

colores_gradiente <- colorRamp2 (c(-2,0,2), c("#1f78b4", "white", "#e31a1c"))


# Dibujar heatmap
Heatmap(matriz_escalada,
        name = "Z-score",
        top_annotation = anotacion_edad,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE,
        column_title = "Perfil de Microbioma por grupo de edad")

# Heatmap ordenado por edad
Heatmap(matriz_escalada,
        name = "Z-score",
        top_annotation = anotacion_edad,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        column_title = "Perfil de Microbioma por grupo de edad")


