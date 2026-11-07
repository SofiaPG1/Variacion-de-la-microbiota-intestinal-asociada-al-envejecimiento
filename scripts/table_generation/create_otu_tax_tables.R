# Carga de paquetes necesarios
library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)
library(phyloseq)
library(rentrez)
library(xml2)

# Lectura del archivo
datos <- read_tsv("tabla_abundancia_combinada.tsv")
# Comprobación primeras filas
head(datos)

# Extracción columnas conteos (num)
datos_num <- datos %>%
  select(name, taxonomy_id, taxonomy_lvl, ends_with("_num")) %>%
  # Renombra las columnas para quitar el sufijo "_num"
  rename_with(~gsub("_num$", "", .x), ends_with("num"))

# Extracción columnas fracciones (frac)
datos_frac <- datos %>%
  select(name, taxonomy_id, taxonomy_lvl, ends_with("_frac")) %>%
  # Renombra las columnas para quitar el sufijo "_frac"
  rename_with(~gsub("_frac$", "", .x), ends_with("_frac"))

# Visualización de los resultados
view(datos_num)
view(datos_frac)

# Conversión formato ancho a largo
# Para los conteos (num)
datos_num_long <- datos_num %>%
  pivot_longer(
    cols = -c(name, taxonomy_id, taxonomy_lvl),
    names_to = "Muestra",
    values_to = "Abundancia"
  ) %>%
  filter(Abundancia > 0)
# Comprobación
view(datos_num_long)

# Para las fracciones (frac)
# Para los conteos (num)
datos_frac_long <- datos_frac %>%
  pivot_longer(
    cols = -c(name, taxonomy_id, taxonomy_lvl),
    names_to = "Muestra",
    values_to = "Fraccion"
  ) %>%
  filter(Fraccion > 0)
# Comprobación
view(datos_frac_long)

# Verificación de los nombres de las columnas
names(datos_num)

# Creación matriz de conteos (taxones en las filas y las muestras en columnas)
matriz_conteos <- datos_num %>%
  column_to_rownames("name") %>%
  select(-taxonomy_id, -taxonomy_lvl) %>%
  rename_with(~sub("\\.bracken$", "", .x)) %>%
  as.matrix()

# Creación de la matriz de fracciones
matriz_fracciones <- datos_frac %>%
  column_to_rownames("name") %>%
  select(-taxonomy_id, -taxonomy_lvl) %>%
  rename_with(~sub("\\.bracken$", "", .x)) %>%
  as.matrix()

# Visualización dimensiones de las matrices
dim(matriz_conteos)
dim(matriz_fracciones)

# Tabla resumen de abundancia total por especie
resumen_taxones <- datos_num_long %>%
  group_by(name) %>%
  summarise(
    Abundancia_total = sum(Abundancia),
    Presencia = n(),
    Media = mean(Abundancia),
    Máximo = max(Abundancia)
  ) %>%
  arrange(desc(Abundancia_total))

# Visualización de los 10 taxones más abundantes
head(resumen_taxones, 10)

# Creación matriz de abundancia (OTU table)

otu_table <- matriz_conteos


# Tabla de taxomía básica
tax_table_basica <- datos %>%
  select(name, taxonomy_id, taxonomy_lvl) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Creación función para obtener taxonomía
# Adaptada de scripts comunitarios de bioinformática en R (GitHub/ Stack Overflow)
# Su función es extraer los 7 rangos taxonómicos principales por medio de la API de NCBI



# Tax table con la estructura jerarquica
tax_matriz <- datos %>%
  mutate(
    Kingdom = "Bacteria",
    Phylum = NA_character_,
    Class = NA_character_,
    Order = NA_character_,
    Family = NA_character_,
    Genus = gsub(" .*$", "", name),
    Species = name,
    NCBI_ID = as.character(taxonomy_id)
  ) %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, NCBI_ID) %>%
  column_to_rownames("Species") %>%
  as.matrix()


# Convertir las tablas a data.frame
otu_table_df <- as.data.frame(otu_table)
tax_table_df <- as.data.frame(tax_matriz)

# Guarda la tabla otu_table en formato .tsv
write.table(otu_table_df, "otu_table.tsv", sep = "\t", quote = FALSE)
write.table(tax_table_df, "tax_table.tsv", sep = "\t", quote = FALSE)
