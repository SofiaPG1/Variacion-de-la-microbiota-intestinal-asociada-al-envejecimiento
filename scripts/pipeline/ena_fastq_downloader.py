#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Importación de las librerias necesarias
import os
import csv

# Función encargada de descargar los archivos FASTQ
def download():

    # Abre el archivo csv que contiene las muestras y URLs de descarga
    with open("./samples.csv", "r") as links:

        # Lectura del csv utilizando los nombres de las columnas como claves
        readlinks= csv.DictReader(links)

        # Recorre cada fila del archivo csv
        for row in readlinks:

            # Extrae los enlaces de la muestra
            link_R1 = row["Enlace_R1"].strip()
            link_R2 = row["Enlace_R2"].strip()

            # Extrae el identificardor de la url
            identificador = os.path.basename(link_R1).split("_")[0]
            identificador_R2 = os.path.basename(link_R2).split("_")[0]

            # Comprobación identificadores iguales tanto en foward como reverse
            if identificador != identificador_R2:
                print("Los archivos R1 y R2 no pertenecen al mismo")
                continue

            # Comprueba que los archivos tengan información
            if not identificador or not link_R1 or not link_R2:
                print(f"Datos incompletos para {identificador}")
                continue

            # Creación carpeta para cada muestra
            carpeta = os.path.join("rawdata", identificador)
            os.makedirs(carpeta, exist_ok=True)

            # Muestra el nombre de la muestra que se está descargando
            print(f"\nDescargando muestra: {identificador}")

            # Descarga archivo FASTQ forward (R1)
            instruccion_R1 = f'wget -c -P "{carpeta}" "{link_R1}"'
            resultado_R1 = os.system(instruccion_R1)
        
            # Comprobación del resultado de la descarga R1
            if resultado_R1 == 0:
                print(f"R1 descargado correctamente: {identificador}")
            else:
                print(f"Error descargando R1: {identificador}")

            # Descarga archivo FASTQ reverse (R2)
            instruccion_R2 = f'wget -c -P "{carpeta}" "{link_R2}"'
            resultado_R2 = os.system(instruccion_R2)

             # Comprobación del resultado de la descarga R1
            if resultado_R2 == 0:
                print(f"R2 descargado correctamente: {identificador}")
            else:
                print(f"Error descargando R2: {identificador}")
      
# Ejecución de la función principal
if __name__ == "__main__":
        download()
        