#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import time

# ============================================================
# CONFIGURACIÓN
# ============================================================
BASE_DIR = os.getcwd()

# Entrada: trimmed de fastp (ajusta si están en otra ruta)
DIR_TRIM = os.path.join(BASE_DIR, "Results_metasofia", "fastp")

# Salida
RESULTS_DIR = os.path.join(BASE_DIR, "Results_kraken")
DIR_KRAKEN  = os.path.join(RESULTS_DIR, "01_Kraken2")
DIR_BRACKEN = os.path.join(RESULTS_DIR, "02_Bracken")

# Base de datos Kraken2
KRAKEN_DB = os.path.abspath(os.path.join(BASE_DIR, "../kraken/k2_standard_20260226"))

THREADS = 16
READ_LEN = 150          # longitud promedio de reads (para Bracken)
BRACKEN_LEVEL = "S"     # S = species, G = genus, etc.
BRACKEN_THRESHOLD = 10  # mínimo de reads para reportar

# ============================================================
# UTILIDADES
# ============================================================
def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}")
    sys.stdout.flush()

def run_cmd(cmd, check=True):
    log(f"CMD: {cmd}")
    ret = os.system(cmd)
    if check and ret != 0:
        log(f"ERROR: comando falló con código {ret}")
    return ret

def ensure_dirs():
    for d in [RESULTS_DIR, DIR_KRAKEN, DIR_BRACKEN]:
        os.makedirs(d, exist_ok=True)

def get_trimmed_pairs():
    """Detecta pares trimmed en DIR_TRIM"""
    pairs = []
    for r1 in sorted(glob.glob(os.path.join(DIR_TRIM, "*_1_trimmed.fastq.gz"))):
        sample = os.path.basename(r1).replace("_1_trimmed.fastq.gz", "")
        r2 = os.path.join(DIR_TRIM, f"{sample}_2_trimmed.fastq.gz")
        if os.path.exists(r2) and os.path.getsize(r1) > 1000 and os.path.getsize(r2) > 1000:
            pairs.append((sample, r1, r2))
        else:
            log(f"ADVERTENCIA: {sample} sin par válido → se omite")
    return pairs

# ============================================================
# 1. Kraken2
# ============================================================
def run_kraken2(sample=None):
    ensure_dirs()

    if not os.path.isdir(KRAKEN_DB):
        log(f"ERROR: No se encuentra la base de datos: {KRAKEN_DB}")
        return

    pairs = get_trimmed_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    if not pairs:
        log("ERROR: No hay samples con trimmed válidos.")
        return

    for samp, r1, r2 in pairs:
        out_dir = os.path.join(DIR_KRAKEN, samp)
        os.makedirs(out_dir, exist_ok=True)

        report  = os.path.join(out_dir, f"{samp}.kraken2.report")
        output  = os.path.join(out_dir, f"{samp}.kraken2.out")
        classified = os.path.join(out_dir, f"{samp}_classified")
        unclassified = os.path.join(out_dir, f"{samp}_unclassified")

        # Si ya existe el report, saltar
        if os.path.exists(report) and os.path.getsize(report) > 100:
            log(f"Kraken2 ya existe, saltando: {samp}")
            continue

        cmd = (
            f"kraken2 "
            f"--db {KRAKEN_DB} "
            f"--threads {THREADS} "
            f"--paired "
            f"--use-names "
            f"--report {report} "
            f"--output {output} "
            f"--classified-out {classified}#.fastq "
            f"--unclassified-out {unclassified}#.fastq "
            f"{r1} {r2}"
        )
        ret = run_cmd(cmd)

        if ret == 0 and os.path.exists(report):
            log(f"Kraken2 OK: {samp}")
        else:
            log(f"ERROR: Kraken2 falló para {samp}")

    log("=== Kraken2 completado ===")

# ============================================================
# 2. Bracken
# ============================================================
def run_bracken(sample=None):
    ensure_dirs()

    # Buscar reports de Kraken2
    report_files = sorted(glob.glob(os.path.join(DIR_KRAKEN, "*", "*.kraken2.report")))
    if sample:
        report_files = [r for r in report_files if sample in r]

    if not report_files:
        log("ERROR: No se encontraron reports de Kraken2. Corre primero run_kraken2()")
        return

    for report in report_files:
        samp = os.path.basename(os.path.dirname(report))
        out_dir = os.path.join(DIR_BRACKEN, samp)
        os.makedirs(out_dir, exist_ok=True)

        bracken_out = os.path.join(out_dir, f"{samp}.bracken")
        bracken_report = os.path.join(out_dir, f"{samp}.bracken.report")

        if os.path.exists(bracken_out) and os.path.getsize(bracken_out) > 100:
            log(f"Bracken ya existe, saltando: {samp}")
            continue

        cmd = (
            f"bracken "
            f"-d {KRAKEN_DB} "
            f"-i {report} "
            f"-o {bracken_out} "
            f"-w {bracken_report} "
            f"-r {READ_LEN} "
            f"-l {BRACKEN_LEVEL} "
            f"-t {BRACKEN_THRESHOLD}"
        )
        ret = run_cmd(cmd)

        if ret == 0 and os.path.exists(bracken_out):
            log(f"Bracken OK: {samp}")
        else:
            log(f"ERROR: Bracken falló para {samp}")

    log("=== Bracken completado ===")

# ============================================================
# Control
# ============================================================
def setup_check_env():
    tools = ["kraken2", "bracken"]
    log("Verificando herramientas...")
    for t in tools:
        ret = os.system(f"which {t} > /dev/null 2>&1")
        log(f"  {t}: {'OK' if ret == 0 else 'FALTA'}")

    if os.path.isdir(KRAKEN_DB):
        log(f"  Base de datos: OK → {KRAKEN_DB}")
    else:
        log(f"  Base de datos: FALTA → {KRAKEN_DB}")

def run_all(sample=None):
    log(">>> KRAKEN2 + BRACKEN <<<")
    setup_check_env()
    run_kraken2(sample)
    run_bracken(sample)
    log(">>> COMPLETADO <<<")

def print_menu():
    print("""
============================================================
 KRAKEN2 + BRACKEN (per-sample)
============================================================
 0. setup_check_env()
 1. run_kraken2([sample])
 2. run_bracken([sample])
-----------------------------------------------------------
 A. run_all()                 - todos los samples
 S. run_all("NOMBRE_SAMPLE")  - solo un sample
 Q. Salir
============================================================
Ejemplos:
  run_kraken2("SRR341694")
  run_bracken("SRR341694")
  run_all("ERR7570166")
""")

if __name__ == "__main__":
    ensure_dirs()
    print_menu()

    while True:
        choice = input("\nElige opción: ").strip().upper()
        if choice == "0":
            setup_check_env()
        elif choice == "1":
            run_kraken2()
        elif choice == "2":
            run_bracken()
        elif choice == "A":
            run_all()
        elif choice == "Q":
            break
        else:
            print("Opción no válida. Usa 0, 1, 2, A o Q.")
            print("Para un sample específico usa Python interactivo.")