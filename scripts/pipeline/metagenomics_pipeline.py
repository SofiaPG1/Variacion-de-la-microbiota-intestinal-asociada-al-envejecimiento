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
DATA_DIR = os.path.join(BASE_DIR, "Data")
RESULTS_DIR = os.path.join(BASE_DIR, "Results_metasofia")

THREADS = 16
MEGAHIT_MEM = 0.7          # fracción de RAM (baja si hay problemas)
MIN_CONTIG_LEN = 1500      # para MEGAHIT
METABAT_MIN_LEN = 1500     # MetaBAT2 exige mínimo 1500

# Subdirectorios base
DIR_FASTQC    = os.path.join(RESULTS_DIR, "01_FastQC")
DIR_TRIM      = os.path.join(RESULTS_DIR, "fastp")
DIR_ASSEMBLY  = os.path.join(RESULTS_DIR, "03_MEGAHIT")
DIR_BWA       = os.path.join(RESULTS_DIR, "05_BWA")
DIR_SAMTOOLS  = os.path.join(RESULTS_DIR, "06_Samtools")
DIR_METABAT   = os.path.join(RESULTS_DIR, "07_MetaBAT2")
DIR_CHECKM    = os.path.join(RESULTS_DIR, "08_CheckM")
DIR_PROKKA    = os.path.join(RESULTS_DIR, "09_Prokka")
DIR_BAT       = os.path.join(RESULTS_DIR, "10_BAT")

# Paths a ajustar
CAT_CMD = "./CAT_pack/CAT_pack/CAT_pack"
CAT_DB_DIR  = os.path.join(BASE_DIR, "databases", "CAT_nr")
CAT_TAX_DIR = os.path.join(BASE_DIR, "databases", "CAT_taxonomy")

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
    for d in [RESULTS_DIR, DIR_FASTQC, DIR_TRIM, DIR_ASSEMBLY,
              DIR_BWA, DIR_SAMTOOLS, DIR_METABAT, DIR_CHECKM,
              DIR_PROKKA, DIR_BAT]:
        os.makedirs(d, exist_ok=True)

def get_sample_pairs():
    """Detecta pares R1/R2 en Data/"""
    pairs = []
    r1_files = sorted(glob.glob(os.path.join(DATA_DIR, "*_R1*.fastq.gz")) +
                      glob.glob(os.path.join(DATA_DIR, "*_R1*.fq.gz")) +
                      glob.glob(os.path.join(DATA_DIR, "*_1.fastq.gz")) +
                      glob.glob(os.path.join(DATA_DIR, "*_1.fq.gz")))

    for r1 in r1_files:
        base = os.path.basename(r1)
        if "_R1" in base:
            r2 = r1.replace("_R1", "_R2").replace("_r1", "_r2")
            sample = base.split("_R1")[0]
        else:
            r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz").replace("_1.fq.gz", "_2.fq.gz")
            sample = base.replace("_1.fastq.gz", "").replace("_1.fq.gz", "")

        if os.path.exists(r2):
            pairs.append((sample, r1, r2))
    return pairs

def get_trimmed_pairs():
    """Devuelve lista de (sample, r1_trimmed, r2_trimmed) que ya existen"""
    pairs = []
    for r1 in sorted(glob.glob(os.path.join(DIR_TRIM, "*_1_trimmed.fastq.gz"))):
        sample = os.path.basename(r1).replace("_1_trimmed.fastq.gz", "")
        r2 = os.path.join(DIR_TRIM, f"{sample}_2_trimmed.fastq.gz")
        if os.path.exists(r2):
            if os.path.getsize(r1) > 1000 and os.path.getsize(r2) > 1000:
                pairs.append((sample, r1, r2))
            else:
                log(f"ADVERTENCIA: {sample} tiene archivos trimmed muy pequeños → se omite")
    return pairs

# ============================================================
# 1. FastQC
# ============================================================
def run_fastqc(sample=None):
    ensure_dirs()
    pairs = get_sample_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    for samp, r1, r2 in pairs:
        outdir = os.path.join(DIR_FASTQC, samp)
        os.makedirs(outdir, exist_ok=True)
        run_cmd(f"fastqc -t {THREADS} -o {outdir} {r1} {r2}")
        log(f"FastQC terminado: {samp}")
    log("=== FastQC completado ===")

# ============================================================
# 2. fastp
# ============================================================
def run_fastp(sample=None):
    """
    Trimming + QC con fastp (reemplaza Trimmomatic).
    Genera:
      - {sample}_1_trimmed.fastq.gz
      - {sample}_2_trimmed.fastq.gz
      - reportes HTML y JSON
    """
    ensure_dirs()
    pairs = get_sample_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    for samp, r1, r2 in pairs:
        out1 = os.path.join(DIR_TRIM, f"{samp}_1_trimmed.fastq.gz")
        out2 = os.path.join(DIR_TRIM, f"{samp}_2_trimmed.fastq.gz")
        html = os.path.join(DIR_TRIM, f"{samp}_fastp.html")
        json = os.path.join(DIR_TRIM, f"{samp}_fastp.json")

        cmd = (
            f"fastp "
            f"-i {r1} -I {r2} "
            f"-o {out1} -O {out2} "
            f"--thread {THREADS} "
            f"--detect_adapter_for_pe "
            f"--cut_front --cut_tail "
            f"--length_required 36 "
            f"--html {html} "
            f"--json {json} "
            f"--report_title {samp}"
        )
        run_cmd(cmd)
        log(f"fastp terminado: {samp}")
    log("=== fastp completado ===")

# ============================================================
# 3. MEGAHIT - INDIVIDUAL por sample
# ============================================================
def run_megahit(sample=None):
    """
    Ensambla cada sample por separado.
    Si sample=None → procesa todos los que tengan trimmed válido.
    """
    ensure_dirs()
    pairs = get_trimmed_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    if not pairs:
        log("ERROR: No hay samples con archivos trimmed válidos.")
        return

    for samp, r1, r2 in pairs:
        outdir = os.path.join(DIR_ASSEMBLY, samp)

        # Limpiar corrida anterior de este sample
        if os.path.exists(outdir):
            run_cmd(f"rm -rf {outdir}")

        cmd = (
            f"megahit -1 {r1} -2 {r2} "
            f"-o {outdir} "
            f"--num-cpu-threads {THREADS} "
            f"-m {MEGAHIT_MEM} "
            f"--min-contig-len {MIN_CONTIG_LEN} "
            f"--presets meta-sensitive"
        )
        ret = run_cmd(cmd)

        contigs = os.path.join(outdir, "final.contigs.fa")
        if os.path.exists(contigs) and os.path.getsize(contigs) > 0:
            log(f"MEGAHIT OK: {samp} → {contigs}")
        else:
            log(f"ERROR: MEGAHIT no generó contigs para {samp}")
    log("=== MEGAHIT (per-sample) completado ===")


# ============================================================
# 5+6. BWA + Samtools (por sample) - versión final
# ============================================================
def run_bwa_samtools(sample=None):
    ensure_dirs()
    pairs = get_trimmed_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    for samp, r1, r2 in pairs:
        contigs = os.path.join(DIR_ASSEMBLY, samp, "final.contigs.fa")
        if not os.path.exists(contigs):
            log(f"Saltando {samp}: no hay contigs")
            continue

        bwa_dir = os.path.join(DIR_BWA, samp)
        sam_dir = os.path.join(DIR_SAMTOOLS, samp)
        os.makedirs(bwa_dir, exist_ok=True)
        os.makedirs(sam_dir, exist_ok=True)

        index_prefix = os.path.join(bwa_dir, "contigs")
        unsorted_bam = os.path.join(sam_dir, f"{samp}.bam")          # intermedio
        sorted_bam   = os.path.join(sam_dir, f"{samp}.sorted.bam")   # final

        # Si ya existe el BAM final, saltar
        if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 1000000:
            log(f"BAM ya existe, saltando: {samp}")
            # por si acaso quedó el intermedio, lo borramos
            if os.path.exists(unsorted_bam):
                run_cmd(f"rm -f {unsorted_bam}")
            continue

        # Crear índice solo si no existe
        if not os.path.exists(index_prefix + ".bwt"):
            run_cmd(f"bwa index -p {index_prefix} {contigs}")

        # 1. Mapear → BAM intermedio (solo si no existe)
        if not os.path.exists(unsorted_bam) or os.path.getsize(unsorted_bam) < 10000:
            cmd1 = (
                f"bash -c '"
                f"set -o pipefail; "
                f"bwa mem -t {THREADS} {index_prefix} {r1} {r2} | "
                f"samtools view -@ {THREADS} -bS - > {unsorted_bam}"
                f"'"
            )
            ret1 = run_cmd(cmd1)

            if ret1 != 0 or not os.path.exists(unsorted_bam) or os.path.getsize(unsorted_bam) < 10000:
                log(f"ERROR: falló la creación del BAM intermedio para {samp}")
                continue

            log(f"BAM intermedio creado: {unsorted_bam} ({os.path.getsize(unsorted_bam)/1e6:.1f} MB)")
        else:
            log(f"Usando BAM intermedio existente: {unsorted_bam}")

        # 2. Ordenar (sintaxis samtools antiguo con -f)
        cmd2 = f"samtools sort -@ {THREADS} -f {unsorted_bam} {sorted_bam}"
        ret2 = run_cmd(cmd2)

        # 3. Indexar y limpiar SOLO el intermedio
        if ret2 == 0 and os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 1000000:
            run_cmd(f"samtools index -@ {THREADS} {sorted_bam}")
            run_cmd(f"rm -f '{unsorted_bam}'")   # solo borra el .bam intermedio
            log(f"BWA+Samtools OK: {samp} → {sorted_bam}")
        else:
            log(f"ERROR: falló samtools sort para {samp}")
            log(f"  → BAM intermedio dejado en: {unsorted_bam}")

    log("=== BWA + Samtools completado ===")

# ============================================================
# 7. MetaBAT2 (por sample)
# ============================================================
def run_metabat2(sample=None):
    ensure_dirs()
    pairs = get_trimmed_pairs()
    if sample:
        pairs = [p for p in pairs if p[0] == sample]

    for samp, _, _ in pairs:
        contigs = os.path.join(DIR_ASSEMBLY, samp, "final.contigs.fa")
        sorted_bam = os.path.join(DIR_SAMTOOLS, samp, f"{samp}.sorted.bam")

        if not os.path.exists(contigs) or not os.path.exists(sorted_bam):
            log(f"Saltando {samp}: faltan contigs o BAM")
            continue

        out_dir = os.path.join(DIR_METABAT, samp)
        os.makedirs(out_dir, exist_ok=True)

        depth_file = os.path.join(out_dir, "depth.txt")
        bins_prefix = os.path.join(out_dir, "bin")

        run_cmd(f"jgi_summarize_bam_contig_depths --outputDepth {depth_file} {sorted_bam}")

        cmd = (
            f"metabat2 -i {contigs} -a {depth_file} "
            f"-o {bins_prefix} "
            f"-m {METABAT_MIN_LEN} "
            f"-t {THREADS} "
            f"--unbinned"
        )
        run_cmd(cmd)

        n_bins = len(glob.glob(os.path.join(out_dir, "bin.*.fa")))
        log(f"MetaBAT2 {samp}: {n_bins} bins")
    log("=== MetaBAT2 completado ===")

# ============================================================
# 8. CheckM (por sample)
# ============================================================
def run_checkm(sample=None):
    ensure_dirs()
    sample_dirs = sorted([d for d in glob.glob(os.path.join(DIR_METABAT, "*")) if os.path.isdir(d)])
    if sample:
        sample_dirs = [d for d in sample_dirs if os.path.basename(d) == sample]

    for sdir in sample_dirs:
        samp = os.path.basename(sdir)
        bins = glob.glob(os.path.join(sdir, "bin.*.fa"))
        if not bins:
            log(f"Saltando CheckM {samp}: no hay bins")
            continue

        out_dir = os.path.join(DIR_CHECKM, samp)
        os.makedirs(out_dir, exist_ok=True)

        cmd = (
            f"checkm lineage_wf "
            f"-t {THREADS} -x fa --tab_table "
            f"-f {os.path.join(out_dir, 'checkm_summary.tsv')} "
            f"{sdir} {out_dir}"
        )
        run_cmd(cmd)
        log(f"CheckM terminado: {samp}")
    log("=== CheckM completado ===")

# ============================================================
# 9. Prokka (por sample / por bin)
# ============================================================
def run_prokka(sample=None):
    ensure_dirs()
    sample_dirs = sorted([d for d in glob.glob(os.path.join(DIR_METABAT, "*")) if os.path.isdir(d)])
    if sample:
        sample_dirs = [d for d in sample_dirs if os.path.basename(d) == sample]

    for sdir in sample_dirs:
        samp = os.path.basename(sdir)
        bins = sorted(glob.glob(os.path.join(sdir, "bin.*.fa")))
        if not bins:
            continue

        for bin_fa in bins:
            bin_name = os.path.basename(bin_fa).replace(".fa", "")
            outdir = os.path.join(DIR_PROKKA, samp, bin_name)
            if os.path.exists(outdir):
                run_cmd(f"rm -rf {outdir}")

            cmd = (
                f"prokka --outdir {outdir} "
                f"--prefix {bin_name} "
                f"--cpus {THREADS} "
                f"--force "
                f"{bin_fa}"
            )
            run_cmd(cmd)
            log(f"Prokka: {samp}/{bin_name}")
    log("=== Prokka completado ===")

# ============================================================
# 10. BAT (por sample)
# ============================================================
def run_bat(sample=None):
    ensure_dirs()
    if not os.path.isdir(CAT_DB_DIR) or not os.path.isdir(CAT_TAX_DIR):
        log("ADVERTENCIA: Bases de datos CAT no encontradas. Saltando BAT.")
        return

    sample_dirs = sorted([d for d in glob.glob(os.path.join(DIR_METABAT, "*")) if os.path.isdir(d)])
    if sample:
        sample_dirs = [d for d in sample_dirs if os.path.basename(d) == sample]

    for sdir in sample_dirs:
        samp = os.path.basename(sdir)
        bins = glob.glob(os.path.join(sdir, "bin.*.fa"))
        if not bins:
            continue

        out_prefix = os.path.join(DIR_BAT, samp, "BAT")
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

        cmd = (
            f"{CAT_CMD} bins "
            f"-b {sdir} "
            f"-d {CAT_DB_DIR} "
            f"-t {CAT_TAX_DIR} "
            f"-o {out_prefix} "
            f"-s .fa "
            f"--nproc {THREADS}"
        )
        run_cmd(cmd)
        log(f"BAT terminado: {samp}")
    log("=== BAT completado ===")

# ============================================================
# Control
# ============================================================
def setup_check_env():
    tools = ["fastqc", "fastp", "megahit", "bwa", "samtools",
             "jgi_summarize_bam_contig_depths", "metabat2", "checkm", "prokka"]
    log("Verificando herramientas...")
    for t in tools:
        ret = os.system(f"which {t} > /dev/null 2>&1")
        log(f"  {t}: {'OK' if ret == 0 else 'FALTA'}")

def run_all(sample=None):
    log(">>> PIPELINE PER-SAMPLE <<<")
    setup_check_env()
    run_fastqc(sample)
    run_fastp(sample)
    run_megahit(sample)
    run_bwa_samtools(sample)
    run_metabat2(sample)
    run_checkm(sample)
    run_prokka(sample)
    run_bat(sample)
    log(">>> PIPELINE COMPLETADO <<<")

def print_menu():
    print("""
============================================================
 PIPELINE METAGENÓMICO - PER SAMPLE (v3.1 - fastp)
============================================================
 0. setup_check_env()
 1. run_fastqc([sample])
 2. run_fastp([sample])
 3. run_megahit([sample])
 5. run_bwa_samtools([sample])
 6. run_metabat2([sample])
 7. run_checkm([sample])
 8. run_prokka([sample])
 9. run_bat([sample])
-----------------------------------------------------------
 A. run_all()                      - todos los samples
 S. run_all("NOMBRE_SAMPLE")       - solo un sample
 Q. Salir
============================================================
Ejemplos desde Python:
  from pipeline_metasofia_v3_per_sample import *
  run_fastp("ERR7570106")
  run_megahit("SRR341725")
  run_all("SRR341725")
""")

if __name__ == "__main__":
    ensure_dirs()
    print_menu()

    while True:
        choice = input("\nElige opción: ").strip().upper()
        if choice == "0":
            setup_check_env()
        elif choice == "1":
            run_fastqc()
        elif choice == "2":
            run_fastp()
        elif choice == "3":
            run_megahit()
        elif choice == "5":
            run_bwa_samtools()
        elif choice == "6":
            run_metabat2()
        elif choice == "7":
            run_checkm()
        elif choice == "8":
            run_prokka()
        elif choice == "9":
            run_bat()
        elif choice == "A":
            run_all()
        elif choice == "Q":
            break
        else:
            print("Opción no válida. Usa 0-9, A o Q.")
            print("Para un sample específico ejecuta desde Python interactivo.")