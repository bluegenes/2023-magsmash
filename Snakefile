import os
import pandas as pd

# write snakefile with the following rules:
#1. download fastq files per sample from csv
#2. sketch fastq sample files with sourmash
#3. sketch MAGs with sourmash
#4. run sourmash gather of fastq files against MAGs + GTDB reference
#5. write csv of % matched per sample

outdir = config.get('outdir', 'output')
logs_dir = os.path.join(outdir, 'logs')
HUMAN = "/group/ctbrowngrp/non-microbial-reference/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"

# read in excel sheet to get sample accessions
xls_file = 'inputs/Table_S01_METADATA_metagenomes.xlsx'
sample_info = pd.read_excel(xls_file, sheet_name='(a) Metagenomes', index_col=1, header=1, usecols= "B:S")
SAMPLES = sample_info['NCBI ACCESSION'].dropna().tolist()[:2]
# map to internal names, in case we need it
sample_namemap = sample_info.set_index('NCBI ACCESSION')['Sample Internal ID'].to_dict()

# mag info
magfile = 'inputs/mags.txt'
MAGS = [x.strip() for x in open(magfile, 'r')]

rule all:
    input:
        # expand(os.path.join(outdir, 'raw', '{sample}_{end}.fastq.gz'), sample = SAMPLES, end = [1,2]),
        #expand(os.path.join(outdir, 'fastp', '{sample}_{end}.trim.fq.gz'), sample = SAMPLES, end= [1,2]),
#        expand(os.path.join(outdir, "sketch", "{sample}.{type}.zip"), sample = SAMPLES, type = ['raw', 'nohost']),
#        os.path.join(outdir, "mags", "mags.zip"),
        expand(os.path.join(outdir, "gather", "{sample}.{type}.gather-k{ksize}.csv"), sample = SAMPLES, type= ['raw', 'nohost'], ksize = [21,31,51]),

####################
# download reads #
####################

rule download_reads:
    input: xls_file, 
    output:
        r1=os.path.join(outdir, 'raw', '{sample}_1.fastq.gz'),
        r2=os.path.join(outdir, 'raw', '{sample}_2.fastq.gz'),
    threads: 8,
    resources:
        mem_mb = 3000,
        time = 100,
        partition = 'low2',
    conda: "conf/env/kingfisher.yml",
    log: os.path.join(logs_dir, 'kingfisher_download', '{sample}.log'),
    benchmark: os.path.join(logs_dir, 'kingfisher_download', '{sample}.benchmark'),
    shell:
        """
        kingfisher get -t {threads} -r {wildcards.sample} -m ena-ftp aws-http prefetch 2> {log}
        """

####################
# preprocess reads #
####################

rule fastp:
    input:
        r1 = os.path.join(outdir, 'raw', '{sample}_1.fastq.gz'),
        r2 = os.path.join(outdir, 'raw', '{sample}_2.fastq.gz')
    output:
        r1 = os.path.join(outdir, 'fastp', '{sample}_1.trim.fq.gz'),
        r2 = os.path.join(outdir, 'fastp', '{sample}_2.trim.fq.gz'),
        html = os.path.join(outdir, 'fastp', '{sample}_fastp.html'),
        json = os.path.join(outdir, 'fastp', '{sample}_fastp.json')
    log: os.path.join(logs_dir, 'fastp', '{sample}.log')
    benchmark: os.path.join(logs_dir, 'fastp', '{sample}.benchmark')
    shell:
        """
        fastp -i {input.r1} -I {input.r2}  \
        -o {output.r1} -O {output.r2} \
        -h {output.html} -j {output.json} > {log}
        """

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    input:
        trim1 = os.path.join(outdir, "fastp", "{sample}_1.trim.fq.gz"),
        trim2 = os.path.join(outdir, "fastp", "{sample}_2.trim.fq.gz"),
        human = HUMAN
    output:
        nohost_r1 = os.path.join(outdir, "bbduk", "{sample}_1.nohost.fq.gz"),
        nohost_r2 = os.path.join(outdir, "bbduk", "{sample}_2.nohost.fq.gz"),
        human_r1 = os.path.join(outdir, "bbduk", "{sample}_1.human.fq.gz"),
        human_r2 = os.path.join(outdir, "bbduk", "{sample}_2.human.fq.gz")
    threads: 4
    resources:
        mem_mb = 68000,
    conda: "conf/env/bbmap.yml"
    log: os.path.join(logs_dir, 'remove_host', '{sample}.log')
    benchmark: os.path.join(logs_dir, 'remove_host', '{sample}.benchmark')
    shell:
        """
        bbduk.sh -Xmx64g t={threads} \
        in={input.trim1} in2={input.trim2}\
        out={output.nohost_r1} out2={output.nohost_r2}\
        outm={output.human_r1} outm2={output.human_r2}\
        k=31 ref={input.human}
        """

rule sketch_raw:
    input: 
        r1 = os.path.join(outdir, "raw", "{sample}_1.fastq.gz"),
        r2 = os.path.join(outdir, "raw", "{sample}_2.fastq.gz")
    output:
        sketch = os.path.join(outdir, "sketch", "{sample}.raw.zip")
    threads: 1
    resources:
        mem_mb = 3000,
        time = 100,
        partition = 'low2',
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, 'sourmash', '{sample}.raw.log')
    benchmark: os.path.join(logs_dir, 'sourmash', '{sample}.raw.benchmark')
    shell:
        """
        sourmash sketch -o {output.sketch} {input.r1} {input.r2} 2> {log}
        """


rule sketch_nohost:
    input: 
        nohost_r1 = os.path.join(outdir, "bbduk", "{sample}_1.nohost.fq.gz"),
        nohost_r2 = os.path.join(outdir, "bbduk", "{sample}_2.nohost.fq.gz") 
    output:
        sketch = os.path.join(outdir, "sketch", "{sample}.nohost.zip")
    threads: 1
    resources:
        mem_mb = 3000,
        time = 100,
        partition = 'low2',
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, 'sourmash', '{sample}.log')
    benchmark: os.path.join(logs_dir, 'sourmash', '{sample}.benchmark')
    shell:
        """
        sourmash sketch -p dna,k=21,k=31,k=51,abund \
                        -o {output.sketch} \
                        {input.nohost_r1} {input.nohost_r2} 2> {log}
        """

rule clone_and_create_anvio_env:
    output: 
        anvio_dir = directory("anvio"),
        anvio_dl=".anvio_download",
    log: os.path.join(logs_dir, 'anvio','clone_and_create_env.log')
    benchmark: os.path.join(logs_dir, 'anvio','clone_and_create_env.benchmark')
    shell:
        """
        git clone --recursive https://github.com/merenlab/anvio.git --depth 1
        mamba env create -n anvio-devs -f anvio/.conda/environment.yaml
        touch {output.anvio_dl}
        """

rule anvio_add_path:
    input:
        anvio_dl=".anvio_download",
    output:
        addpath=".anvio_addpath",
    params:
        anvio_dir= "anvio",
    log: os.path.join(logs_dir, 'anvio','add_path.log')
    benchmark: os.path.join(logs_dir, 'anvio','add_path.benchmark')
    conda: "anvio-devs"
    shell:
        """
        pip install -r {params.anvio_dir}/requirements.txt
        set -e  # Abort the script if any command fails
        SHELL_NAME=$(basename "$SHELL")

        # check if 'anvio-devs is in CONDA PREFIX
        if [[ "$CONDA_PREFIX" != *"/anvio-devs"* ]]; then
            echo "anvio-devs is not in CONDA_PREFIX"
            exit 1
        fi

        # Create an activation script for the conda environment (add to path when opening the env)
        mkdir -p "${{CONDA_PREFIX}}/etc/conda/activate.d/"

        SCRIPT_PATH="${{CONDA_PREFIX}}/etc/conda/activate.d/anvio.sh"
        touch "$SCRIPT_PATH"
        chmod +x "$SCRIPT_PATH"
        touch {output.addpath}
        cat <<EOF > "$SCRIPT_PATH"
        export PYTHONPATH=\$PYTHONPATH:$PWD/anvio
        export PATH=\$PATH:$PWD/anvio/bin:$PWD/anvio/sandbox
        echo -e "Updating from anvi'o GitHub (press CTRL+C to cancel) ..."
        cd $PWD/anvio && git pull && cd -
        EOF
        """


rule mags_anviodb_to_fasta:
    input:
        mag_db= 'mags/{mag}/CONTIGS.db',
        addpath=".anvio_addpath",
    output:
        mag_fasta = os.path.join(outdir, "mags", "{mag}.fa")
    threads: 1
    resources:
        mem_mb = 1000,
        time = 100,
        partition = 'low2',
    conda: "anvio-devs"
    log: os.path.join(logs_dir, 'anvi-export-contigs', '{mag}.log')
    benchmark: os.path.join(logs_dir, 'anvi-export-contigs', '{mag}.benchmark')
    shell:
        """
        anvi-export-contigs -c {input.mag_db} -o {output.mag_fasta} 2> {log}
        """

rule write_fromfile_csv:
    input:
        mag_fastas = expand(os.path.join(outdir, "mags", "{mag}.fa"), mag = MAGS)
    output:
        mag_fromfile = os.path.join(outdir, "mags", "mags.fromfile.csv")
    run:
        with open(output.mag_fromfile, 'w') as f:
            f.write('name,genome_filename,protein_filename\n')
            for magfile in mag_fastas:
                mag_name = os.path.basename(magfile).split('.')[0]
                f.write(f'{mag_name},{magfile},\n')


rule sketch_mags:
    input:
        mag_fastas = expand(os.path.join(outdir, "mags", "{mag}.fa"), mag = MAGS),
        mag_fromfile = os.path.join(outdir, "mags", "mags.fromfile.csv")
    output:
        mag_zip = os.path.join(outdir, "mags", "mags.zip")
    threads: 1
    resources:
        mem_mb = 1000,
        time = 100,
        partition = 'low2',
    conda: "conf/env/sourmash.yml"
    log: os.path.join(logs_dir, 'sketch', 'mags.log')
    benchmark: os.path.join(logs_dir, 'sketch', 'mags.benchmark')
    shell:
        """
        sourmash sketch fromfile {input.mag_fromfile} \
                        -p dna,k=21,k=31,k=51,abund \
                        -o {output.mag_zip} 2> {log}
        """


rule gather:
    input:
        sample=os.path.join(outdir, "sketch", "{sample}.{type}.zip"),
        magdb = os.path.join(outdir, "mags", "mags.zip"),
    output:
        gather_csv = os.path.join(outdir, "gather", "{sample}.{type}.gather-k{ksize}.csv"),
        gather_txt = os.path.join(outdir, "gather", "{sample}.{type}.gather-k{ksize}.txt"),
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = 240,
        partition = 'bmm',
    conda:
        "conf/env/sourmash.yml"
    shell:
        """
        sourmash gather --ksize {wildcards.ksize} -o {output.gather_csv} \
                        {input.sample} {input.magdb} > {output.gather_txt} 2> {log}
        """

# rule aggregate_results:
#     input:
#         expand(os.path.join(outdir, "gather", "{sample}.gather-k{{ksize}}.csv"), sample = SAMPLES),
#     output:
#         csv = os.path.join(outdir, "gather", "{basename}.gather-k{{ksize}}.csv")
#     threads: 1
#     resources:
#         mem_mb = 3000,
#         time = 100,
#         partition = 'low2',
#     shell:
#         """
#         python aggregate.py {input.gather} > {output.csv} 2> {log}
#         """