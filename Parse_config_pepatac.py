from pathlib import Path
import snakemake.io
from collections import defaultdict
import itertools
import os
import subprocess

class ParsedConfig:

    def __init__(self, config):
        self.ROOT = Path(config["ROOT"]).resolve()
        self.genome_folder = Path(config['PEPATAC_params']['genome_folder']).resolve()
        self.genome = config['PEPATAC_params']["genome_version"]
        self.DATASET = config["DATASET"]
        self.fastq_suffix1 = config[self.DATASET]["fasq_suffix1"]
        self.fastq_suffix2 = config[self.DATASET]["fasq_suffix2"]

        os.makedirs(os.path.join(config["ROOT"], config["DATASET"]), exist_ok=True)
        self.OUTDIR = os.path.join(config["ROOT"], config["DATASET"])

        print("self.DATASET:", repr(self.DATASET))
        print("config keys:", list(config.keys()))
        print("config[ self.DATASET ]:", config.get(self.DATASET))
        # Add samples names to each cohort
        if config[self.DATASET]["is_geo"] == True:
            samples_IDs = []
            with open(os.path.join(config["snakemake_path"], config[self.DATASET]["SRA_IDs_list"])) as f:
                samples_IDs.extend(f.read().strip().split('\n'))
            self.Sample_IDs = samples_IDs
            self.fastq_path = Path(config[self.DATASET]["sra_folder"]).resolve()
        else:
            print("No SRA IDs for this dataset.")
            samples_IDs = []
            with open(os.path.join(config["snakemake_path"], config[self.DATASET]["Sample_IDs_list"])) as f:
                samples_IDs.extend(f.read().strip().split('\n'))
            self.Sample_IDs = samples_IDs
            self.fastq_path = Path(config[self.DATASET]["fastq_files"]).resolve()

        # Define variables
        path_adapters = config['PEPATAC_params']["adapters"]
        original_config = config['PEPATAC_params']["yaml_file"]
        updated_config = config["ROOT"] + self.DATASET + "_pepatac.yaml"
        command = f"sed 's|adapters: null|adapters: {path_adapters}|' {original_config} > {updated_config}"
        subprocess.run(command, shell=True, check=True)
        self.pepatac_config = updated_config

    def generate_targets_pepatac_inputs(self, config):
        targets = []
        if self.genome == 'mm10':
            # targets.append(f"{self.genome_folder}/alias/mm10/fasta/default/mm10.fa")
            # targets.append(f"{self.genome_folder}/alias/mm10/fasta/default/mm10.fai")
            # targets.append(f"{self.genome_folder}/alias/mm10/fasta/default/mm10.chrom.sizes")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.3.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.4.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.rev.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/bowtie2_index/default/mm10.rev.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/mm10/refgene_anno/default/mm10_TSS.bed")
            # targets.append(f"{self.genome_folder}/alias/mm10/refgene_anno/default/mm10_exons.bed")
            # targets.append(f"{self.genome_folder}/alias/mm10/refgene_anno/default/mm10_introns.bed")
            # targets.append(f"{self.genome_folder}/alias/mm10/refgene_anno/default/mm10_pre-mRNA.bed")
            # targets.append(f"{self.genome_folder}/alias/mm10/refgene_anno/default/mm10_refGene.txt.gz")
            # targets.append(f"{self.genome_folder}/alias/mm10/ensembl_gtf/default/mm10.gtf.gz")
            # targets.append(f"{self.genome_folder}/alias/mm10/ensembl_gtf/default/mm10_ensembl_TSS.bed")
            # targets.append(f"{self.genome_folder}/alias/mm10/ensembl_gtf/default/mm10_ensembl_gene_body.bed")
            targets.append(f"{self.genome_folder}/alias/mm10/ensembl_rb/default/mm10.gff.gz")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/fasta/default/mouse_chrM2x.chrom.sizes")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/fasta/default/mouse_chrM2x.fa")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/fasta/default/mouse_chrM2x.fa.fai")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.3.bt2")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.4.bt2")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.rev.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x.rev.2.bt2")
        if self.genome == 'hg38':
            # targets.append(f"{self.genome_folder}/alias/hg38/fasta/default/hg38.fa")
            # targets.append(f"{self.genome_folder}/alias/hg38/fasta/default/hg38.fai")
            # targets.append(f"{self.genome_folder}/alias/hg38/fasta/default/hg38.chrom.sizes")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.3.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.4.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.rev.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/bowtie2_index/default/hg38.rev.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/hg38/refgene_anno/default/hg38_TSS.bed")
            # targets.append(f"{self.genome_folder}/alias/hg38/refgene_anno/default/hg38_exons.bed")
            # targets.append(f"{self.genome_folder}/alias/hg38/refgene_anno/default/hg38_introns.bed")
            # targets.append(f"{self.genome_folder}/alias/hg38/refgene_anno/default/hg38_pre-mRNA.bed")
            # targets.append(f"{self.genome_folder}/alias/hg38/refgene_anno/default/hg38_refGene.txt.gz")
            # targets.append(f"{self.genome_folder}/alias/hg38/ensembl_gtf/default/hg38.gtf.gz")
            # targets.append(f"{self.genome_folder}/alias/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed")
            # targets.append(f"{self.genome_folder}/alias/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed")
            targets.append(f"{self.genome_folder}/alias/hg38/ensembl_rb/default/hg38.gff.gz")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/fasta/default/human_repeats.chrom.sizes")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/fasta/default/human_repeats.fa")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/fasta/default/human_repeats.fa.fai")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.3.bt2")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.4.bt2")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.rev.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/human_repeats/bowtie2_index/default/human_repeats.rev.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/fasta/default/rCRSd.chrom.sizes")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/fasta/default/rCRSd.fa")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/fasta/default/rCRSd.fa.fai")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.2.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.3.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.4.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.rev.1.bt2")
            # targets.append(f"{self.genome_folder}/alias/rCRSd/bowtie2_index/default/rCRSd.rev.2.bt2")
        return targets

    def get_sra_output_files(self, config):
        targets = []
        for sample in self.Sample_IDs:
            targets.append(f"{self.fastq_path}/{sample}_{self.fastq_suffix1}"),
            targets.append(f"{self.fastq_path}/{sample}_{self.fastq_suffix2}")
        return targets

    def get_pepatac_output_files(self, config):
        targets = []
        for sample in self.Sample_IDs:
            targets.append(f"{self.OUTDIR}/{sample}/aligned_{self.genome}/{sample}_sort_dedup.bam"),
            targets.append(f"{self.OUTDIR}/{sample}/aligned_{self.genome}/{sample}_sort_dedup_shifted.bam"),
            targets.append(f"{self.OUTDIR}/{sample}/peak_calling_bampe_{self.genome}/{sample}_peaks_fixedWidth.narrowPeak"),
            targets.append(f"{self.OUTDIR}/{sample}/peak_calling_{self.genome}/{sample}_peaks_rmBlacklist.narrowPeak"),
            targets.append(f"{self.OUTDIR}/fragments/{sample}_fragments.gz"),
            targets.append(f"{self.OUTDIR}/fragments2/{sample}_fragments.gz")
        return targets
