## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS pipeline
## Authors: Floris Barthel, Samir Amin, Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import os
import pandas as pd
import itertools

## Import manifest processing functions
from python.glassfunc import dbconfig, locate
from python.PostgreSQLManifestHandler import PostgreSQLManifestHandler
from python.JSONManifestHandler import JSONManifestHandler

## Connect to database
## dbconf = dbconfig("/home/barthf/.odbc.ini", "VerhaakDB")
dbconf = dbconfig(config["db"]["configfile"], config["db"]["configsection"])

#print("Cohort set to ", str(config["cohort"]))
by_cohort = None
if len(str(config["cohort"])) > 0:
    by_cohort = str(config["cohort"]).zfill(2)

## Instantiate manifest
manifest = PostgreSQLManifestHandler(host = dbconf["servername"], port = dbconf["port"], user = dbconf["username"], password = dbconf["password"], database = dbconf["database"],
    source_file_basepath = config["data"]["source_path"], aligned_file_basepath = config["data"]["realn_path"], from_source = config["from_source"], by_cohort = by_cohort)
print(manifest)

## Set working directory based on configuration file
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE = config["gdc_token"]

## Cluster metadata (memory, CPU, etc)
CLUSTER_META = json.load(open(config["cluster_json"]))

## List of scatterlist items to iterate over
## Each Mutect2 run spawns 50 jobs based on this scatterlist

WGS_SCATTERLIST = ["temp_{num}_of_50".format(num=str(j+1).zfill(4)) for j in range(50)]

## Load modules

## We do not want the additional DAG processing if not from source
#if(config["from_source"]):

#DNA modules
#include: "snakemake/download.smk"
#include: "snakemake/align.smk"
#include: "snakemake/haplotype-map.smk"
#include: "snakemake/fingerprinting.smk"
#include: "snakemake/telseq.smk"
#include: "snakemake/mutect2.smk"
#include: "snakemake/mutect2-post.smk"
#include: "snakemake/varscan2.smk"
#include: "snakemake/cnvnator.smk"
#include: "snakemake/lumpy.smk"
#include: "snakemake/delly.smk"
#include: "snakemake/manta.smk"
include: "snakemake/cnv.smk"
#include: "snakemake/sequenza.smk"
#include: "snakemake/optitype.smk"
#include: "snakemake/pvacseq.smk"
#include: "snakemake/lohhla.smk"
#include: "snakemake/cnv-post.smk"		(Deprecated)
#include: "snakemake/titan.smk"
#include: "snakemake/pyclone.smk"

#RNA modules
#include: "snakemake/kallisto.smk"
#include: "snakemake/kallisto_nc.smk"
#include: "snakemake/rna_fingerprinting.smk"
#include: "snakemake/mixcr.smk"
#include: "snakemake/prada.smk"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Upload coverage to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#This is currently done manually as new cohorts are added, as this rule does not support cohort-specific uploads

# rule cov2db:
#     input:
#         metrics = lambda wildcards: expand("results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt", aliquot_barcode = manifest.getSelectedAliquots())
#     output:
#         tsv = "results/align/wgsmetrics.merged.tsv"
#     params:
#         mem = CLUSTER_META["cov2db"]["mem"]
#     threads:
#         CLUSTER_META["cov2db"]["ppn"]
#     #conda:
#     #    "../envs/r.yaml"
#     log:
#         "logs/align/cov2db/cov2db.log"
#     benchmark:
#         "benchmarks/align/cov2db/cov2db.txt"
#     message:
#         "Merge coverage file and convert to TSV for easy database upload"
#     script:
#         "R/snakemake/cov2db.R"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Haplotype map creation rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#rule build_haplotype_map:
#    input:
#        "data/ref/fingerprint.filtered.map"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Alignment rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule align:
    input:
        expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/align/wgsmetrics/{aliquot_barcode}.WgsMetrics.txt", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/align/validatebam/{aliquot_barcode}.ValidateSamFile.txt", aliquot_barcode = manifest.getSelectedAliquots()),
        lambda wildcards: ["results/align/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample = aliquot_barcode, rg = readgroup)
          for aliquot_barcode, readgroups in manifest.getSelectedReadgroupsByAliquot().items()
          for readgroup in readgroups]

rule gencode:
    input:
        expand("results/align/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv", aliquot_barcode = manifest.getAliquotsByBatch('GLSS-H2-WXS'))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download only rule
## Run snakemake with 'snakemake download_only' to activate
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#rule download_only:
#   input: expand("{file}", file = ALIQUOT_TO_BAM_PATH.values())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## QC rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule qc:
    input: 
        "results/align/multiqc/multiqc_report.html"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (Mutect2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mutect2:
    input:
        expand("results/mutect2/m2filter/{case_barcode}.filtered.vcf.gz", case_barcode = manifest.getSelectedCases()),
        expand("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz", pair_barcode = manifest.getSelectedPairs())

rule ssmutect2:
    input:
        expand("results/mutect2/ssm2filter/{pair_barcode}.filtered.vcf.gz", pair_barcode = manifest.getSelectedPairs())

rule m2db:
    input:
        "results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv",
        "results/mutect2/consensusvcf/consensus.normalized.sorted.vep.maf",
        expand("results/mutect2/geno2db/{case_barcode}.info.tsv", case_barcode = manifest.getSelectedCases()),
        expand("results/mutect2/geno2db/{case_barcode}.geno.tsv", case_barcode = manifest.getSelectedCases())

rule sequenza:
    input:
        expand("results/sequenza/mergeseqz/{pair_barcode}.small.seqz.gz", pair_barcode = manifest.getSelectedPairs()),
        expand("results/sequenza/seqzR/{pair_barcode}/{pair_barcode}_cellularity.ploidy.txt", pair_barcode = manifest.getSelectedPairs())

rule titancna:
    input:
        expand("results/cnv/titanfinal/seg/{pair_barcode}.seg.txt", pair_barcode = manifest.getSelectedPairs())

# rule mutect2post:
# 	input:
# 		expand("results/mutect2/m2post/{pair_barcode}.normalized.sorted.vcf.gz", pair_barcode = manifest.getSelectedPairs())

# rule genotypefreebayes:
#     input:
#         expand("results/mutect2/freebayes/{aliquot_barcode}.normalized.sorted.vcf.gz", aliquot_barcode = manifest.getSelectedAliquots())

# rule massfreebayes:
#     input:
#         expand("results/mutect2/freebayes/batch{batch}/{aliquot_barcode}.normalized.sorted.vcf.gz", batch = [str(i) for i in range(2,6)], aliquot_barcode = manifest.getSelectedAliquots())

# rule genotypevcf2vcf:
#     input:
#         expand("results/mutect2/genotypes/{aliquot_barcode}.normalized.sorted.vcf.gz", aliquot_barcode = manifest.getSelectedAliquots())

# rule finalfreebayes:
#     input:
#         expand("results/mutect2/batches2db/{aliquot_barcode}.normalized.sorted.tsv", aliquot_barcode = manifest.getSelectedAliquots())

# rule preparem2pon:
#     input:
#         expand("results/mutect2/mergepon/{aliquot_barcode}.pon.vcf", aliquot_barcode = manifest.getPONAliquots())

# rule genodb:
#     input:
#         expand("results/mutect2/geno2db/{aliquot_barcode}.normalized.sorted.tsv", aliquot_barcode = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## PON rule (Mutect2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# rule mutect2pon:
#     input:
#     	"results/mutect2/pon/pon.vcf"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## HLAtyping rule (OptiType)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule call_hla:
	input:
		expand("results/optitype/HLA_calls/{aliquot_barcode}/{aliquot_barcode}_result.tsv", aliquot_barcode = manifest.getSelectedAliquots()),
		expand("results/optitype/HLA_calls/{aliquot_barcode}/{aliquot_barcode}_extended_result.tsv", aliquot_barcode = manifest.getSelectedAliquots())
	
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Neoantigen rule (pVACseq)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule call_neoag:
	input:
		expand("results/pvacseq/neoantigens/{case_barcode}/MHC_Class_I/{case_barcode}.final.tsv", case_barcode = manifest.getSelectedCases())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (VarScan2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule varscan2:
    input:
        expand("results/varscan2/fpfilter/{pair_barcode}.{type}.Somatic.hc.final.vcf", pair_barcode = manifest.getSelectedPairs(), type = ["snp", "indel"])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling pipeline
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnv:
    input:
        expand("results/cnv/plots/{aliquot_barcode}.pdf", aliquot_barcode = manifest.getSelectedAliquots()),
        expand("results/cnv/callsegments/{aliquot_barcode}.called.seg", aliquot_barcode = manifest.getSelectedAliquots()),
        "results/cnv/callsegments.merged.tsv"		#May need to remove this if it breaks the pipeline

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call loss of heterozygosity in HLA regions (LOHHLA)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lohhla:
    input:
        expand("results/lohhla/final/{pair_barcode}/{minCoverageFilter}/{pair_barcode}.{minCoverageFilter}.DNA.HLAlossPrediction_CI.xls", pair_barcode = manifest.getSelectedPairs(), minCoverageFilter=[5,10,20,30])
        #expand("results/lohhla/final/{pair_barcode}/{pair_barcode}.DNA.HLAlossPrediction_CI.txt", pair_barcode = manifest.getSelectedPairs())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate TL using telseq
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule telseq:
    input:
        expand("results/telseq/{aliquot_barcode}.telseq.txt", aliquot_barcode = manifest.getSelectedAliquots())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run PyClone
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule pyclone:
    input:
        #lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/plots/loci/{plot_type}.pdf", pyclone_short_name = manifest.getPyCloneCases(), plot_type = ['density','parallel_coordinates','scatter','similarity_matrix','vaf_parallel_coordinates','vaf_scatter']),
        #lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/plots/clusters/{plot_type}.pdf", pyclone_short_name = manifest.getPyCloneCases(), plot_type = ['density','parallel_coordinates','scatter']),
        lambda wildcards: expand("results/pyclone/run/{pyclone_short_name}/tables/{table_type}.tsv", pyclone_short_name = manifest.getPyCloneCases(), table_type = ['cluster','loci'])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Fingerprinting pipeline
## Check sample and case fingerprints
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprint:
   input:
       expand("results/fingerprinting/sample/{aliquot_barcode}.crosscheck_metrics", aliquot_barcode = manifest.getSelectedAliquots()),
       expand("results/fingerprinting/case/{case_barcode}.crosscheck_metrics", case_barcode = manifest.getSelectedCases()) ## ,
##       "results/fingerprinting/GLASS_HF.crosscheck_metrics"

# Fingerprint DNA and RNA samples together
rule full_fingerprint:
    input:
       lambda wildcards: ["results/rnafingerprint/star/{sample}/{sample}.{rg}.ReadsPerGene.out.tab".format(sample = aliquot_barcode, rg = readgroup)
       	for aliquot_barcode, readgroups in manifest.getSelectedReadgroupsByAliquot(analyte='R').items()
       	for readgroup in readgroups],
       #expand("results/rnafingerprint/sample/{aliquot_barcode}.crosscheck_metrics", aliquot_barcode = manifest.getSelectedAliquots(analyte='R')),
       expand("results/rnafingerprint/case/{case_barcode}.crosscheck_metrics", case_barcode = manifest.getSelectedCases()),
       expand("results/rnafingerprint/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv", aliquot_barcode = manifest.getSelectedAliquots(analyte='R'))

#         expand("results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getSelectedAliquots(analyte = 'R')),
#         expand("results/rnafingerprint/validatebam/{aliquot_barcode}.ValidateSamFile.txt", aliquot_barcode = manifest.getSelectedAliquots(analyte = 'R')),
#         lambda wildcards: ["results/rnafingerprint/fastqc/{sample}/{sample}.{rg}.unaligned_fastqc.html".format(sample = aliquot_barcode, rg = readgroup)
#           for aliquot_barcode, readgroups in manifest.getSelectedReadgroupsByAliquot().items()
#           for readgroup in readgroups]
       
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run kallisto pipeline
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule quant_tpm:
    input:
       "results/kallisto/kallisto/final/transcript_tpms_all_samples.tsv",
       "results/kallisto/kallisto/final/transcript_tpm_matrix_all_samples.tsv",
       "results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv",
       "results/kallisto/kallisto/final/p_result_gene_tpm_matrix_all_samples.gct.txt",
       "results/kallisto/pizzly/final/fusions_all_samples.tsv"

rule quant_nctpm:
    input:
       "results/kallisto/noncoding/final/nc_transcript_tpms_all_samples.tsv",
       "results/kallisto/noncoding/final/nc_transcript_tpm_matrix_all_samples.tsv",
       "results/kallisto/noncoding/final/noncoding_tpm_matrix_all_samples.tsv"
       
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run PRADA for fusions and EGFRvIII variants
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule prada:
    input:
        expand("results/prada/{aliquot_barcode}/fusion/prada.fus.summary.txt", aliquot_barcode = manifest.getSelectedAliquots(analyte = 'R')),
        expand("results/prada/{aliquot_barcode}/fusion/prada.fus.summary.taf.txt", aliquot_barcode = manifest.getSelectedAliquots(analyte = 'R')),
        expand("results/prada/{aliquot_barcode}/guess-if/{gene}/{gene}.GUESS-IF.summary.txt", aliquot_barcode = manifest.getSelectedAliquots(analyte = 'R'), gene=['EGFR'])
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run MiXCR
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mixcr:
    input:
        "results/mixcr/final/merged_clones_tcr.txt",
        "results/mixcr/final/merged_clones_bcr.txt"
                
          
## END ##
