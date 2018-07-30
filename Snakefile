## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Snakefile for GLASS-WG pipeline
## Authors: Floris Barthel, Samir Amin
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

import pandas as pd
import itertools
import os

## Touch function taken from stackoverflow
## Link: https://stackoverflow.com/questions/1158076/implement-touch-using-python
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)

## Turn an unnamed list of dicts into a nammed list of dicts
## Taken from stackoverflow
## https://stackoverflow.com/questions/4391697/find-the-index-of-a-dict-within-a-list-by-matching-the-dicts-value
def build_dict(seq, key):
    return dict((d[key], dict(d, index=index)) for (index, d) in enumerate(seq))

## Set working directory based on configuration file
workdir: config["workdir"]

## GDC token file for authentication
KEYFILE     = config["gdc_token"]

## Cluster metadata (memory, CPU, etc)
CLUSTER_META    = json.load(open(config["cluster_json"]))

## JSON data
CASES 		= json.load(open(config["cases_json"]))
SAMPLES 	= json.load(open(config["samples_json"]))
ALIQUOTS 	= json.load(open(config["aliquots_json"]))
FILES 		= json.load(open(config["files_json"]))
READGROUPS 	= json.load(open(config["readgroups_json"]))
PAIRS 		= json.load(open(config["pairs_json"]))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## JSON processing
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# @sbamin Unless already implemented, we should be explicitly checking json input for 
# 1. non-empty variables, and 
# 2. unique RGID and RGPU but an identical RGSM tags, e.g., https://github.com/TheJacksonLaboratory/glass_wgs_alignment/blob/d72fb20659bd20fddf952d331533b9ffd88d446e/runner/preprocess_fqs.R#L25 
# We can either check it upfront while making json or more preferable to check just before snakemake submits a workflow per case or sample.
# That way, snakemake should STOP with error or emit WARN for non-compliant RG format. 
# This is more practical if input is FQ and not BAM unless we already have RG info for BAM file.
## @barthf : TO-DO

### NOTE NEED TO SPERATE BAM AND FASTQ READGROUPS

## Validate CASES JSON
## CASES should be unique
## CHECK THAT ALL CASE_ID VALUES IN CASES ARE UNIQUE
## TO-DO


## Validate FILES JSON
## FILES -> FILE_UUID should be unique
## CHECK THAT ALL FILE_UUID VALUES IN FILES ARE UNIQUE
## Check that input files exist
## TO-DO


## Validate PAIRS JSON
## PAIR -> PAIR_ID should be unique
## CHECK THAT ALL PAIR_ID VALUES IN PAIR ARE UNIQUE
## IN PROGRESS


## CASES -> DICT
CASES_DICT = build_dict(CASES, "case_id")


## FILES -> DICT
FILES_DICT = build_dict(FILES, "file_uuid")

## ALIQUOTS -> DICT
ALIQUOTS_DICT = build_dict(ALIQUOTS, "aliquot_id")

## SAMPLES -> DICT
SAMPLES_DICT = build_dict(SAMPLES, "sample_id")


## Pair IDs are unique, PAIRS -> DICT
PAIRS_DICT = build_dict(PAIRS, "pair_id")


## Aliquot IDs and BAM/FQ files map 1:1 (or 1:2 for FQ)

ALIQUOT_TO_BAM_PATH = {}
for file in FILES:
    if file["file_format"] == "BAM":
        ALIQUOT_TO_BAM_PATH[ file["aliquot_id"] ] = file["file_path"]


## Dict of aliquots per case
## Dict of aliquots per batch

BATCH_TO_ALIQUOT = {}
CASE_TO_ALIQUOT = {}
for aliquot in ALIQUOTS:
    aliquot["case_id"] = SAMPLES_DICT[ aliquot["sample_id"] ]["case_id"]
    aliquot["project_id"] = CASES_DICT[ aliquot["case_id"] ]["project_id"]
    
    if aliquot["case_id"] not in CASE_TO_ALIQUOT:
        CASE_TO_ALIQUOT[ aliquot["case_id"] ] = [ aliquot["aliquot_id"] ]
    elif aliquot["aliquot_id"] not in CASE_TO_ALIQUOT[ aliquot["case_id"] ]:
        CASE_TO_ALIQUOT[ aliquot["case_id"] ].append(aliquot["aliquot_id"])
    
    if aliquot["project_id"] not in BATCH_TO_ALIQUOT:
        BATCH_TO_ALIQUOT[ aliquot["project_id"] ] = [ aliquot["aliquot_id"] ]
    elif aliquot["aliquot_id"] not in BATCH_TO_ALIQUOT:
        BATCH_TO_ALIQUOT[ aliquot["project_id"] ].append(aliquot["aliquot_id"])



## Aliquots and RGIDs map 1:many

ALIQUOT_TO_RGID = {}        
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_RGID:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ] = [ readgroup["RGID"] ]
    else:
        ALIQUOT_TO_RGID[ readgroup["aliquot_id"] ].append(readgroup["RGID"])


## Batches and normal aliquot IDs map 1:many
## Normal aliquot IDs are repeated across multiple pairs from same case
## Each pair has one normal and one tumor

BATCH_TO_NORMAL = {}
for pair in PAIRS:
    pair["project_id"] = CASES_DICT[ pair["case_id"] ]["project_id"]
    PAIRS_DICT[ pair["pair_id"] ]["project_id"] = pair["project_id"]
    if pair["project_id"] not in BATCH_TO_NORMAL:
        BATCH_TO_NORMAL[ pair["project_id"] ] = [ pair["normal_aliquot_id"] ]
    elif pair["normal_aliquot_id"] not in BATCH_TO_NORMAL[ pair["project_id"] ]:
        BATCH_TO_NORMAL[ pair["project_id"] ].append(pair["normal_aliquot_id"])
        

## Readgroup information and 
## Aliquots and RGIDs map 1:many
## RGIDs are unique within an aliquot
## Aliquot IDs and fastQ files map 1:many

ALIQUOT_TO_READGROUP = {}
ALIQUOT_TO_FQ_PATH = {}
for readgroup in READGROUPS:
    if readgroup["aliquot_id"] not in ALIQUOT_TO_READGROUP:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ] = { readgroup["RGID"] : readgroup }
    else:
        ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ] = readgroup
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_path"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_uuid"] ]["file_path"]
    ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_format"] = FILES_DICT[ ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_uuid"] ]["file_format"]
    if ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_format"] == "FQ":
        if readgroup["aliquot_id"] not in ALIQUOT_TO_FQ_PATH:
            ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ] = {}
        ALIQUOT_TO_FQ_PATH[ readgroup["aliquot_id"] ][ readgroup["RGID"] ] = ALIQUOT_TO_READGROUP[ readgroup["aliquot_id"] ][ readgroup["RGID"] ]["file_path"].split(",")

## List of scatterlist items to iterate over
## Each Mutect2 run spawns 50 jobs based on this scatterlist

WGS_SCATTERLIST = ["temp_{num}_of_50".format(num=str(j+1).zfill(4)) for j in range(50)]

## Load modules
include: "snakemake/download.smk"
include: "snakemake/align.smk"
include: "snakemake/mutect2.smk"
include: "snakemake/vep.smk"
include: "snakemake/lumpy.smk"
include: "snakemake/cnv-gatk.smk"
include: "snakemake/varscan2.smk"
include: "snakemake/fingerprinting.smk"
include: "snakemake/delly.smk"
include: "snakemake/manta.smk"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Master rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule all:
    input: "results/qc/multiqc/multiqc_report.html"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Alignment rule
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule align:
    input: expand("results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam", aliquot_id=ALIQUOTS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Download only rule
## Run snakemake with 'snakemake download_only' to activate
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule download_only:
   input: expand("{file}", file=ALIQUOT_TO_BAM_PATH.values())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mt2:
    input: expand("results/mutect2/vep/{pair_id}.filtered2.anno.maf", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SNV rule (VarScan2)
## Run snakemake with target 'snv'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule vs2:
    input:
        expand("results/varscan2/final/{pair_id}.somatic.hc.filtered.final.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SV preprocessing rule
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule svprepare:
    input:
    	expand("results/lumpy/{aliquot_id}.realn.mdup.bqsr.splitters.sorted.bam", aliquot_id=ALIQUOT_TO_READGROUP.keys()),
    	expand("results/lumpy/{aliquot_id}.realn.mdup.bqsr.discordant.sorted.bam", aliquot_id=ALIQUOT_TO_READGROUP.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling pipeline
## Run snakemake with target 'svprepare'
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnv:
    input:
        expand("results/cnv/callsegments/{pair_id}.called.seg", pair_id=PAIRS_DICT.keys()),
        expand("results/cnv/plotmodeledsegments/{pair_id}/{pair_id}.modeled.png", pair_id=PAIRS_DICT.keys()),
        expand("results/cnv/plotcr/{aliquot_id}/{aliquot_id}.denoised.png", aliquot_id=ALIQUOT_TO_READGROUP.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Delly
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule delly:
    input:
        expand("results/delly/call/{pair_id}.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using LUMPY-SV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule lumpy:
    input:
        expand("results/lumpy/call/{pair_id}.dict.sorted.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SV using Manta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule manta:
    input:
        expand("results/manta/{pair_id}/somaticSV.vcf.gz", pair_id=PAIRS_DICT.keys())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Fingerprinting pipeline
## Check sample and case fingerprints
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprint:
    input:
        expand("results/fingerprinting/sample/{aliquot_id}.crosscheck_metrics", aliquot_id=ALIQUOT_TO_READGROUP.keys()),
        expand("results/fingerprinting/case/{case_id}.crosscheck_metrics", case_id=CASES_DICT.keys()),
        expand("results/fingerprinting/batch/{batch}.crosscheck_metrics", batch=BATCH_TO_ALIQUOT.keys()),
        "results/fingerprinting/GLASS-WG.crosscheck_metrics",
        expand("results/fingerprinting/batch/{batch}.clustered.crosscheck_metrics", batch=BATCH_TO_ALIQUOT.keys()),
        "results/fingerprinting/GLASS-WG.clustered.crosscheck_metrics"

## END ##
