## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert from FASTQ pair to uBAM
## This step eases follow up steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
## Note for this pipeline, GLSS-CU-RNA fastq files need to first be cleaned using the bin/seqkit_cleanup.sh script
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fq2ubam:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        temp("results/rnafingerprint/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam")
    params:
        RGID = lambda wildcards: wildcards.readgroup,                                                     ## ID
        RGPL = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform"), ## Platform
        RGPU = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_platform_unit"), ## Platform unit
        RGLB = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_library"), ## Library
        RGDT = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_timestamp"), ## Date
        RGSM = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_sample_id"), ## Sample ID
        RGCN = lambda wildcards: manifest.getRGTag(wildcards.aliquot_barcode, wildcards.readgroup, "readgroup_center"), ## Center
        mem = CLUSTER_META["fq2ubam"]["mem"]
    threads:
        CLUSTER_META["fq2ubam"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/fq2ubam/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/rnafingerprint/fq2ubam/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Converting FASTQ file to uBAM format\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem} FastqToSam \
            --FASTQ={input.R1} \
            --FASTQ2={input.R2} \
            --OUTPUT={output} \
            --READ_GROUP_NAME=\"{params.RGID}\" \
            --PLATFORM_UNIT=\"{params.RGPU}\" \
            --SAMPLE_NAME=\"{params.RGSM}\" \
            --PLATFORM=\"{params.RGPL}\" \
            --LIBRARY_NAME=\"{params.RGLB}\" \
            --SEQUENCING_CENTER=\"{params.RGCN}\" \
            --SORT_ORDER=queryname \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on uBAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## 6/11/18
## Added --extract parameter, which includes a "summary.txt" file which gives WARN, FAIL, PASS
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fastqc:
    input:
        "results/rnafingerprint/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam"
    output:
        "results/rnafingerprint/fastqc/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned_fastqc.html"
    params:
        dir = "results/rnafingerprint/fastqc/{aliquot_barcode}",
        mem = CLUSTER_META["fastqc"]["mem"],
        walltime = CLUSTER_META["fastqc"]["time"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["fastqc"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/fastqc/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/rnafingerprint/fastqc/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Running FASTQC\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Illumina Adapters
## Add XT tag to read records to mark the 5' start position of adapter sequences
## Adapter sequences are then removed by subsequent steps
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markadapters:
    input:
        "results/rnafingerprint/ubam/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.unaligned.bam"
    output:
        bam = temp("results/rnafingerprint/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.bam"),
        metric = "results/rnafingerprint/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.metrics.txt"
    params:
        mem = CLUSTER_META["markadapters"]["mem"],
        walltime = CLUSTER_META["revertsam"]["time"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["markadapters"]["cpus-per-task"]
    params:
        mem = CLUSTER_META["markadapters"]["mem"]
    log: 
        dynamic("logs/rnafingerprint/markadapters/{aliquot_barcode}.{readgroup}.log")
    benchmark:
        "benchmarks/rnafingerprint/markadapters/{aliquot_barcode}.{readgroup}.txt"
    message:
        "Adding XT tags. This marks Illumina Adapters and allows them to be removed in later steps\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options -Xmx{params.mem} MarkIlluminaAdapters \
            --INPUT={input} \
            --OUTPUT={output.bam} \
            --METRICS={output.metric} \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## (1) BAM to FASTQ
## Converts cleaned-up uBAM to FASTQ format, one FASTQ pair per readgroup
##
## (2) Align reads using BWA-MEM
## This is optimzed for >70 bp PE reads and this pipeline will need update if we intend
## to use it with shorter reads
##
## (3) Merge BAM Alignment
## Restore altered data and apply and adjust meta information lost during alignment
## ie. restore all the original readgroup tags
##
## See: https://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
##
## Update 06/01: Added markadapter metrics as input even though not required. Metrics file
## is saved all the way at the end of markadapters step, and adding it makes sure that the
## input data is good
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule samtofastq:
    input:
        bam = "results/rnafingerprint/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.bam",
        metric = "results/rnafingerprint/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.metrics.txt"
    output:
        r1 = temp("results/rnafingerprint/samtofastq/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.r1.fastq"),
        r2 = temp("results/rnafingerprint/samtofastq/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.r2.fastq")
    threads:
        CLUSTER_META["samtofastq"]["cpus-per-task"]
    conda:
        "../envs/align.yaml"
    params:
        mem = CLUSTER_META["samtofastq"]["mem"],
        walltime = CLUSTER_META["samtofastq"]["time"]
    log: 
        "logs/rnafingerprint/samtofastq/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/rnafingerprint/samtofastq/{aliquot_barcode}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "gatk --java-options {config[samtofastq_java_opt]} SamToFastq \
            --INPUT={input.bam} \
			--FASTQ={output.r1} \
			--SECOND_END_FASTQ={output.r2} \
            --CLIPPING_ATTRIBUTE=XT \
            --CLIPPING_ACTION=2 \
            --INTERLEAVE=false \
            --NON_PF=true \
            --TMP_DIR={config[tempdir]}  \
            > {log} 2>&1"


rule star:
    input:
        r1 = "results/rnafingerprint/samtofastq/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.r1.fastq",
        r2 = "results/rnafingerprint/samtofastq/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.r2.fastq"
    output:
        bam = temp("results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.Aligned.sortedByCoord.out.bam"),
        rpg = "results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.ReadsPerGene.out.tab"
    threads:
        CLUSTER_META["star"]["cpus-per-task"]
    conda:
        "../envs/star.yaml"
    params:
    	outdir = "results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}."
    log: 
        "logs/rnafingerprint/star/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/rnafingerprint/star/{aliquot_barcode}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
        "STAR \
            --runThreadN {threads} \
            --genomeDir /projects/verhaak-lab/varnf/data/ref/star/hg19 \
            --readFilesIn {input.r1} {input.r2} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFileNamePrefix {params.outdir}  \
            > {log} 2>&1"

rule mergebamalignment:
    input:
        aln = "results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.Aligned.sortedByCoord.out.bam",
        unm = "results/rnafingerprint/markadapters/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.markadapters.bam"
    output:
        bam = temp("results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.realn.bam"),
        bai = temp("results/rnafingerprint/star/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.realn.bai")
    threads:
        CLUSTER_META["mergebamalignment"]["cpus-per-task"]
    conda:
        "../envs/align.yaml"
    params:
        mem = CLUSTER_META["mergebamalignment"]["mem"],
        walltime = CLUSTER_META["mergebamalignment"]["time"]
    log: 
        "logs/rnafingerprint/mergebamalignment/{aliquot_barcode}.{readgroup}.log"
    benchmark:
        "benchmarks/rnafingerprint/mergebamalignment/{aliquot_barcode}.{readgroup}.txt"
    message:
        "BAM to FASTQ --> BWA-MEM --> Merge BAM Alignment.\n"
        "The first step converts the reverted BAM to an interleaved FASTQ, removing Illumina "
        "adapters. The output is then piped to BWA-MEM and aligned. Aligned reads are merged "
        "with the original pre-aligned BAM to preserve original metadata, including read groups.\n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Readgroup: {wildcards.readgroup}"
    shell:
         "gatk --java-options {config[mergebamalignment_java_opt]} MergeBamAlignment \
            --ALIGNED_BAM={input.aln} \
            --UNMAPPED_BAM={input.unm} \
            --OUTPUT={output.bam} \
            --REFERENCE_SEQUENCE={config[reference_fasta]} \
            --CREATE_INDEX=true \
            --ADD_MATE_CIGAR=true \
            --CLIP_ADAPTERS=false \
            --CLIP_OVERLAPPING_READS=true \
            --INCLUDE_SECONDARY_ALIGNMENTS=true \
            --MAX_INSERTIONS_OR_DELETIONS=-1 \
            --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            --ATTRIBUTES_TO_RETAIN=XS \
            --TMP_DIR={config[tempdir]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Mark Duplicates & merge readgroups
## This step marks duplicate reads
## See: https://gatkforums.broadinstitute.org/gatk/discussion/2799
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule markduplicates:
    input:
        lambda wildcards: expand("results/rnafingerprint/star/{sample}/{sample}.{rg}.realn.bam", sample=wildcards.aliquot_barcode, rg=manifest.getRGIDs(wildcards.aliquot_barcode))
    output:
        bam = temp("results/rnafingerprint/markduplicates/{aliquot_barcode}.realn.mdup.bam"),
        bai = temp("results/rnafingerprint/markduplicates/{aliquot_barcode}.realn.mdup.bai"),
        metrics = "results/rnafingerprint/markduplicates/{aliquot_barcode}.metrics.txt"
    params:
        max_records = 6000000,
        walltime = lambda wildcards: CLUSTER_META["markduplicates"]["time"],
        mem = lambda wildcards: CLUSTER_META["markduplicates"]["mem"]
    #resources:
    # 	mem = lambda wildcards, attempt: CLUSTER_META["markduplicates"]["mem"] if attempt == 1 else CLUSTER_META["markduplicates"]["mem_if_fail"]
    threads:
        CLUSTER_META["markduplicates"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/markduplicates/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/rnafingerprint/markduplicates/{aliquot_barcode}.txt"
    message:
        "Readgroup-specific BAM files are combined into a single BAM. "
        "Potential PCR duplicates are marked.\n"
        "Sample: {wildcards.aliquot_barcode}"
    run:
        multi_input = " ".join(["--INPUT=" + s for s in input])
        shell("gatk --java-options -Xmx{params.mem} MarkDuplicates \
            {multi_input} \
            --OUTPUT={output.bam} \
            --METRICS_FILE={output.metrics} \
            --CREATE_INDEX=true \
            --TMP_DIR={config[tempdir]} \
            --MAX_RECORDS_IN_RAM={params.max_records} \
            > {log} 2>&1")


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Recalibrate base quality scores
## This steps computes a bare recalibration table
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule baserecalibrator:
    input:
        "results/rnafingerprint/markduplicates/{aliquot_barcode}.realn.mdup.bam"
    output:
        "results/rnafingerprint/bqsr/{aliquot_barcode}.bqsr.txt"
    params:
        mem = CLUSTER_META["baserecalibrator"]["mem"],
        walltime = CLUSTER_META["baserecalibrator"]["time"]
    threads:
        CLUSTER_META["baserecalibrator"]["cpus-per-task"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/rnafingerprint/bqsr/{aliquot_barcode}.recal.log"
    benchmark:
        "benchmarks/rnafingerprint/bqsr/{aliquot_barcode}.recal.txt"
    message:
        "Calculating base recalibration scores.\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} BaseRecalibrator \
            -R {config[reference_fasta]} \
            -I {input} \
            -O {output} \
            --known-sites {config[gnomad_vcf]} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Apply BQSR
## This step applies a base recalibration table to an input BAM file
## Formerly "PrintReads"
## See: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule applybqsr:
    input:
        bam = "results/rnafingerprint/markduplicates/{aliquot_barcode}.realn.mdup.bam",
        bqsr = "results/rnafingerprint/bqsr/{aliquot_barcode}.bqsr.txt"
    output:
        protected("results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    params:
        mem = CLUSTER_META["applybqsr"]["mem"],
        walltime = CLUSTER_META["applybqsr"]["time"]
    threads:
        CLUSTER_META["applybqsr"]["cpus-per-task"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/rnafingerprint/bqsr/{aliquot_barcode}.apply.log"
    benchmark:
        "benchmarks/rnafingerprint/bqsr/{aliquot_barcode}.apply.txt"
    message:
        "Applying base recalibration scores and generating final BAM file\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} ApplyBQSR \
            -R {config[reference_fasta]} \
            -I {input.bam} \
            -OQ true \
            -O {output} \
            -bqsr {input.bqsr} \
            --create-output-bam-md5 true \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Validate BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Final check to ensure no errors in final analysis-ready BAM file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Jun 28: Added "|| true" to for exit code zero. Snakemake deletes output if exit code
## != zero, and ValidateSamFile returns exit code 2 (errors) or 3 (warnings) if notable
## events are found
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule validatebam:
    input:
        ancient("results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        "results/rnafingerprint/validatebam/{aliquot_barcode}.ValidateSamFile.txt"
    params:
        mem = CLUSTER_META["validatebam"]["mem"],
        walltime = CLUSTER_META["validatebam"]["time"]
    threads:
        CLUSTER_META["validatebam"]["cpus-per-task"]
    conda:
        "../envs/align.yaml"
    log:
        "logs/rnafingerprint/validatebam/{aliquot_barcode}.ValidateSamFile.log"
    benchmark:
        "benchmarks/rnafingerprint/validatebam/{aliquot_barcode}.ValidateSamFile.txt"
    message:
        "Validating BAM file\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} ValidateSamFile \
            -I {input} \
            -O {output} \
            -M SUMMARY \
            > {log} 2>&1 \
            || true"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run FASTQC on aligned BAM
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fastqc_bam:
    input:
        "results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/rnafingerprint/fastqc/{aliquot_barcode}/{aliquot_barcode}.aligned_fastqc.html"
    params:
        dir = "results/rnafingerprint/fastqc/{aliquot_barcode}",
        mem = CLUSTER_META["fastqc_bam"]["mem"]
    conda:
        "../envs/align.yaml"
    threads:
        CLUSTER_META["fastqc_bam"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/fastqc-bam/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/rnafingerprint/fastqc-bam/{aliquot_barcode}.txt"
    message:
        "Running FASTQC\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "fastqc \
            --extract \
            -o {params.dir} \
            -f bam \
            {input} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintSample
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This rule uses CrosscheckFingerprints to check that all readgroups in a singlee sample
## come from the same individual
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintsample_rna:
    input:
        "results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
    output:
        "results/rnafingerprint/sample/{aliquot_barcode}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintsample"]["mem"]
    threads:
        CLUSTER_META["fingerprintsample"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/{aliquot_barcode}.fingerprintsample.log"
    benchmark:
        "benchmarks/rnafingerprint/{aliquot_barcode}.fingerprintsample.txt"
    message:
        "Running Picard CrosscheckFingerprints to check that all readgroups in a sample come from the same individual\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map][file]} \
            --INPUT {input} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## FingerprintCase
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Running Picard CrosscheckFingerprints across multiple samples from the sample individual 
## to check for mismatches
## This differs from the fingerprinting.smk file by accounting for RNA and DNA samples
## See: 
## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_fingerprint_CrosscheckFingerprints.php
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule fingerprintcase_all:
    input:
        lambda wildcards: expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getAliquotsByCase(wildcards.case_barcode, analyte='D')), 
        lambda wildcards: expand("results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getAliquotsByCase(wildcards.case_barcode, analyte='R'))
    output:
        "results/rnafingerprint/case/{case_barcode}.crosscheck_metrics"
    params:
        mem = CLUSTER_META["fingerprintcase"]["mem"],
        samples = lambda _, input: " ".join(["--INPUT " + s for s in input])
    threads:
        CLUSTER_META["fingerprintcase"]["cpus-per-task"]
    log:
        "logs/rnafingerprint/{case_barcode}.fingerprintcase.log"
    benchmark:
        "benchmarks/rnafingerprint/{case_barcode}.fingerprintcase.txt"
    message:
        "Running Picard CrosscheckFingerprints across multiple samples from the sample individual to check for mismatches\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} CrosscheckFingerprints \
            --HAPLOTYPE_MAP {config[haplotype_map][file]} \
            --CROSSCHECK_BY SAMPLE \
            {params.samples} \
            --OUTPUT {output} \
            > {log} 2>&1 \
            || true"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Get genic coverage given a BAM file
## URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# 
# rule gencode_coverage:
#     input:
#         "results/rnafingerprint/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"
#     output:
#         "results/rnafingerprint/gencode-coverage/{aliquot_barcode}.gencode-coverage.tsv"
#     params:
#         mem = CLUSTER_META["gencode_coverage"]["mem"]
#     conda:
#         "../envs/bedsam.yaml"
#     threads:
#         CLUSTER_META["gencode_coverage"]["cpus-per-task"]
#     log:
#         "logs/rnafingerprint/gencode-coverage/{aliquot_barcode}.log"
#     benchmark:
#         "benchmarks/rnafingerprint/gencode-coverage/{aliquot_barcode}.txt"
#     message:
#         "Computing coverage using flattened gencode GTF\n"
#         "Sample: {wildcards.aliquot_barcode}"
#     shell:"""
#         set +o pipefail;
#         samtools view -q 10 -b {input} | 
#             bedtools coverage -a {config[gencode_gtf_flat]} -b stdin -d -sorted -g {config[bedtools_genome]} | 
#             bedtools groupby -i stdin -g 1,2,3,4,5 -c 7 -o sum | 
#             sort -k5,5 |
#             bedtools groupby -i stdin -g 5 -c 4,6 -o sum,sum |
#             awk -F\"[+\\t]\" 'BEGIN {{OFS=\"\\t\"}}{{for(i=1;i<(NF-1);i++){{split($i,g,\".\"); print g[1],$(NF-1),$NF}}}}' \
#             > {output} 2> {log}
#         """


## END ##