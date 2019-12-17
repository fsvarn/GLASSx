CREATE MATERIALIZED VIEW  analysis.gene_tpm AS
SELECT tpm.aliquot_barcode,
etm.gene_symbol,
sum(tpm.est_counts) AS est_counts,
sum(tpm.tpm) AS tpm
FROM analysis.transcript_tpm tpm
LEFT JOIN ref.ensembl_transcript_mapping etm ON tpm.target_id = etm.ensembl_transcript_id
GROUP BY tpm.aliquot_barcode, etm.gene_symbol;