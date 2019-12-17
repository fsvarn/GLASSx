WITH pten AS
(
	SELECT *
	FROM variants.pgeno pg
	WHERE pg.gene_symbol = 'PTEN' AND pg.variant_classification NOT IN ('INTRON','SILENT','THREE_PRIME_UTR','FIVE_PRIME_UTR','SPLICE_SITE')

)
SELECT ps.*, im1.signature_name, 
im1.enrichment_score AS enrichment_score_a, 
im2.enrichment_score AS enrichment_score_b, 
pg.gene_symbol, 
pg.variant_classification, 
pg.chrom, 
pg.pos, 
pg.ref, 
pg.alt, 
pg.mutect2_call_a, 
pg.mutect2_call_b
FROM analysis.platinum_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im2.signature_name = im1.signature_name
LEFT JOIN pten pg ON ps.dna_barcode_a = pg.tumor_barcode_a 
ORDER BY 1,2

