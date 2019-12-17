WITH 
cn_counts AS
(
	SELECT gene_symbol, count(*) AS ct
	FROM variants.gene_copy_number cn
	GROUP BY gene_symbol
),
cn_max AS
(
	SELECT gene_symbol
	FROM cn_counts cn
	INNER JOIN
		(SELECT MAX(ct) AS mx
		FROM cn_counts) mv 
	ON cn.ct = mv.mx
),
genes AS
(
	SELECT gene_symbol
	FROM analysis.gene_tpm 
	
	INTERSECT 
	
	SELECT gene_symbol
	FROM cn_max
),
rand_genes AS
(
	SELECT gene_symbol, 'random' AS signature_name, 'random' AS signature_set
	FROM genes
	ORDER BY RANDOM()
	LIMIT 1000
),
sig_genes AS
(
	SELECT *
	FROM ref.immune_signatures
	WHERE signature_name != 'LGG.purity.correlated'
	
	UNION
	
	SELECT *
	FROM rand_genes
),
tpm_sig_genes AS
(
	SELECT aliquot_barcode, gt.gene_symbol, sg.signature_name, sg.signature_set, gt.tpm
	FROM analysis.gene_tpm gt
	JOIN sig_genes sg ON gt.gene_symbol = sg.gene_symbol
),
cn_sig_genes AS
(
	SELECT aliquot_barcode, cn.gene_symbol,  sg.signature_name, sg.signature_set, cn.wcr
	FROM variants.gene_copy_number cn
	JOIN sig_genes sg ON cn.gene_symbol = sg.gene_symbol
)
SELECT gt.aliquot_barcode, gt.gene_symbol, gt.signature_name, gt.signature_set,  gt.tpm, se.cellularity, cn.wcr
FROM tpm_sig_genes gt
JOIN analysis.analyte_sets an ON an.rna_barcode = gt.aliquot_barcode
JOIN analysis.pairs pa ON pa.tumor_barcode = an.dna_barcode
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = an.dna_barcode OR gs.tumor_barcode_b = an.dna_barcode
JOIN variants.seqz_params se ON se.pair_barcode = pa.pair_barcode
JOIN cn_sig_genes cn ON cn.gene_symbol = gt.gene_symbol AND cn.aliquot_barcode = an.dna_barcode AND cn.signature_set = gt.signature_set AND cn.signature_name = gt.signature_name
ORDER BY 1,2,3