WITH selected_tumor_pairs AS 
(
	SELECT ps.tumor_pair_barcode,
	ps.case_barcode,
	ps.tumor_barcode_a,
	ps.tumor_barcode_b,
	row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
	FROM analysis.rnaseq_pairs ps
	LEFT JOIN analysis.rna_blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a::text
	LEFT JOIN analysis.rna_blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b::text
	WHERE ps.comparison_type = 'longitudinal'::text AND ps.sample_type_b <> 'M1'::bpchar AND b1.complexity_exclusion = 'allow'::bpchar::text AND b2.complexity_exclusion = 'allow'::bpchar::text
)
SELECT selected_tumor_pairs.tumor_pair_barcode,
selected_tumor_pairs.case_barcode,
selected_tumor_pairs.tumor_barcode_a,
selected_tumor_pairs.tumor_barcode_b
FROM selected_tumor_pairs
WHERE selected_tumor_pairs.priority = 1;