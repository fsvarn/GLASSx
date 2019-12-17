WITH dna_pairs AS
(
	SELECT ps.tumor_pair_barcode,
	ps.case_barcode,
	ps.tumor_barcode_a,
	ps.tumor_barcode_b,
	row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
	FROM analysis.tumor_pairs ps
	WHERE ps.comparison_type = 'longitudinal'::text AND ps.sample_type_b <> 'M1'::bpchar
),
selected_dna_pairs AS
(
	SELECT dna_pairs.tumor_pair_barcode,
	dna_pairs.case_barcode,
	dna_pairs.tumor_barcode_a,
	dna_pairs.tumor_barcode_b
	FROM dna_pairs
	WHERE dna_pairs.priority = 1
),
rna_pairs AS
(
	SELECT ps.tumor_pair_barcode,
	ps.case_barcode,
	ps.tumor_barcode_a,
	ps.tumor_barcode_b,
	row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
	FROM analysis.rnaseq_pairs ps
	WHERE ps.comparison_type = 'longitudinal'::text AND ps.sample_type_b <> 'M1'::bpchar
),
selected_rna_pairs AS
(
	SELECT rna_pairs.tumor_pair_barcode,
	rna_pairs.case_barcode,
	rna_pairs.tumor_barcode_a,
	rna_pairs.tumor_barcode_b
	FROM rna_pairs
	WHERE rna_pairs.priority = 1
)
SELECT ss.case_barcode,
ss.tumor_barcode_a AS rna_barcode_a,
ss.tumor_barcode_b AS rna_barcode_b,
an1.dna_barcode AS dna_barcode_a,
an2.dna_barcode AS dna_barcode_b
FROM selected_rna_pairs ss
JOIN analysis.analyte_sets an1 ON ss.tumor_barcode_a = an1.rna_barcode
JOIN analysis.analyte_sets an2 ON ss.tumor_barcode_b = an2.rna_barcode
JOIN selected_dna_pairs gs ON gs.tumor_barcode_a = an1.dna_barcode AND gs.tumor_barcode_b = an2.dna_barcode
ORDER BY ss.case_barcode