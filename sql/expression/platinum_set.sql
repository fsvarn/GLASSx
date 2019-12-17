SELECT ss.case_barcode,
ss.tumor_barcode_a AS rna_barcode_a,
ss.tumor_barcode_b AS rna_barcode_b,
an1.dna_barcode AS dna_barcode_a,
an2.dna_barcode AS dna_barcode_b
FROM analysis.rna_silver_set ss
JOIN analysis.analyte_sets an1 ON ss.tumor_barcode_a = an1.rna_barcode
JOIN analysis.analyte_sets an2 ON ss.tumor_barcode_b = an2.rna_barcode
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = an1.dna_barcode AND gs.tumor_barcode_b = an2.dna_barcode
ORDER BY ss.case_barcode;