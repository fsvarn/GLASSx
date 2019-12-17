SELECT "substring"(a1.aliquot_barcode::text, 1, 18) || "substring"(a2.aliquot_barcode::text, 13, 11) AS tumor_pair_barcode,
s1.case_barcode,
a1.aliquot_barcode AS tumor_barcode_a,
a2.aliquot_barcode AS tumor_barcode_b,
s1.sample_type AS sample_type_a,
s2.sample_type AS sample_type_b,
a1.aliquot_portion AS portion_a,
a2.aliquot_portion AS portion_b,
CASE
	WHEN u1.surgery_number = u2.surgery_number THEN 'multisector'::text
	ELSE 'longitudinal'::text
END AS comparison_type,
u2.surgical_interval_mo - u1.surgical_interval_mo AS surgical_interval_mo
FROM biospecimen.aliquots a1
JOIN clinical.surgeries u1 ON u1.sample_barcode = a1.sample_barcode
JOIN clinical.surgeries u2 ON u2.case_barcode = u1.case_barcode
JOIN biospecimen.aliquots a2 ON a2.sample_barcode = u2.sample_barcode AND a2.aliquot_analysis_type = a1.aliquot_analysis_type
JOIN biospecimen.samples s1 ON a1.sample_barcode = s1.sample_barcode
JOIN biospecimen.samples s2 ON a2.sample_barcode = s2.sample_barcode
WHERE a1.aliquot_analyte_type = 'R'::bpchar AND a2.aliquot_analyte_type = 'R'::bpchar AND (u1.surgery_number < u2.surgery_number OR a1.aliquot_portion < a2.aliquot_portion AND u1.surgery_number = u2.surgery_number)
ORDER BY s1.case_barcode, (u2.surgical_interval_mo - u1.surgical_interval_mo);