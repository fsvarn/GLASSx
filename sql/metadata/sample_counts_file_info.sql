SELECT al.aliquot_barcode, al.sample_barcode, al.aliquot_uuid_short, al.aliquot_analyte_type, al.aliquot_analysis_type, al.aliquot_portion, al.aliquot_batch, f.file_name, f.file_size, f.file_md5sum, f.file_format
FROM biospecimen.aliquots al
LEFT JOIN analysis.files f ON al.aliquot_barcode = f.aliquot_barcode
WHERE aliquot_analyte_type = 'R' AND al.aliquot_barcode NOT LIKE 'GLSS-MD-%'
ORDER BY 1,8

SELECT al.aliquot_barcode, al.sample_barcode, al.aliquot_uuid_short, al.aliquot_analyte_type, al.aliquot_analysis_type, al.aliquot_portion, al.aliquot_batch, f.file_name, f.file_size, f.file_md5sum, f.file_format
FROM biospecimen.aliquots al
LEFT JOIN analysis.files f ON al.aliquot_barcode = f.aliquot_barcode
WHERE aliquot_analyte_type = 'D'
ORDER BY 1,8