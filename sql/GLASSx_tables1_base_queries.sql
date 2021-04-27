--DNA base query

WITH subset AS
(
	SELECT DISTINCT substring(aliquot_barcode, 1, 12) AS sample_barcode
	FROM analysis.blocklist
	WHERE fingerprint_exclusion = 'allow' AND aliquot_barcode NOT LIKE '%-NB-%' --AND aliquot_barcode LIKE 'GLSS-HF-%' AND aliquot_barcode NOT LIKE 'GLSS-HF-2%' AND aliquot_barcode NOT LIKE 'GLSS-HF-3%'
)

SELECT substring(sample_barcode,1,7), count(*) 
FROM subset
GROUP BY 1
ORDER BY 1


--RNA base query
WITH subset AS
(
	SELECT DISTINCT substring(aliquot_barcode, 1, 15
							 ) AS sample_barcode
	FROM analysis.rna_blocklist
	WHERE fingerprint_exclusion = 'allow' AND aliquot_barcode NOT LIKE '%-NB-%' --AND aliquot_barcode LIKE 'GLSS-LU-%'
)

SELECT substring(sample_barcode,1,7), count(*) 
FROM subset
GROUP BY 1
ORDER BY 1