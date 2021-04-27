WITH sample_tab AS
(
	SELECT DISTINCT sample_barcode
	FROM biospecimen.aliquots
	WHERE aliquot_barcode NOT LIKE 'TCGA%'	
),
sample_center AS
(
	SELECT SUBSTRING(sample_barcode, 6,2) AS center
	FROM sample_tab 
),
sample_count AS
(
	SELECT cs.case_source_description, COUNT(*) AS sample_count
	FROM sample_center cc
	JOIN clinical.case_sources cs ON cs.case_source = cc.center
	GROUP BY center, cs.case_source_description
	ORDER BY 1
),
case_tab AS
(
	SELECT DISTINCT substring(sample_barcode,1,12) AS case_barcode
	FROM biospecimen.aliquots
	WHERE aliquot_barcode NOT LIKE 'TCGA%'	
),
case_center AS
(
	SELECT SUBSTRING(case_barcode, 6,2) AS center
	FROM case_tab 
),
case_count AS
(
	SELECT cs.case_source_description, COUNT(*) AS case_count
	FROM case_center cc
	JOIN clinical.case_sources cs ON cs.case_source = cc.center
	GROUP BY center, cs.case_source_description
	ORDER BY 1
)
SELECT sample_count.*, case_count.case_count
FROM sample_count
JOIN case_count ON case_count.case_source_description = sample_count.case_source_description