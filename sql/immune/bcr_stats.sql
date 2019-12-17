/*
Query for making the bcr_stats materialized view
Calculates Shannon entropy, evenness, and richness from MiXCR BCR output
*/

WITH bcr_stats AS (
	 SELECT bc.aliquot_barcode,
		'-1'::integer::double precision * sum(bc.clonefraction * ln(bc.clonefraction)) AS shannon,
		'-1'::integer::double precision * sum(bc.clonefraction * ln(bc.clonefraction)) / NULLIF(ln(count(*)::double precision), 0::double precision) AS evenness,
		count(bc.clonecount)::integer AS richness,
		COALESCE(sum(bc.clonecount)::integer, 0) AS total_bcr
	   FROM analysis.mixcr_bcr bc
	  GROUP BY bc.aliquot_barcode
	  ORDER BY bc.aliquot_barcode
	)
SELECT bcr_stats.aliquot_barcode,
	CASE
		WHEN bcr_stats.shannon > 0::double precision THEN bcr_stats.shannon
		ELSE 0::double precision
	END AS shannon,
bcr_stats.evenness,
bcr_stats.richness,
bcr_stats.total_bcr
FROM bcr_stats
ORDER BY bcr_stats.aliquot_barcode