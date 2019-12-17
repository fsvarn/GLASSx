/*
Query for making the tcr_stats materialized view
Calculates Shannon entropy, evenness, and richness from MiXCR TCR output
*/

WITH tcr_stats AS (
	 SELECT tc.aliquot_barcode,
		'-1'::integer::double precision * sum(tc.clonefraction * ln(tc.clonefraction)) AS shannon,
		'-1'::integer::double precision * sum(tc.clonefraction * ln(tc.clonefraction)) / NULLIF(ln(count(*)::double precision), 0::double precision) AS evenness,
		count(tc.clonecount)::integer AS richness,
		COALESCE(sum(tc.clonecount)::integer, 0) AS total_tcr
	   FROM analysis.mixcr_tcr tc
	  GROUP BY tc.aliquot_barcode
	  ORDER BY tc.aliquot_barcode
	)
SELECT tcr_stats.aliquot_barcode,
	CASE
		WHEN tcr_stats.shannon > 0::double precision THEN tcr_stats.shannon
		ELSE 0::double precision
	END AS shannon,
tcr_stats.evenness,
tcr_stats.richness,
tcr_stats.total_tcr
FROM tcr_stats
ORDER BY tcr_stats.aliquot_barcode;