/*
Query for making the tcr_stats materialized view
Calculates Shannon entropy, evenness, and richness from MiXCR TCR output
New version converts evenness NAs to 1, and thus makes 1 the default value for evenness. 
	-- This is necessary as Shannon values of 0 occur when the number of species is 1 or 0. When the value is 1, the evenness should be 1 as well but instead comes out as NA due to a 0 in the denominator.
	-- This table is filtered to only include aliquots with total_tcr > 0. These statistics are meaningless when total_tcr = 0.
*/
WITH tcr_stats AS (
         SELECT tc.aliquot_barcode,
            '-1'::integer::double precision * sum(tc.clonefraction * ln(tc.clonefraction)) AS shannon,
            COALESCE('-1'::integer::double precision * sum(tc.clonefraction * ln(tc.clonefraction)) / NULLIF(ln(count(*)::double precision), 0::double precision),1) AS evenness,
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
   WHERE total_tcr > 0
  ORDER BY tcr_stats.aliquot_barcode;