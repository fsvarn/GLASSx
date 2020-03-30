WITH 
partitioned_aliquots AS
(
	SELECT bl.aliquot_barcode, sa.case_barcode, cs.idh_codel_subtype,
	row_number() OVER (PARTITION BY \"substring\"(bl.aliquot_barcode, 1, 20) ORDER BY (\"substring\"(bl.aliquot_barcode, 21, 3)) DESC) AS priority
	FROM biospecimen.aliquots al
	JOIN analysis.blocklist bl ON bl.aliquot_barcode = al.aliquot_barcode
	JOIN biospecimen.samples sa ON sa.sample_barcode = al.sample_barcode
	JOIN clinical.subtypes cs ON cs.case_barcode = sa.case_barcode
	WHERE bl.coverage_exclusion = 'allow'::bpchar AND cnv_exclusion <> 'block'::bpchar AND al.aliquot_barcode NOT LIKE '%-NB-%' AND al.aliquot_barcode NOT LIKE '%-NM-%' AND al.aliquot_barcode NOT LIKE '%-M1-%'
),
selected_aliquots AS 
(
	SELECT aliquot_barcode, case_barcode, idh_codel_subtype 
	FROM partitioned_aliquots
	WHERE priority = 1
),
selected_genes AS 
(
	SELECT DISTINCT sn.gene_symbol,
	sn.chrom,
	sn.pos,
	sn.alt,
	sn.variant_classification,
	vc.variant_classification_priority,
	sn.protein_change
	FROM variants.anno sn
	JOIN ref.driver_genes ds ON ds.gene_symbol::text = sn.gene_symbol::text
	LEFT JOIN variants.variant_classifications vc ON sn.variant_classification::text = vc.variant_classification::text
	WHERE ds.has_mut IS TRUE AND ((sn.gene_symbol::text <> ALL (ARRAY['TERT'::character varying, 'IDH1'::character varying]::text[])) AND vc.variant_classification_priority IS NOT NULL OR sn.gene_symbol::text = 'TERT'::text AND sn.variant_classification::text = '5''Flank'::text AND (lower(sn.pos) = ANY (ARRAY[1295228, 1295250])))
	ORDER BY 1
), 
selected_genes_samples AS 
(
	SELECT selected_aliquots.aliquot_barcode,
	selected_aliquots.case_barcode,
	selected_genes.gene_symbol,
	selected_genes.chrom,
	lower(selected_genes.pos) AS start_pos,
	upper(selected_genes.pos) - 1 AS end_pos,
	selected_genes.alt
	FROM selected_aliquots,
	selected_genes
),
gene_sample_coverage AS 
(
	SELECT sg.aliquot_barcode,
	sg.case_barcode,
	sg.gene_symbol,
	gc.gene_coverage::double precision / eg.gene_size::double precision AS gene_cov
	FROM ref.ensembl_gene_mapping gm
	JOIN analysis.gencode_coverage gc ON gc.ensembl_gene_id::text = gm.ensembl_gene_id::text
	JOIN ref.ensembl_genes eg ON eg.ensembl_gene_id::text = gm.ensembl_gene_id::text
	JOIN selected_genes_samples sg ON sg.aliquot_barcode = gc.aliquot_barcode AND sg.gene_symbol::text = gm.gene_symbol::text
), 
variants_by_case_and_gene AS 
(
	SELECT anno.gene_symbol,
	gtc.case_barcode,
	gtc.aliquot_barcode,
	gtc.chrom,
	gtc.pos,
	gtc.alt,
	anno.variant_classification,
	sg.protein_change,
	gtc.ssm2_pass_call,
	gtc.ad_alt,
	lag(gtc.ad_alt > 0) OVER w = (gtc.ad_alt > 0) AS is_same_variant,
	row_number() OVER w AS priority
	FROM variants.passgeno gtc
	JOIN selected_aliquots stp ON stp.aliquot_barcode = gtc.aliquot_barcode
	JOIN selected_genes sg ON sg.chrom = gtc.chrom AND sg.pos = gtc.pos AND sg.alt::text = gtc.alt::text
    JOIN variants.anno anno ON anno.variant_id = gtc.variant_id
	WHERE (gtc.ad_ref + gtc.ad_alt) >= 15
	WINDOW w AS (PARTITION BY anno.gene_symbol, gtc.aliquot_barcode ORDER BY (gtc.ssm2_pass_call::integer) DESC, sg.variant_classification_priority, (gtc.ad_ref + gtc.ad_alt) DESC)
),
variants_case_select AS
(
	SELECT vcg.*,
	SUM(ssm2_pass_call::int) OVER w AS case_call
	FROM variants_by_case_and_gene vcg
	WHERE vcg.priority = 1 OR vcg.priority = 2 AND NOT vcg.is_same_variant
	WINDOW w AS (PARTITION BY vcg.case_barcode, vcg.gene_symbol, vcg.chrom, vcg.pos, vcg.alt)
),
sign_genes_by_subtype AS 
(
	SELECT DISTINCT snv_drivers_subtype.gene_symbol,
	snv_drivers_subtype.idh_codel_subtype
	FROM ref.snv_drivers_subtype
), 
variants_by_case AS 
(
	SELECT vg.case_barcode,
	vg.aliquot_barcode,
	cs.idh_codel_subtype,
	btrim(string_agg(
	CASE
		WHEN vg.case_call > 0 AND vg.ad_alt > 0 THEN ((vg.gene_symbol::text || ' '::text) || vg.protein_change::text) || ', '::text
		ELSE ''::text
	END, ''::text), ', '::text) AS driver,
	count(DISTINCT vg.gene_symbol) AS driver_count
	FROM variants_case_select vg
	LEFT JOIN clinical.subtypes cs ON cs.case_barcode = vg.case_barcode
	LEFT JOIN sign_genes_by_subtype sg ON vg.gene_symbol::text = sg.gene_symbol::text AND sg.idh_codel_subtype::text = cs.idh_codel_subtype::text
	GROUP BY vg.case_barcode, vg.aliquot_barcode,  cs.idh_codel_subtype
), 
driver_status AS 
(
	SELECT 
	variants_by_case.case_barcode,
	variants_by_case.aliquot_barcode,
	variants_by_case.idh_codel_subtype,
	variants_by_case.driver,
	variants_by_case.driver_count
	FROM variants_by_case
	),
driver_stability AS 
(
	SELECT
	sa.case_barcode,
	sa.aliquot_barcode,
	sa.idh_codel_subtype,
	CASE
	WHEN ds.driver_count > 0 THEN ds.driver_count
	ELSE 0::bigint
	END AS snv_driver_count,
	CASE
	WHEN ds.driver <> ''::text THEN ds.driver
	ELSE NULL::text
	END AS snv_driver
	FROM driver_status ds
    RIGHT JOIN selected_aliquots sa ON ds.aliquot_barcode = sa.aliquot_barcode
),
selected_cnv_genes AS 
(
	SELECT dg.gene_symbol,
	rg.chrom,
	rg.pos,
	dg.cnv_direction,
	CASE
	WHEN dg.cnv_direction = 1 THEN 'amp'::text
	WHEN dg.cnv_direction = '-1'::integer THEN 'del'::text
	ELSE NULL::text
	END AS effect
	FROM ref.driver_genes dg
	LEFT JOIN ref.genes rg ON rg.gene_symbol::text = dg.gene_symbol::text
	WHERE dg.has_cnv IS TRUE
),
selected_cnv_genes_samples AS 
(
	SELECT selected_aliquots.aliquot_barcode,
	selected_aliquots.case_barcode,
	selected_cnv_genes.gene_symbol,
	selected_cnv_genes.chrom,
	selected_cnv_genes.pos,
	selected_cnv_genes.cnv_direction,
	selected_cnv_genes.effect
	FROM selected_aliquots,
	selected_cnv_genes
),
seg_stats_optimized AS 
(
	SELECT gs.aliquot_barcode,
	LEAST(0.9, gs.neu_fwmean - 2::numeric * gs.neu_fwsd) AS del_thres,
	GREATEST(1.1, gs.neu_fwmean + 2::numeric * gs.neu_fwsd) AS amp_thres,
	CASE
	WHEN gsa.max_loss_arm_wmean < 0.9 AND gsa.max_loss_arm_n >= 3 THEN GREATEST(0::numeric, gsa.max_loss_arm_wmean - 2::numeric * gsa.max_loss_arm_wsd)
	WHEN gs.del_fwmean < 0.9 AND gs.del_n >= 3 THEN GREATEST(0::numeric, gs.del_fwmean - 2::numeric * gs.del_fwsd)
	ELSE NULL::numeric
	END AS hldel_thres,
	CASE
	WHEN gsa.max_loss_arm_wmean < 0.9 AND gsa.max_loss_arm_n >= 3 THEN gsa.max_loss_arm_wmean
	WHEN gs.del_fwmean < 0.9 AND gs.del_n >= 3 THEN gs.del_fwmean
	ELSE NULL::numeric
	END AS hldel_fwmean,
	CASE
	WHEN gsa.max_loss_arm_wmean < 0.9 AND gsa.max_loss_arm_n >= 3 THEN gsa.max_loss_arm_wsd
	WHEN gs.del_fwmean < 0.9 AND gs.del_n >= 3 THEN gs.del_fwsd
	ELSE NULL::numeric
	END AS hldel_fwsd,
	CASE
	WHEN gsa.max_gain_arm_wmean > 1.1 AND gsa.max_gain_arm_n >= 3 THEN gsa.max_gain_arm_wmean + 2::numeric * gsa.max_gain_arm_wsd
	WHEN gs.amp_fwmean > 1.1 AND gs.amp_n >= 3 THEN gs.amp_fwmean + 2::numeric * gs.amp_fwsd
	ELSE NULL::numeric
	END AS hlamp_thres,
	CASE
	WHEN gsa.max_gain_arm_wmean > 1.1 AND gsa.max_gain_arm_n >= 3 THEN gsa.max_gain_arm_wmean
	WHEN gs.amp_fwmean > 1.1 AND gs.amp_n >= 3 THEN gs.amp_fwmean
	ELSE NULL::numeric
	END AS hlamp_fwmean,
	CASE
	WHEN gsa.max_gain_arm_wmean > 1.1 AND gsa.max_gain_arm_n >= 3 THEN gsa.max_gain_arm_wsd
	WHEN gs.amp_fwmean > 1.1 AND gs.amp_n >= 3 THEN gs.amp_fwsd
	ELSE NULL::numeric
	END AS hlamp_fwsd
	FROM analysis.gatk_seg_stats gs
	LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
), 
cnv_by_aliquot_gene AS 
(
	SELECT sgs.aliquot_barcode,
	sgs.case_barcode,
	sgs.gene_symbol,
	sgs.effect,
	round(c1.wcr, 6) AS cr,
	CASE
	WHEN sgs.cnv_direction = '-1'::integer AND (c1.hlvl_call = '-2'::integer) THEN 'HLDEL'::text
	WHEN sgs.cnv_direction = 1 AND (c1.hlvl_call = 2) THEN 'HLAMP'::text
	ELSE NULL::text
	END AS cnv_state,
	round((c1.wcr - ss1.hldel_fwmean) / ss1.hldel_fwsd, 6) AS hldel_zs,
	round((c1.wcr - ss1.hlamp_fwmean) / ss1.hlamp_fwsd, 6) AS hlamp_zs,
	c1.hlvl_call AS cn
	FROM selected_cnv_genes_samples sgs
	LEFT JOIN analysis.gatk_cnv_by_gene c1 ON c1.aliquot_barcode = sgs.aliquot_barcode AND c1.gene_symbol::text = sgs.gene_symbol::text
	LEFT JOIN seg_stats_optimized ss1 ON ss1.aliquot_barcode = sgs.aliquot_barcode
	WHERE (sgs.cnv_direction = '-1'::integer AND c1.hlvl_call = '-2') OR (sgs.cnv_direction = '1'::integer AND c1.hlvl_call = '2')
),
cnv_by_aliquot AS 
(
	SELECT sa.case_barcode,
	sa.aliquot_barcode,
	sa.idh_codel_subtype,
	count(DISTINCT cpg.gene_symbol) AS cnv_driver_count,
	string_agg(cpg.gene_symbol::text || ' '::text || cpg.effect, ', '::text) AS cnv_driver
	FROM cnv_by_aliquot_gene cpg
	RIGHT JOIN selected_aliquots sa ON sa.aliquot_barcode = cpg.aliquot_barcode
	LEFT JOIN ref.gistic_genes gg ON gg.gene_symbol = cpg.gene_symbol::text AND gg.idh_codel_subtype = sa.idh_codel_subtype::text
	GROUP BY sa.case_barcode, sa.aliquot_barcode, sa.idh_codel_subtype
)
SELECT ds.case_barcode,
ds.aliquot_barcode,
ds.idh_codel_subtype,
ds.snv_driver_count,
ds.snv_driver,
ca.cnv_driver_count,
ca.cnv_driver
FROM driver_stability ds
JOIN cnv_by_aliquot ca ON ca.aliquot_barcode = ds.aliquot_barcode;