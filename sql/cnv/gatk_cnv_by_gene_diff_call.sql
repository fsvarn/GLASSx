CREATE MATERIALIZED VIEW analysis.gatk_cnv_by_gene_diff_call
AS
WITH 
seg_diff AS
(
	SELECT pa.tumor_pair_barcode,
	pa.case_barcode,
	pa.tumor_barcode_a,
	pa.tumor_barcode_b,
	s1.chrom,
	s1.pos * s2.pos AS pos,
	s2.log2_copy_ratio - s1.log2_copy_ratio AS delta_log2_copy_ratio
	FROM analysis.tumor_pairs pa
	JOIN variants.gatk_seg s1 ON s1.aliquot_barcode = pa.tumor_barcode_a
	JOIN variants.gatk_seg s2 ON s2.aliquot_barcode = pa.tumor_barcode_b AND s1.chrom = s2.chrom AND s1.pos && s2.pos
),
unfiltered_seg_wmean_wsd AS 
(
	SELECT gs.tumor_pair_barcode,
	gs.case_barcode,
	gs.tumor_barcode_a,
	gs.tumor_barcode_b,
	count(*) AS num_seg,
	sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * (2::numeric ^ gs.delta_log2_copy_ratio)) / sum(upper(gs.pos) - lower(gs.pos) - 1)::numeric AS wmean,
	sqrt((sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * ((2::numeric ^ gs.delta_log2_copy_ratio) ^ 2::numeric)) - (sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * (2::numeric ^ gs.delta_log2_copy_ratio)) ^ 2::numeric) / sum(upper(gs.pos) - lower(gs.pos) - 1)::numeric) / (sum(upper(gs.pos) - lower(gs.pos) - 1) - 1)::numeric) AS wsd
	FROM seg_diff gs
	WHERE (2::numeric ^ gs.delta_log2_copy_ratio) >= 0.9 AND (2::numeric ^ gs.delta_log2_copy_ratio) <= 1.1
	GROUP BY gs.tumor_pair_barcode, gs.case_barcode, gs.tumor_barcode_a, gs.tumor_barcode_b
), 
filtered_seg_wmean_wsd AS 
(
	SELECT gs.tumor_pair_barcode,
	gs.case_barcode,
	gs.tumor_barcode_a,
	gs.tumor_barcode_b,
	us.num_seg,
	us.wmean,
	us.wsd,
	count(*) AS fnum_seg,
	sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * (2::numeric ^ gs.delta_log2_copy_ratio)) / sum(upper(gs.pos) - lower(gs.pos) - 1)::numeric AS fwmean,
	sqrt((sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * ((2::numeric ^ gs.delta_log2_copy_ratio) ^ 2::numeric)) - (sum((upper(gs.pos) - lower(gs.pos) - 1)::numeric * (2::numeric ^ gs.delta_log2_copy_ratio)) ^ 2::numeric) / sum(upper(gs.pos) - lower(gs.pos) - 1)::numeric) / (sum(upper(gs.pos) - lower(gs.pos) - 1) - 1)::numeric) AS fwsd
	FROM seg_diff gs
	JOIN unfiltered_seg_wmean_wsd us ON us.tumor_barcode_a = gs.tumor_barcode_a AND us.tumor_barcode_b = gs.tumor_barcode_b
	WHERE (2::numeric ^ gs.delta_log2_copy_ratio) >= 0.9 AND (2::numeric ^ gs.delta_log2_copy_ratio) <= 1.1 AND ((2::numeric ^ gs.delta_log2_copy_ratio) - us.wmean) > ('-2.0'::numeric * us.wsd) AND ((2::numeric ^ gs.delta_log2_copy_ratio) - us.wmean) < (2.0 * us.wsd)
	GROUP BY gs.tumor_pair_barcode, gs.case_barcode, gs.tumor_barcode_a, gs.tumor_barcode_b, us.num_seg, us.wmean, us.wsd
),
genes AS
(
	SELECT gene_symbol, ge.chrom, ge.pos
	FROM ref.genes ge
	WHERE (ge.chrom <> ALL (ARRAY[23, 24]))
), 
join_genes AS 
(
	SELECT gs.aliquot_barcode,
	ge.gene_symbol,
	ge.chrom,
	ge.pos AS gene_pos,
	ge.pos * gs.pos AS seg_pos,
	upper(ge.pos * gs.pos)::numeric - lower(ge.pos * gs.pos)::numeric - 1::numeric AS seg_size,
	2::numeric ^ gs.log2_copy_ratio AS cr
	FROM genes ge
	LEFT JOIN variants.gatk_seg gs ON gs.chrom = ge.chrom AND gs.pos && ge.pos
), 
weighted_gene_mean AS
(
	SELECT t1.aliquot_barcode,
	t1.gene_symbol,
	t1.chrom,
	t1.seg_pos,
	int4range(min(lower(t1.seg_pos)), max(upper(t1.seg_pos))) AS pos,
	count(*) AS num_seg,
	sum(t1.seg_size * t1.cr) / sum(t1.seg_size) AS wcr,
	CASE
	WHEN count(*) > 1 AND sum(t1.seg_size * (t1.cr ^ 2::numeric)) > ((sum(t1.seg_size * t1.cr) ^ 2::numeric) / sum(t1.seg_size)) THEN sqrt((sum(t1.seg_size * (t1.cr ^ 2::numeric)) - (sum(t1.seg_size * t1.cr) ^ 2::numeric) / sum(t1.seg_size)) / (sum(t1.seg_size) - 1::numeric))
	ELSE 0::numeric
	END AS wsd,
	sum(t1.seg_size) AS seg_size
	FROM join_genes t1
	WHERE t1.cr > 0::numeric AND t1.seg_size > 0::numeric
	GROUP BY t1.aliquot_barcode, t1.gene_symbol, t1.chrom, t1.seg_pos
	ORDER BY t1.aliquot_barcode, t1.gene_symbol, t1.chrom, t1.seg_pos
), 
bin_diff AS 
(
	SELECT pa.tumor_pair_barcode,
	pa.case_barcode,
	pa.tumor_barcode_a,
	pa.tumor_barcode_b,
	s1.gene_symbol,
	s1.chrom,
	s1.pos,
	s1.seg_size,
	log(2::numeric, s2.wcr) - log(2::numeric, s1.wcr) AS delta_log2_copy_ratio
	FROM analysis.tumor_pairs pa
	JOIN weighted_gene_mean s1 ON s1.aliquot_barcode = pa.tumor_barcode_a
	JOIN weighted_gene_mean s2 ON s2.aliquot_barcode = pa.tumor_barcode_b AND s1.chrom = s2.chrom AND s1.pos = s2.pos AND s1.gene_symbol = s2.gene_symbol
),
call_cnv AS 
(
	SELECT gs.tumor_pair_barcode,
	gs.case_barcode,
	gs.tumor_barcode_a,
	gs.tumor_barcode_b,
	gs.gene_symbol,
	gs.chrom,
	gs.pos,
	gs.delta_log2_copy_ratio,
	CASE
		WHEN (2::numeric ^ gs.delta_log2_copy_ratio) >= 0.9 AND (2::numeric ^ gs.delta_log2_copy_ratio) <= 1.1 THEN 0
		WHEN ((2::numeric ^ gs.delta_log2_copy_ratio) - fis.fwmean) < ('-2.0'::numeric * fis.fwsd) THEN '-1'::integer
		WHEN ((2::numeric ^ gs.delta_log2_copy_ratio) - fis.fwmean) > (2.0 * fis.fwsd) THEN 1
	ELSE 0
	END AS cnv_call
	FROM bin_diff gs
	JOIN filtered_seg_wmean_wsd fis ON fis.tumor_barcode_a = gs.tumor_barcode_a AND fis.tumor_barcode_b = gs.tumor_barcode_b
)
SELECT call_cnv.tumor_pair_barcode,
call_cnv.case_barcode,
call_cnv.tumor_barcode_a,
call_cnv.tumor_barcode_b,
call_cnv.gene_symbol,
call_cnv.chrom,
call_cnv.pos,
call_cnv.delta_log2_copy_ratio,
call_cnv.cnv_call
FROM call_cnv
ORDER BY call_cnv.tumor_pair_barcode, call_cnv.chrom, call_cnv.pos;
