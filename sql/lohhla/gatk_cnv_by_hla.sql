 WITH selected_genes AS (
         SELECT dr.gene_symbol,
            dr.chrom,
            dr.pos
           FROM ref.hla_genes dr
        ), gene_seg_intersect AS (
         SELECT gs.aliquot_barcode,
            t0.gene_symbol,
            gs.chrom,
            upper(t0.pos * gs.pos) - lower(t0.pos * gs.pos) - 1 AS w,
            2::numeric ^ gs.log2_copy_ratio AS cr
           FROM variants.gatk_seg gs
             JOIN selected_genes t0 ON t0.chrom = gs.chrom AND t0.pos && gs.pos
        ), gene_sample_call AS (
         SELECT gene_seg_intersect.aliquot_barcode,
            gene_seg_intersect.gene_symbol,
            sum(gene_seg_intersect.w::numeric * gene_seg_intersect.cr) / sum(gene_seg_intersect.w)::numeric AS wcr
           FROM gene_seg_intersect
          GROUP BY gene_seg_intersect.aliquot_barcode, gene_seg_intersect.gene_symbol
        ), seg_stats_optimized AS (
         SELECT gs.aliquot_barcode,
            LEAST(0.9, gs.neu_fwmean - 2::numeric * gs.neu_fwsd) AS del_thres,
            GREATEST(1.1, gs.neu_fwmean + 2::numeric * gs.neu_fwsd) AS amp_thres,
                CASE
                    WHEN gsa.max_loss_arm_wmean < 0.9 AND gsa.max_loss_arm_n >= 3 THEN GREATEST(0::numeric, gsa.max_loss_arm_wmean - 2::numeric * gsa.max_loss_arm_wsd)
                    WHEN gs.del_fwmean < 0.9 AND gs.del_n >= 3 THEN GREATEST(0::numeric, gs.del_fwmean - 2::numeric * gs.del_fwsd)
                    ELSE NULL::numeric
                END AS hldel_thres,
                CASE
                    WHEN gsa.max_gain_arm_wmean > 1.1 AND gsa.max_gain_arm_n >= 3 THEN gsa.max_gain_arm_wmean + 2::numeric * gsa.max_gain_arm_wsd
                    WHEN gs.amp_fwmean > 1.1 AND gs.amp_n >= 3 THEN gs.amp_fwmean + 2::numeric * gs.amp_fwsd
                    ELSE NULL::numeric
                END AS hlamp_thres
           FROM analysis.gatk_seg_stats gs
             LEFT JOIN analysis.gatk_aneuploidy gsa ON gsa.aliquot_barcode = gs.aliquot_barcode
        ), gene_cp AS (
         SELECT ts.aliquot_barcode,
            t0.gene_symbol,
            ts.chrom,
            upper(t0.pos * ts.pos) - lower(t0.pos * ts.pos) - 1 AS w,
            ts.cellular_prevalence AS cp
           FROM variants.titan_seg ts
             JOIN selected_genes t0 ON t0.chrom = ts.chrom AND t0.pos && ts.pos
        ), gene_cp_agg AS (
         SELECT gene_cp.aliquot_barcode,
            gene_cp.gene_symbol,
            COALESCE(sum(gene_cp.w::double precision * gene_cp.cp) / NULLIF(sum(gene_cp.w), 0)::double precision, NULL::double precision) AS wcp
           FROM gene_cp
          GROUP BY gene_cp.aliquot_barcode, gene_cp.gene_symbol
        )
 SELECT gc.aliquot_barcode,
    gc.gene_symbol,
        CASE
            WHEN gc.wcr >= ss.del_thres AND gc.wcr <= ss.amp_thres THEN 0
            WHEN gc.wcr < ss.hldel_thres THEN '-2'::integer
            WHEN gc.wcr < ss.del_thres THEN '-1'::integer
            WHEN gc.wcr > ss.hlamp_thres THEN 2
            WHEN gc.wcr > ss.amp_thres THEN 1
            ELSE NULL::integer
        END AS hlvl_call,
    gc.wcr,
    cp.wcp AS cellular_prevalence
   FROM gene_sample_call gc
     LEFT JOIN seg_stats_optimized ss ON ss.aliquot_barcode = gc.aliquot_barcode
     LEFT JOIN gene_cp_agg cp ON cp.aliquot_barcode = gc.aliquot_barcode AND cp.gene_symbol::text = gc.gene_symbol::text
  ORDER BY (
        CASE
            WHEN gc.wcr >= ss.del_thres AND gc.wcr <= ss.amp_thres THEN 0
            WHEN gc.wcr < ss.hldel_thres THEN '-2'::integer
            WHEN gc.wcr < ss.del_thres THEN '-1'::integer
            WHEN gc.wcr > ss.hlamp_thres THEN 2
            WHEN gc.wcr > ss.amp_thres THEN 1
            ELSE NULL::integer
        END), gc.aliquot_barcode;