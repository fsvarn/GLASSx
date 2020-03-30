 WITH hla_loss AS (
         SELECT ps.case_barcode,
            pa1.tumor_barcode AS tumor_barcode_a,
            pa2.tumor_barcode AS tumor_barcode_b,
            lh1.hla_type1,
            lh1.hla_type2,
            lh1.pval AS lohhla_pval_a,
            lh2.pval AS lohhla_pval_b,
                CASE
                    WHEN lh1.pval < 0.1::double precision AND (lh1.hla_type1_copy_number < 0.5::double precision OR lh1.hla_type2_copy_number < 0.5::double precision) THEN lh1.loss_allele
                    ELSE NULL::text
                END AS loss_a,
                CASE
                    WHEN lh2.pval < 0.1::double precision AND (lh2.hla_type1_copy_number < 0.5::double precision OR lh2.hla_type2_copy_number < 0.5::double precision) THEN lh2.loss_allele
                    ELSE NULL::text
                END AS loss_b,
            cs.idh_codel_subtype
           FROM analysis.lohhla_set ps
             JOIN analysis.pairs pa1 ON pa1.tumor_barcode::text = ps.tumor_barcode_a
             JOIN analysis.pairs pa2 ON pa2.tumor_barcode::text = ps.tumor_barcode_b
             JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
             JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh1.hla_type1 = lh2.hla_type1 AND lh1.hla_type2 = lh2.hla_type2
             JOIN clinical.subtypes cs ON ps.case_barcode = cs.case_barcode::text
          WHERE lh1.coverage_filter = 20::double precision AND lh2.coverage_filter = 20::double precision
        )
 SELECT hla_loss.case_barcode,
    hla_loss.tumor_barcode_a,
    hla_loss.tumor_barcode_b,
    hla_loss.hla_type1,
    hla_loss.hla_type2,
    hla_loss.lohhla_pval_a,
    hla_loss.lohhla_pval_b,
    hla_loss.loss_a,
    hla_loss.loss_b,
        CASE
            WHEN hla_loss.loss_a IS NULL AND hla_loss.loss_b IS NOT NULL THEN 'gain'::text
            WHEN hla_loss.loss_a IS NOT NULL AND hla_loss.loss_b IS NULL THEN 'loss'::text
            WHEN hla_loss.loss_a IS NULL AND hla_loss.loss_b IS NULL THEN 'none'::text
            WHEN hla_loss.loss_a IS NOT NULL AND hla_loss.loss_b IS NOT NULL AND hla_loss.loss_a = hla_loss.loss_b THEN 'stable'::text
            WHEN hla_loss.loss_a IS NOT NULL AND hla_loss.loss_b IS NOT NULL AND hla_loss.loss_a <> hla_loss.loss_b THEN 'switch'::text
            ELSE NULL::text
        END AS hla_loh_change,
    hla_loss.idh_codel_subtype
   FROM hla_loss
  WHERE hla_loss.lohhla_pval_a IS NOT NULL AND hla_loss.lohhla_pval_b IS NOT NULL
  ORDER BY hla_loss.case_barcode, hla_loss.hla_type1;