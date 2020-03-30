/* analysis.gene_tpm materialized view */

 SELECT tpm.aliquot_barcode,
    etm.gene_symbol,
    sum(tpm.est_counts) AS est_counts,
    sum(tpm.tpm) AS tpm
   FROM analysis.transcript_tpm tpm
     LEFT JOIN ref.ensembl_transcript_mapping etm ON tpm.target_id = etm.ensembl_transcript_id
  GROUP BY tpm.aliquot_barcode, etm.gene_symbol;
  
/* analysis.analyte_sets view */
   SELECT s1.case_barcode,
    s1.sample_barcode,
    a1.aliquot_barcode AS dna_barcode,
    a2.aliquot_barcode AS rna_barcode,
    s1.sample_type,
    a1.aliquot_portion
   FROM biospecimen.aliquots a1
     JOIN biospecimen.aliquots a2 ON a2.sample_barcode = a1.sample_barcode AND a2.aliquot_portion = a1.aliquot_portion
     JOIN biospecimen.samples s1 ON a1.sample_barcode = s1.sample_barcode
     JOIN biospecimen.samples s2 ON a2.sample_barcode = s2.sample_barcode
     LEFT JOIN analysis.blocklist b1 ON b1.aliquot_barcode = a1.aliquot_barcode
  WHERE a1.aliquot_analyte_type = 'D'::bpchar AND a2.aliquot_analyte_type = 'R'::bpchar AND b1.fingerprint_exclusion = 'allow'::bpchar AND b1.coverage_exclusion = 'allow'::bpchar
  ORDER BY s1.case_barcode;
  
/* analysis.platinum_set view */
 SELECT ss.case_barcode,
    ss.tumor_barcode_a AS rna_barcode_a,
    ss.tumor_barcode_b AS rna_barcode_b,
    an1.dna_barcode AS dna_barcode_a,
    an2.dna_barcode AS dna_barcode_b
   FROM analysis.rna_silver_set ss
     JOIN analysis.analyte_sets an1 ON ss.tumor_barcode_a = an1.rna_barcode
     JOIN analysis.analyte_sets an2 ON ss.tumor_barcode_b = an2.rna_barcode
     JOIN analysis.gold_set gs ON gs.tumor_barcode_a = an1.dna_barcode AND gs.tumor_barcode_b = an2.dna_barcode
  ORDER BY ss.case_barcode;
  
/* analysis.rna_dna_pairs view */
 WITH selected_tumor_pairs AS (
         SELECT ps.tumor_pair_barcode,
            ps.case_barcode,
            ps.tumor_barcode_a,
            ps.tumor_barcode_b,
            row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
           FROM analysis.rnaseq_pairs ps
          WHERE ps.comparison_type = 'longitudinal'::text
        )
 SELECT sp.case_barcode,
    sp.tumor_barcode_a AS rna_barcode_a,
    sp.tumor_barcode_b AS rna_barcode_b,
    an1.dna_barcode AS dna_barcode_a,
    an2.dna_barcode AS dna_barcode_b,
    sp.tumor_pair_barcode AS rna_pair_barcode,
    ss.tumor_pair_barcode AS dna_pair_barcode
   FROM selected_tumor_pairs sp
     JOIN analysis.analyte_sets an1 ON sp.tumor_barcode_a = an1.rna_barcode
     JOIN analysis.analyte_sets an2 ON sp.tumor_barcode_b = an2.rna_barcode
     JOIN analysis.silver_set ss ON ss.tumor_barcode_a = an1.dna_barcode AND ss.tumor_barcode_b = an2.dna_barcode
  WHERE sp.priority = 1
  ORDER BY sp.case_barcode;
  
/* analysis.rna_silver_set view */
 WITH selected_tumor_pairs AS (
         SELECT ps.tumor_pair_barcode,
            ps.case_barcode,
            ps.tumor_barcode_a,
            ps.tumor_barcode_b,
            row_number() OVER (PARTITION BY ps.case_barcode ORDER BY ps.surgical_interval_mo DESC, ps.portion_a, ps.portion_b, ("substring"(ps.tumor_pair_barcode, 27, 3))) AS priority
           FROM analysis.rnaseq_pairs ps
             LEFT JOIN analysis.rna_blocklist b1 ON b1.aliquot_barcode = ps.tumor_barcode_a::text
             LEFT JOIN analysis.rna_blocklist b2 ON b2.aliquot_barcode = ps.tumor_barcode_b::text
          WHERE ps.comparison_type = 'longitudinal'::text AND ps.sample_type_b <> 'M1'::bpchar AND b1.complexity_exclusion = 'allow'::bpchar::text AND b2.complexity_exclusion = 'allow'::bpchar::text
        )
 SELECT selected_tumor_pairs.tumor_pair_barcode,
    selected_tumor_pairs.case_barcode,
    selected_tumor_pairs.tumor_barcode_a,
    selected_tumor_pairs.tumor_barcode_b
   FROM selected_tumor_pairs
  WHERE selected_tumor_pairs.priority = 1;
  
/* analysis.rnaseq pairs view */
 SELECT "substring"(a1.aliquot_barcode::text, 1, 18) || "substring"(a2.aliquot_barcode::text, 13, 11) AS tumor_pair_barcode,
    s1.case_barcode,
    a1.aliquot_barcode AS tumor_barcode_a,
    a2.aliquot_barcode AS tumor_barcode_b,
    s1.sample_type AS sample_type_a,
    s2.sample_type AS sample_type_b,
    a1.aliquot_portion AS portion_a,
    a2.aliquot_portion AS portion_b,
        CASE
            WHEN u1.surgery_number = u2.surgery_number THEN 'multisector'::text
            ELSE 'longitudinal'::text
        END AS comparison_type,
    u2.surgical_interval_mo - u1.surgical_interval_mo AS surgical_interval_mo
   FROM biospecimen.aliquots a1
     JOIN clinical.surgeries u1 ON u1.sample_barcode = a1.sample_barcode
     JOIN clinical.surgeries u2 ON u2.case_barcode = u1.case_barcode
     JOIN biospecimen.aliquots a2 ON a2.sample_barcode = u2.sample_barcode AND a2.aliquot_analysis_type = a1.aliquot_analysis_type
     JOIN biospecimen.samples s1 ON a1.sample_barcode = s1.sample_barcode
     JOIN biospecimen.samples s2 ON a2.sample_barcode = s2.sample_barcode
  WHERE a1.aliquot_analyte_type = 'R'::bpchar AND a2.aliquot_analyte_type = 'R'::bpchar AND (u1.surgery_number < u2.surgery_number OR a1.aliquot_portion < a2.aliquot_portion AND u1.surgery_number = u2.surgery_number)
  ORDER BY s1.case_barcode, (u2.surgical_interval_mo - u1.surgical_interval_mo);
