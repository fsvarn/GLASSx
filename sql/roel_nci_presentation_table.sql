WITH dna_surg AS
(
	SELECT gs.case_barcode, SUBSTRING(gs.tumor_barcode_a, 1, 15) AS sample_barcode, surgical_interval, received_tmz, received_rt
	FROM analysis.tumor_clinical_comparison tc
	JOIN analysis.silver_set gs ON gs.tumor_barcode_a = tc.tumor_barcode_a AND gs.tumor_barcode_b = tc.tumor_barcode_b
),
rna_surg AS
(
	SELECT gs.case_barcode, SUBSTRING(gs.tumor_barcode_a, 1, 15) AS sample_barcode, surgical_interval, received_tmz, received_rt
	FROM analysis.tumor_rna_clinical_comparison tc
	JOIN analysis.rna_silver_set gs ON gs.tumor_barcode_a = tc.tumor_barcode_a AND gs.tumor_barcode_b = tc.tumor_barcode_b
),
all_surg AS
(
	SELECT * FROM rna_surg
	UNION
	SELECT * FROM dna_surg
),
uni_surg AS
(
	SELECT DISTINCT * FROM all_surg
),
full_data AS
(
	SELECT us.*, grade, cc.case_overall_survival_mo, case_age_diagnosis_years, case_sex, cs.idh_codel_subtype
	FROM uni_surg us
	JOIN clinical.cases cc ON cc.case_barcode = us.case_barcode
	JOIN clinical.subtypes cs ON cs.case_barcode = us.case_barcode
	JOIN clinical.surgeries su ON su.sample_barcode = us.sample_barcode
)
SELECT idh_codel_subtype, 
COUNT(*) AS patients, 
COUNT(CASE WHEN grade = 'II' THEN 1 END) AS grade2,
COUNT(CASE WHEN grade = 'III' THEN 1 END) AS grade3,
COUNT(CASE WHEN grade = 'IV' THEN 1 END) AS grade4,
AVG(case_age_diagnosis_years) AS age, 
stddev(case_age_diagnosis_years) AS age_sd, 
COUNT(CASE WHEN case_sex = 'male' THEN 1 END)::float/COUNT(CASE WHEN case_sex IS NOT NULL THEN 1 END )::float AS male,
AVG(case_overall_survival_mo) AS surv, 
stddev(case_overall_survival_mo) AS surv_sd, 
AVG(surgical_interval) AS surgical_interval, 
stddev(surgical_interval) AS surgical_interval_sd,
COUNT(CASE WHEN received_tmz THEN 1 END)::float/COUNT(CASE WHEN received_tmz IS NOT NULL THEN 1 END )::float AS temozolomide_treated,
COUNT(CASE WHEN received_rt THEN 1 END)::float/COUNT(CASE WHEN received_rt IS NOT NULL THEN 1 END )::float AS male
FROM full_data
GROUP BY idh_codel_subtype
ORDER BY 1 DESC