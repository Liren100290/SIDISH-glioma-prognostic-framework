#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 5B: Current 9-gene CNS/glioma translational drug tiering v1.1
#
# This step is strictly restricted to the current Step 2 v2.6 direction-stable
# 9-gene SIDISH model and the current Step 5 v1.0 sensitivity + CLUE reversal
# outputs. It does not read, audit, mention, or reuse any previous locked 10-gene
# or old drug-candidate results.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
set.seed(20260523)

CONFIG <- list(
  output_dir = "G:/NEW/The2/SIDISH_Step5B_Current9_CNS_Glioma_Translational_Tiering_v1_1",
  step5_dir = "G:/NEW/The2/SIDISH_Step5_Therapeutic_Vulnerability_Drug_Prioritization_v1_0",
  # Current-9 only priority universe.
  # These drugs must be supported by the current Step 5 sensitivity/CLUE outputs
  # or serve as the current reference drug. Previous locked-10 or legacy drug
  # candidates are intentionally excluded.
  priority_drugs = c(
    "Dasatinib", "Selumetinib", "AZD7762", "Bortezomib", "Trametinib",
    "Dabrafenib", "AZD5363", "PD0325901", "MK-1775", "Dinaciclib",
    "XAV939", "Crizotinib", "Ibrutinib", "Taselisib", "AZD6738",
    "Luminespib", "PLX-4720", "OSI-027", "OTX015", "Teniposide",
    "Topotecan", "Docetaxel", "LCL161", "VE-821", "Temozolomide"
  )
)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
safe_dir_create <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); invisible(x) }
load_or_stop <- function(pkg) { if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg); invisible(TRUE) }
load_or_stop("data.table")
fifelse <- data.table::fifelse
copy <- data.table::copy

OUT <- list(
  root = CONFIG$output_dir,
  audit = file.path(CONFIG$output_dir, "00_audit"),
  annotation = file.path(CONFIG$output_dir, "01_CNS_glioma_annotation"),
  tiering = file.path(CONFIG$output_dir, "02_final_tiering"),
  figures = file.path(CONFIG$output_dir, "04_figures"),
  final = file.path(CONFIG$output_dir, "99_final_outputs")
)
invisible(lapply(OUT, safe_dir_create))

write_tsv <- function(dt, file) {
  safe_dir_create(dirname(file))
  data.table::fwrite(data.table::as.data.table(dt), file, sep = "\t", quote = FALSE, na = "NA")
}
write_text <- function(x, file) { safe_dir_create(dirname(file)); writeLines(x, file, useBytes = TRUE) }
read_tsv <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  data.table::fread(file, sep = "\t", data.table = TRUE, showProgress = FALSE, fill = TRUE)
}
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
drug_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("_[0-9]+$", "", x)
  x <- gsub("\\s*\\([^)]*\\)", "", x)
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

###############################################################################
# Inputs
###############################################################################
step5_final_file <- file.path(CONFIG$step5_dir, "99_final_drug_prioritization", "Step5_final_drug_prioritization_table.tsv")
step5_overlap_file <- file.path(CONFIG$step5_dir, "99_final_drug_prioritization", "Step5_integrated_sensitivity_CLUE_overlap_candidates.tsv")
step5_sens_file <- file.path(CONFIG$step5_dir, "03_drug_sensitivity_prediction", "Step5_drug_sensitivity_cohort_consistency_summary.tsv")
step5_clue_file <- file.path(CONFIG$step5_dir, "02_sidish_high_program_reversal", "Step5_CLUE_trt_cp_compound_reversal_summary.tsv")

full <- read_tsv(step5_final_file)
overlap <- read_tsv(step5_overlap_file)
sens <- read_tsv(step5_sens_file)
clue <- read_tsv(step5_clue_file)
for (dt in list(full, overlap, sens, clue)) {
  if ("drug_name" %in% names(dt) && !"drug_key" %in% names(dt)) dt[, drug_key := drug_key(drug_name)]
}

input_status <- data.table::data.table(
  input = c("Step5_final", "Step5_sensitivity_CLUE_overlap", "Step5_GDSC2_sensitivity", "Step5_CLUE_trt_cp"),
  file = c(step5_final_file, step5_overlap_file, step5_sens_file, step5_clue_file),
  exists = c(file.exists(step5_final_file), file.exists(step5_overlap_file), file.exists(step5_sens_file), file.exists(step5_clue_file)),
  n_rows = c(nrow(full), nrow(overlap), nrow(sens), nrow(clue))
)
write_tsv(input_status, file.path(OUT$audit, "Step5B_input_status.tsv"))

###############################################################################
# Curated translational annotation
###############################################################################
ann <- data.table::rbindlist(list(
  data.table::data.table(drug_name="Temozolomide", drug_key="temozolomide", drug_class="alkylating chemotherapy", known_targets="DNA alkylation / MGMT-repair axis", FDA_EMA_or_clinical_status="standard-of-care glioma chemotherapy", BBB_or_CNS_penetration_evidence="strong; established CNS/glioma clinical use", glioma_preclinical_evidence="established glioma reference agent", glioma_clinical_trial_evidence="standard glioma regimen/reference", toxicity_or_feasibility_note="hematologic toxicity; used as reference not novel nomination", CNS_score=3, glioma_score=3, readiness_score=3, feasibility_score=2, mechanism_score=2, translational_note="Positive-control/reference drug for SIDISH-high sensitivity, not a new repositioning claim."),
  data.table::data.table(drug_name="Trametinib", drug_key="trametinib", drug_class="MEK inhibitor", known_targets="MAP2K1/MAP2K2", FDA_EMA_or_clinical_status="approved oncology drug; used with dabrafenib for BRAF-altered tumors; brain tumor trials/clinical reports exist", BBB_or_CNS_penetration_evidence="moderate; CNS activity mostly in molecularly selected MAPK/BRAF-driven tumors", glioma_preclinical_evidence="MAPK pathway glioma rationale", glioma_clinical_trial_evidence="clinical brain tumor experience, strongest when MAPK/BRAF altered", toxicity_or_feasibility_note="rash, cardiomyopathy, ocular toxicity; mutation/pathway selection needed", CNS_score=2, glioma_score=2, readiness_score=3, feasibility_score=2, mechanism_score=3, translational_note="High-confidence conditional MAPK-axis candidate; not universal for all SIDISH-high patients."),
  data.table::data.table(drug_name="Selumetinib", drug_key="selumetinib", drug_class="MEK inhibitor", known_targets="MAP2K1/MAP2K2", FDA_EMA_or_clinical_status="approved/clinical MEK inhibitor with pediatric low-grade glioma experience", BBB_or_CNS_penetration_evidence="moderate; clinically used in CNS tumor settings", glioma_preclinical_evidence="MAPK pathway glioma rationale", glioma_clinical_trial_evidence="brain tumor clinical experience, especially low-grade glioma/NF1-MAPK biology", toxicity_or_feasibility_note="requires MAPK biology context; not established broad adult high-grade glioma treatment", CNS_score=2, glioma_score=2, readiness_score=3, feasibility_score=2, mechanism_score=3, translational_note="High-confidence conditional MEK-axis candidate."),
  data.table::data.table(drug_name="Dabrafenib", drug_key="dabrafenib", drug_class="RAF inhibitor", known_targets="BRAF", FDA_EMA_or_clinical_status="approved oncology drug; relevant to BRAF V600-altered glioma, often with MEK inhibition", BBB_or_CNS_penetration_evidence="moderate; CNS activity in BRAF-mutant tumors", glioma_preclinical_evidence="BRAF/MAPK glioma rationale", glioma_clinical_trial_evidence="BRAF-mutant glioma/brain tumor clinical evidence exists", toxicity_or_feasibility_note="only rational for BRAF-altered subgroup; resistance and paradoxical MAPK activation concerns", CNS_score=2, glioma_score=2, readiness_score=3, feasibility_score=2, mechanism_score=2, translational_note="Conditional candidate requiring BRAF/MAPK alteration evidence."),
  data.table::data.table(drug_name="Topotecan", drug_key="topotecan", drug_class="topoisomerase I inhibitor", known_targets="TOP1/TOP1MT", FDA_EMA_or_clinical_status="approved oncology drug; CNS delivery strategies have been explored", BBB_or_CNS_penetration_evidence="moderate; CNS/CSF delivery evidence stronger than systemic broad BBB efficacy", glioma_preclinical_evidence="DNA damage/topoisomerase rationale", glioma_clinical_trial_evidence="brain tumor studies exist but not standard-of-care", toxicity_or_feasibility_note="myelosuppression; delivery/therapeutic window constraints", CNS_score=2, glioma_score=2, readiness_score=2, feasibility_score=1, mechanism_score=2, translational_note="Translational but toxicity/delivery-limited."),
  data.table::data.table(drug_name="Ibrutinib", drug_key="ibrutinib", drug_class="BTK inhibitor", known_targets="BTK/BLK/BMX", FDA_EMA_or_clinical_status="approved hematologic oncology drug", BBB_or_CNS_penetration_evidence="moderate to strong; CNS lymphoma activity supports CNS exposure", glioma_preclinical_evidence="BTK/myeloid and glioma microenvironment rationale", glioma_clinical_trial_evidence="glioma evidence limited compared with CNS lymphoma", toxicity_or_feasibility_note="bleeding/infection/cardiac risk; glioma efficacy uncertain", CNS_score=2, glioma_score=1, readiness_score=3, feasibility_score=2, mechanism_score=2, translational_note="CNS-feasible repurposing candidate, glioma evidence still limited."),
  data.table::data.table(drug_name="MK-1775", drug_key="mk1775", drug_class="WEE1 inhibitor", known_targets="WEE1", FDA_EMA_or_clinical_status="clinical investigational drug (adavosertib)", BBB_or_CNS_penetration_evidence="limited to moderate; brain tumor/RT combination rationale", glioma_preclinical_evidence="DDR/radiation-sensitization glioma rationale", glioma_clinical_trial_evidence="not established as glioma standard; investigational", toxicity_or_feasibility_note="myelosuppression/GI toxicity; schedule-sensitive", CNS_score=1, glioma_score=2, readiness_score=2, feasibility_score=1, mechanism_score=3, translational_note="DDR candidate for combination/radiation contexts."),
  data.table::data.table(drug_name="AZD6738", drug_key="azd6738", drug_class="ATR inhibitor", known_targets="ATR", FDA_EMA_or_clinical_status="clinical investigational drug (ceralasertib)", BBB_or_CNS_penetration_evidence="limited; CNS evidence not strong enough for Tier 1", glioma_preclinical_evidence="DDR/radiation-sensitization rationale", glioma_clinical_trial_evidence="investigational; glioma-specific readiness limited", toxicity_or_feasibility_note="hematologic toxicity; combination scheduling important", CNS_score=1, glioma_score=1, readiness_score=2, feasibility_score=1, mechanism_score=3, translational_note="Mechanistically attractive DDR candidate, CNS filter not yet sufficient."),
  data.table::data.table(drug_name="Dasatinib", drug_key="dasatinib", drug_class="SRC/ABL/KIT/PDGFR kinase inhibitor", known_targets="SRC/ABL/KIT/PDGFR/EPH axis", FDA_EMA_or_clinical_status="approved oncology drug; tested in glioma settings", BBB_or_CNS_penetration_evidence="limited by efflux/BBB concerns despite some CNS exposure contexts", glioma_preclinical_evidence="SRC/PDGFR glioma rationale", glioma_clinical_trial_evidence="glioma trials have not established broad efficacy", toxicity_or_feasibility_note="strong computation but major CNS/clinical caveat", CNS_score=1, glioma_score=1, readiness_score=3, feasibility_score=1, mechanism_score=3, translational_note="Strong computational candidate but downgraded by BBB/clinical limitations."),
  data.table::data.table(drug_name="Bortezomib", drug_key="bortezomib", drug_class="proteasome/NF-kB inhibitor", known_targets="proteasome subunits/NFKB axis", FDA_EMA_or_clinical_status="approved oncology drug", BBB_or_CNS_penetration_evidence="poor to limited CNS penetration", glioma_preclinical_evidence="proteasome/NF-kB glioma rationale", glioma_clinical_trial_evidence="not established for glioma", toxicity_or_feasibility_note="neuropathy and CNS exposure limitations", CNS_score=0, glioma_score=1, readiness_score=3, feasibility_score=0, mechanism_score=2, translational_note="Strong computation but CNS limitation prevents high tier."),
  data.table::data.table(drug_name="AZD7762", drug_key="azd7762", drug_class="CHK1/CHK2 inhibitor", known_targets="CHEK1/CHEK2", FDA_EMA_or_clinical_status="discontinued/early investigational", BBB_or_CNS_penetration_evidence="not adequate for CNS nomination", glioma_preclinical_evidence="DDR rationale only", glioma_clinical_trial_evidence="no strong glioma translational basis", toxicity_or_feasibility_note="development/toxicity limitations; not a practical candidate", CNS_score=0, glioma_score=1, readiness_score=0, feasibility_score=0, mechanism_score=3, translational_note="Computational hit but development feasibility weak."),
  data.table::data.table(drug_name="MG-132", drug_key="mg132", drug_class="proteasome inhibitor tool compound", known_targets="proteasome", FDA_EMA_or_clinical_status="research tool, not clinically deployable", BBB_or_CNS_penetration_evidence="not clinically relevant", glioma_preclinical_evidence="mechanistic tool evidence only", glioma_clinical_trial_evidence="not a clinical candidate", toxicity_or_feasibility_note="tool compound; do not nominate clinically", CNS_score=0, glioma_score=0, readiness_score=0, feasibility_score=0, mechanism_score=2, translational_note="Use only as pathway support for proteasome vulnerability."),
  data.table::data.table(drug_name="XAV939", drug_key="xav939", drug_class="tankyrase/Wnt pathway inhibitor tool compound", known_targets="TNKS/TNKS2", FDA_EMA_or_clinical_status="research tool/preclinical", BBB_or_CNS_penetration_evidence="not established", glioma_preclinical_evidence="Wnt/TCF axis rationale linked to TCF7L2", glioma_clinical_trial_evidence="not clinical", toxicity_or_feasibility_note="tool compound; pathway-level support only", CNS_score=0, glioma_score=1, readiness_score=0, feasibility_score=0, mechanism_score=3, translational_note="Mechanistically useful because TCF/Wnt axis maps to locked model, but not a clinical drug."),
  data.table::data.table(drug_name="Crizotinib", drug_key="crizotinib", drug_class="ALK/MET/ROS1 inhibitor", known_targets="ALK/MET/ROS1", FDA_EMA_or_clinical_status="approved oncology drug", BBB_or_CNS_penetration_evidence="limited compared with newer CNS-active ALK/ROS1 inhibitors", glioma_preclinical_evidence="MET/ALK fusion or amplification rationale only", glioma_clinical_trial_evidence="conditional target-driven use, not broad glioma", toxicity_or_feasibility_note="requires target alteration; CNS-active alternatives may be better", CNS_score=1, glioma_score=1, readiness_score=3, feasibility_score=1, mechanism_score=2, translational_note="Conditional kinase candidate; not broad SIDISH-high nomination."),
  data.table::data.table(drug_name="AZD5363", drug_key="azd5363", drug_class="AKT inhibitor", known_targets="AKT1/AKT2/AKT3", FDA_EMA_or_clinical_status="related clinical class; capivasertib oncology use", BBB_or_CNS_penetration_evidence="limited CNS evidence", glioma_preclinical_evidence="PI3K/AKT glioma pathway rationale", glioma_clinical_trial_evidence="not established in glioma", toxicity_or_feasibility_note="PI3K/AKT class toxicity; CNS uncertainty", CNS_score=1, glioma_score=1, readiness_score=2, feasibility_score=1, mechanism_score=2, translational_note="Pathway candidate but CNS/glioma evidence insufficient for Tier 1."),
  data.table::data.table(drug_name="PD0325901", drug_key="pd0325901", drug_class="MEK inhibitor", known_targets="MAP2K1/MAP2K2", FDA_EMA_or_clinical_status="investigational/older MEK inhibitor", BBB_or_CNS_penetration_evidence="limited compared with clinically used MEK inhibitors", glioma_preclinical_evidence="MAPK pathway rationale", glioma_clinical_trial_evidence="less compelling than trametinib/selumetinib", toxicity_or_feasibility_note="not preferred over more developed MEK inhibitors", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=3, translational_note="MEK class support; lower priority than trametinib/selumetinib."),
  data.table::data.table(drug_name="Dinaciclib", drug_key="dinaciclib", drug_class="CDK inhibitor", known_targets="CDK1/CDK2/CDK5/CDK9", FDA_EMA_or_clinical_status="investigational", BBB_or_CNS_penetration_evidence="limited", glioma_preclinical_evidence="cell-cycle/transcriptional CDK rationale", glioma_clinical_trial_evidence="not established", toxicity_or_feasibility_note="myelosuppression and narrow window", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=2, translational_note="Computational and mechanistic candidate; translational evidence limited."),
  data.table::data.table(drug_name="Taselisib", drug_key="taselisib", drug_class="PI3K inhibitor", known_targets="PIK3CA/PIK3CB/PIK3CD/PIK3CG", FDA_EMA_or_clinical_status="investigational; development limited by toxicity", BBB_or_CNS_penetration_evidence="limited", glioma_preclinical_evidence="PI3K pathway rationale", glioma_clinical_trial_evidence="not established in glioma", toxicity_or_feasibility_note="toxicity/development caveats", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=0, mechanism_score=2, translational_note="Pathway hit, not high-confidence translational candidate."),
  data.table::data.table(drug_name="Luminespib", drug_key="luminespib", drug_class="HSP90 inhibitor", known_targets="HSP90 family", FDA_EMA_or_clinical_status="investigational", BBB_or_CNS_penetration_evidence="uncertain/limited", glioma_preclinical_evidence="stress/chaperone pathway rationale", glioma_clinical_trial_evidence="not established", toxicity_or_feasibility_note="ocular/GI and class toxicity concerns", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=2, translational_note="Computational candidate with limited CNS maturity."),
  data.table::data.table(drug_name="OTX015", drug_key="otx015", drug_class="BET/bromodomain inhibitor", known_targets="BRD2/BRD3/BRD4", FDA_EMA_or_clinical_status="investigational", BBB_or_CNS_penetration_evidence="limited; some BET inhibitors have brain-tumor rationale but agent-specific evidence required", glioma_preclinical_evidence="epigenetic/transcriptional glioma rationale", glioma_clinical_trial_evidence="not established", toxicity_or_feasibility_note="hematologic/GI toxicity; agent-specific CNS unknown", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=2, translational_note="Epigenetic program-reversal candidate; needs CNS validation."),
  data.table::data.table(drug_name="Teniposide", drug_key="teniposide", drug_class="topoisomerase II inhibitor", known_targets="TOP2A/TOP2B", FDA_EMA_or_clinical_status="older chemotherapy", BBB_or_CNS_penetration_evidence="limited/systemic CNS feasibility uncertain", glioma_preclinical_evidence="DNA damage rationale", glioma_clinical_trial_evidence="not modern standard glioma therapy", toxicity_or_feasibility_note="myelosuppression; lower translational attractiveness", CNS_score=1, glioma_score=1, readiness_score=2, feasibility_score=1, mechanism_score=2, translational_note="DNA-damage candidate but less attractive than TMZ/topotecan/DDR targeted options."),
  data.table::data.table(drug_name="Docetaxel", drug_key="docetaxel", drug_class="taxane", known_targets="microtubules", FDA_EMA_or_clinical_status="approved oncology drug", BBB_or_CNS_penetration_evidence="poor BBB penetration/efflux-limited", glioma_preclinical_evidence="nonspecific cytotoxic rationale", glioma_clinical_trial_evidence="not established glioma therapy", toxicity_or_feasibility_note="poor CNS delivery; deprioritize", CNS_score=0, glioma_score=0, readiness_score=3, feasibility_score=0, mechanism_score=1, translational_note="Not prioritized for SIDISH-high glioma despite possible cytotoxicity."),
  data.table::data.table(drug_name="LCL161", drug_key="lcl161", drug_class="SMAC/IAP antagonist", known_targets="IAP family", FDA_EMA_or_clinical_status="investigational", BBB_or_CNS_penetration_evidence="limited", glioma_preclinical_evidence="apoptosis sensitization rationale", glioma_clinical_trial_evidence="not established", toxicity_or_feasibility_note="combination-dependent; CNS unknown", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=2, translational_note="Exploratory apoptosis-sensitization candidate."),
  data.table::data.table(drug_name="VE-821", drug_key="ve821", drug_class="ATR inhibitor tool compound", known_targets="ATR", FDA_EMA_or_clinical_status="preclinical tool compound", BBB_or_CNS_penetration_evidence="not clinically established", glioma_preclinical_evidence="DDR/radiation rationale", glioma_clinical_trial_evidence="not clinical; related ATR inhibitors should be considered instead", toxicity_or_feasibility_note="tool compound, not final clinical candidate", CNS_score=0, glioma_score=1, readiness_score=0, feasibility_score=0, mechanism_score=3, translational_note="Use as DDR pathway support only."),
  data.table::data.table(drug_name="PLX-4720", drug_key="plx4720", drug_class="BRAF inhibitor research analog", known_targets="BRAF", FDA_EMA_or_clinical_status="preclinical/research analog; clinical analogs preferred", BBB_or_CNS_penetration_evidence="not preferred for CNS nomination", glioma_preclinical_evidence="BRAF pathway rationale", glioma_clinical_trial_evidence="clinical use should focus on approved BRAF inhibitors", toxicity_or_feasibility_note="do not nominate PLX-4720 itself", CNS_score=0, glioma_score=1, readiness_score=0, feasibility_score=0, mechanism_score=2, translational_note="Pathway support for BRAF/RAF, not clinical candidate."),
  data.table::data.table(drug_name="OSI-027", drug_key="osi027", drug_class="mTOR inhibitor", known_targets="MTOR", FDA_EMA_or_clinical_status="investigational", BBB_or_CNS_penetration_evidence="limited", glioma_preclinical_evidence="PI3K/AKT/mTOR glioma rationale", glioma_clinical_trial_evidence="not established", toxicity_or_feasibility_note="class limitations and CNS uncertainty", CNS_score=1, glioma_score=1, readiness_score=1, feasibility_score=1, mechanism_score=2, translational_note="Exploratory pathway candidate.")
), fill = TRUE)

ann[, total_translational_score := CNS_score + glioma_score + readiness_score + feasibility_score + mechanism_score]
write_tsv(ann, file.path(OUT$annotation, "Step5B_curated_CNS_glioma_translational_annotation.tsv"))

source_audit <- data.table::data.table(
  evidence_scope = c(
    "method_CLUE_L1000",
    "method_oncoPredict",
    "reference_temozolomide",
    "tier1_MAPK_BRAF",
    "tier1_selumetinib",
    "tier2_dasatinib_caveat",
    "tier2_ibrutinib_CNS_exposure",
    "tier2_MK1775_glioblastoma_window",
    "tier3_AZD6738_DDR"
  ),
  drug_or_method = c(
    "CLUE/L1000",
    "oncoPredict/GDSC2",
    "Temozolomide",
    "Dabrafenib + Trametinib",
    "Selumetinib",
    "Dasatinib",
    "Ibrutinib",
    "MK-1775/adavosertib",
    "AZD6738/ceralasertib"
  ),
  evidence_role = c(
    "Perturbational reversal resource for the current 9-gene SIDISH-high transcriptional program.",
    "Cell-line pharmacogenomic sensitivity prediction resource; hypothesis-generating only.",
    "Glioblastoma reference chemotherapy and positive-control drug.",
    "CNS/glioma clinical actionability for BRAF V600E low-grade glioma; supports conditional MAPK-axis tiering.",
    "Clinical low-grade glioma MEK-inhibitor evidence; supports conditional MEK-axis tiering.",
    "Explains downgrade despite strong computation: recurrent GBM trial and brain efflux limitation.",
    "Supports CNS feasibility but not glioma-specific efficacy.",
    "Supports glioblastoma combination/window evidence for WEE1 inhibition.",
    "Supports DDR mechanism and clinical investigational status; glioma/CNS evidence remains incomplete."
  ),
  source_label = c(
    "Subramanian et al., Cell 2017 / clue.io",
    "Maeser et al., Brief Bioinform 2021",
    "Stupp et al., NEJM 2005",
    "FDA approval notice, 2023-03-16",
    "Fangusaro et al., Lancet Oncol 2019",
    "RTOG 0627 recurrent GBM trial",
    "PCNSL ibrutinib clinical/CNS literature",
    "Adult Brain Tumor Consortium phase I/window study",
    "NCI formulary and ceralasertib clinical literature"
  ),
  url = c(
    "https://pubmed.ncbi.nlm.nih.gov/29195078/",
    "https://pubmed.ncbi.nlm.nih.gov/34260682/",
    "https://pubmed.ncbi.nlm.nih.gov/15758009/",
    "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-dabrafenib-trametinib-pediatric-patients-low-grade-glioma-braf-v600e-mutation",
    "https://pubmed.ncbi.nlm.nih.gov/31151904/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC5762006/",
    "https://pubmed.ncbi.nlm.nih.gov/32677095/",
    "https://pubmed.ncbi.nlm.nih.gov/39820282/",
    "https://dctd.cancer.gov/drug-discovery-development/reagents-materials/formulary/about/agents/ceralasertib"
  )
)
write_tsv(source_audit, file.path(OUT$annotation, "Step5B_CNS_glioma_annotation_source_audit.tsv"))

###############################################################################
# Merge current Step 5 evidence with curated translational annotation
###############################################################################
priority <- data.table::data.table(drug_name_input = CONFIG$priority_drugs)
priority[, drug_key := drug_key(drug_name_input)]

best_current <- full[drug_key %in% priority$drug_key]
if (nrow(best_current) > 0) {
  best_current[, sensitivity_rank_score := fifelse(!is.na(mean_standardized_delta), -mean_standardized_delta, NA_real_)]
  data.table::setorder(best_current, drug_key, -n_high_more_sensitive, best_connectivity_score, min_wilcox_p)
  best_current <- best_current[, .SD[1], by = drug_key]
}
missing_current <- priority[!drug_key %in% best_current$drug_key]
if (nrow(missing_current) > 0) {
  missing_rows <- missing_current[, .(
    drug_key,
    drug_name = drug_name_input,
    sensitivity_evidence = "not_available",
    reversal_evidence = "not_available",
    high_risk_sensitivity_direction = NA_character_,
    n_high_more_sensitive = NA_integer_,
    n_high_less_sensitive = NA_integer_,
    mean_standardized_delta = NA_real_,
    min_wilcox_p = NA_real_,
    best_connectivity_score = NA_real_,
    clue_moa = NA_character_,
    clue_target_name = NA_character_
  )]
  best_current <- data.table::rbindlist(list(best_current, missing_rows), fill = TRUE)
}

tier <- merge(best_current, ann, by = "drug_key", all.x = TRUE, suffixes = c("", "_annotation"))
tier[, curated_drug_name := fifelse(!is.na(drug_name_annotation), drug_name_annotation, drug_name)]
tier[, has_sensitivity_support := sensitivity_evidence == "supportive" & high_risk_sensitivity_direction == "SIDISH_high_predicted_more_sensitive"]
tier[, has_strong_sensitivity := has_sensitivity_support & !is.na(n_high_more_sensitive) & n_high_more_sensitive >= 4]
tier[, has_reversal_support := reversal_evidence == "supportive_negative_connectivity" & !is.na(best_connectivity_score) & best_connectivity_score < 0]
tier[, has_strong_reversal := has_reversal_support & best_connectivity_score <= -1.25]
tier[, has_CNS_glioma_support := CNS_score >= 2 & glioma_score >= 2]
tier[, has_clinical_feasibility := readiness_score >= 2 & feasibility_score >= 1]

tier[, final_tier_step5B := fifelse(
  drug_key == "temozolomide",
  "Tier_1R_reference_standard_of_care",
  fifelse(
    has_sensitivity_support & has_reversal_support & has_CNS_glioma_support & has_clinical_feasibility & total_translational_score >= 11,
    "Tier_1_conditional_translational_candidate",
    fifelse(
      has_sensitivity_support & has_reversal_support & (CNS_score >= 1 | glioma_score >= 1) & readiness_score >= 1,
      "Tier_2_computational_plus_partial_translational_support",
      fifelse(
        (has_sensitivity_support | has_reversal_support | total_translational_score >= 8) & !(readiness_score == 0 & feasibility_score == 0),
        "Tier_3_exploratory_or_external_rationale",
        "Tier_4_deprioritized_or_not_clinically_actionable"
      )
    )
  )
)]

tier[drug_key %in% c("azd7762", "mg132", "xav939", "ve821", "plx4720"),
     final_tier_step5B := "Tier_4_deprioritized_or_not_clinically_actionable"]
tier[drug_key %in% c("bortezomib", "taselisib") & final_tier_step5B == "Tier_2_computational_plus_partial_translational_support",
     final_tier_step5B := "Tier_3_exploratory_or_external_rationale"]
tier[drug_key == "docetaxel",
     final_tier_step5B := "Tier_4_deprioritized_or_not_clinically_actionable"]

tier[, final_recommendation := fifelse(
  final_tier_step5B == "Tier_1R_reference_standard_of_care",
  "Use as positive-control/reference drug; not a new SIDISH-specific repositioning claim.",
  fifelse(
    final_tier_step5B == "Tier_1_conditional_translational_candidate",
    "Prioritize for manuscript nomination with explicit molecular/CNS caveats.",
    fifelse(
      final_tier_step5B == "Tier_2_computational_plus_partial_translational_support",
      "Keep as computational translational candidate; requires CNS/glioma validation before Tier 1.",
      fifelse(
        final_tier_step5B == "Tier_3_exploratory_or_external_rationale",
        "Mention only as exploratory or external-rationale candidate, not final SIDISH nomination.",
        "Do not prioritize in final therapeutic nomination."
      )
    )
  )
)]
tier[, final_nomination_flag := fifelse(
  final_tier_step5B %in% c("Tier_1R_reference_standard_of_care", "Tier_1_conditional_translational_candidate"),
  "main_text_candidate_or_reference",
  fifelse(
    final_tier_step5B == "Tier_2_computational_plus_partial_translational_support",
    "supplementary_candidate",
    fifelse(final_tier_step5B == "Tier_3_exploratory_or_external_rationale", "exploratory_only", "do_not_nominate")
  )
)]

data.table::setorder(tier, final_tier_step5B, -has_strong_sensitivity, best_connectivity_score, -total_translational_score)
write_tsv(tier, file.path(OUT$tiering, "Step5B_priority_drug_final_tiering.tsv"))

overlap_tiered <- merge(overlap, ann[, .(drug_key, CNS_score, glioma_score, readiness_score, feasibility_score, mechanism_score, total_translational_score, BBB_or_CNS_penetration_evidence, glioma_clinical_trial_evidence, toxicity_or_feasibility_note, translational_note)], by = "drug_key", all.x = TRUE)
overlap_tiered[, final_tier_step5B := fifelse(
  CNS_score >= 2 & glioma_score >= 2 & readiness_score >= 2 & feasibility_score >= 1 & n_high_more_sensitive >= 4 & best_connectivity_score <= -1.25,
  "Tier_1_conditional_translational_candidate",
  fifelse(!is.na(CNS_score), "Tier_2_or_3_after_curated_annotation", "Tier_2_unannotated_computational_candidate")
)]
write_tsv(overlap_tiered, file.path(OUT$tiering, "Step5B_all_sensitivity_CLUE_overlap_with_CNS_tiering.tsv"))


summary <- data.table::data.table(
  metric = c(
    "n_priority_drugs_reviewed",
    "n_current_overlap_candidates",
    "n_Tier_1R_reference",
    "n_Tier_1_conditional",
    "n_Tier_2",
    "n_Tier_3",
    "n_Tier_4",
    "current9_only_guardrail",
    "final_guardrail"
  ),
  value = c(
    nrow(tier),
    nrow(overlap),
    tier[final_tier_step5B == "Tier_1R_reference_standard_of_care", .N],
    tier[final_tier_step5B == "Tier_1_conditional_translational_candidate", .N],
    tier[final_tier_step5B == "Tier_2_computational_plus_partial_translational_support", .N],
    tier[final_tier_step5B == "Tier_3_exploratory_or_external_rationale", .N],
    tier[final_tier_step5B == "Tier_4_deprioritized_or_not_clinically_actionable", .N],
    "Only current Step 2 v2.6 9-gene SIDISH model and current Step 5 sensitivity/CLUE evidence are used.",
    "Tier 1 candidates are conditional translational hypotheses, not clinical treatment recommendations."
  )
)
write_tsv(summary, file.path(OUT$final, "SIDISH_Step5B_final_tiering_summary.tsv"))

tier_final <- tier[, .(
  final_tier_step5B,
  curated_drug_name,
  drug_key,
  sensitivity_evidence,
  high_risk_sensitivity_direction,
  n_high_more_sensitive,
  n_high_less_sensitive,
  mean_standardized_delta,
  min_wilcox_p,
  reversal_evidence,
  best_connectivity_score,
  clue_moa,
  clue_target_name,
  drug_class,
  known_targets,
  CNS_score,
  glioma_score,
  readiness_score,
  feasibility_score,
  mechanism_score,
  total_translational_score,
  BBB_or_CNS_penetration_evidence,
  glioma_clinical_trial_evidence,
  toxicity_or_feasibility_note,
  final_nomination_flag,
  final_recommendation,
  translational_note
)]
write_tsv(tier_final, file.path(OUT$final, "SIDISH_Step5B_final_curated_drug_nomination_table.tsv"))
main_manuscript_candidates <- tier_final[final_nomination_flag %in% c("main_text_candidate_or_reference", "supplementary_candidate")]
write_tsv(main_manuscript_candidates, file.path(OUT$final, "SIDISH_Step5B_manuscript_candidate_table.tsv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  tier_levels <- c(
    "Tier_1R_reference_standard_of_care",
    "Tier_1_conditional_translational_candidate",
    "Tier_2_computational_plus_partial_translational_support",
    "Tier_3_exploratory_or_external_rationale",
    "Tier_4_deprioritized_or_not_clinically_actionable"
  )
  fig_tier <- copy(tier)
  fig_tier[, final_tier_step5B := factor(final_tier_step5B, levels = tier_levels)]
  fig_tier[, tier_display := factor(
    fifelse(final_tier_step5B == "Tier_1R_reference_standard_of_care", "Reference",
      fifelse(final_tier_step5B == "Tier_1_conditional_translational_candidate", "Tier 1",
        fifelse(final_tier_step5B == "Tier_2_computational_plus_partial_translational_support", "Tier 2",
          fifelse(final_tier_step5B == "Tier_3_exploratory_or_external_rationale", "Tier 3", "Tier 4")
        )
      )
    ),
    levels = c("Reference", "Tier 1", "Tier 2", "Tier 3", "Tier 4")
  )]
  tier_counts <- fig_tier[, .N, by = final_tier_step5B]
  tier_counts[, tier_display := factor(c("Reference", "Tier 1", "Tier 2", "Tier 3", "Tier 4")[as.integer(final_tier_step5B)], levels = c("Reference", "Tier 1", "Tier 2", "Tier 3", "Tier 4"))]
  p_count <- ggplot2::ggplot(tier_counts, ggplot2::aes(x = tier_display, y = N, fill = tier_display)) +
    ggplot2::geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 3.3) +
    ggplot2::scale_fill_manual(values = c("#4A6FA5", "#2E8B57", "#C49A3A", "#8A6FAD", "#8C4A4A"), drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Number of reviewed drugs", title = "SIDISH Step 5B Final Drug Tiering") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 12, face = "bold"))
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_final_tier_counts.png"), p_count, width = 8, height = 4.8, dpi = 300)
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_final_tier_counts.pdf"), p_count, width = 8, height = 4.8)

  fig_tier[, computational_score := fifelse(!is.na(mean_standardized_delta), pmax(0, -mean_standardized_delta), 0) +
             fifelse(!is.na(best_connectivity_score), pmax(0, -best_connectivity_score), 0) +
             fifelse(!is.na(n_high_more_sensitive) & n_high_more_sensitive >= 4, 0.5, 0)]
  fig_tier[, translational_score := total_translational_score]
  label_dt <- fig_tier[final_tier_step5B %in% tier_levels[1:2] | curated_drug_name %in% c("Dasatinib", "Ibrutinib", "MK-1775", "Topotecan")]
  p_scatter <- ggplot2::ggplot(fig_tier, ggplot2::aes(x = computational_score, y = translational_score, color = tier_display)) +
    ggplot2::geom_point(size = 2.8, alpha = 0.85) +
    ggplot2::geom_text(data = label_dt, ggplot2::aes(label = curated_drug_name), size = 2.7, vjust = -0.75, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = c("#4A6FA5", "#2E8B57", "#C49A3A", "#8A6FAD", "#8C4A4A"), drop = FALSE) +
    ggplot2::labs(x = "Integrated computational evidence score", y = "CNS/glioma translational score", color = "Final tier", title = "Computational Evidence vs Translational Plausibility") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(size = 12, face = "bold"))
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_computational_vs_translational_score.png"), p_scatter, width = 8.8, height = 5.6, dpi = 300)
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_computational_vs_translational_score.pdf"), p_scatter, width = 8.8, height = 5.6)

  heat_dt <- fig_tier[, .(
    curated_drug_name,
    final_tier_step5B,
    Sensitivity = fifelse(has_strong_sensitivity, 3, fifelse(has_sensitivity_support, 2, 0)),
    Reversal = fifelse(has_strong_reversal, 3, fifelse(has_reversal_support, 2, 0)),
    CNS = CNS_score,
    Glioma = glioma_score,
    Readiness = readiness_score,
    Feasibility = feasibility_score,
    Mechanism = mechanism_score
  )]
  data.table::setorder(heat_dt, final_tier_step5B, -Sensitivity, -Reversal, -CNS, -Glioma)
  heat_long <- data.table::melt(heat_dt, id.vars = c("curated_drug_name", "final_tier_step5B"), variable.name = "evidence_axis", value.name = "score")
  heat_long[, curated_drug_name := factor(curated_drug_name, levels = rev(unique(heat_dt$curated_drug_name)))]
  p_heat <- ggplot2::ggplot(heat_long, ggplot2::aes(x = evidence_axis, y = curated_drug_name, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(colors = c("#F2F2F2", "#C8D7E8", "#5D8AB8", "#1E5A8A"), limits = c(0, 3), breaks = 0:3) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Evidence score", title = "SIDISH Step 5B Drug Evidence Matrix") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), plot.title = ggplot2::element_text(size = 12, face = "bold"))
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_drug_evidence_matrix.png"), p_heat, width = 7.2, height = 8.2, dpi = 300)
  ggplot2::ggsave(file.path(OUT$figures, "Step5B_drug_evidence_matrix.pdf"), p_heat, width = 7.2, height = 8.2)
}

readme <- c(
  "SIDISH Step 5B current 9-gene CNS/glioma tiering v1.1",
  "",
  "Purpose: current 9-gene SIDISH-only CNS/glioma translational annotation and conservative drug tiering.",
  "Current-9 guardrail: no previous locked 10-gene drug candidate or legacy legacy drug result is read, audited, mentioned, or reused.",
  "",
  "Tier logic:",
  "Tier_1R_reference_standard_of_care: reference/control drug, not a new nomination.",
  "Tier_1_conditional_translational_candidate: Step 5 sensitivity + CLUE reversal + CNS/glioma evidence + clinical feasibility.",
  "Tier_2_computational_plus_partial_translational_support: strong computation but CNS/glioma or feasibility caveat.",
  "Tier_3_exploratory_or_external_rationale: partial computation or external CNS/glioma rationale only.",
  "Tier_4_deprioritized_or_not_clinically_actionable: direction conflict, tool compound, poor CNS feasibility, or weak development path.",
  "",
  "Main output:",
  "99_final_outputs/SIDISH_Step5B_final_curated_drug_nomination_table.tsv",
  "99_final_outputs/SIDISH_Step5B_manuscript_candidate_table.tsv",
  "",
  "Evidence/source audit:",
  "01_CNS_glioma_annotation/Step5B_CNS_glioma_annotation_source_audit.tsv",
  "",
  "Figure-ready outputs:",
  "04_figures/Step5B_final_tier_counts.png",
  "04_figures/Step5B_computational_vs_translational_score.png",
  "04_figures/Step5B_drug_evidence_matrix.png"
)
write_text(readme, file.path(OUT$root, "README_SIDISH_Step5B_Current9_CNS_Glioma_Translational_Tiering_v1_1.txt"))
write_text(capture.output(utils::sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))
msg("SIDISH Step 5B final CNS/glioma translational tiering finished.")
msg("Output directory: ", OUT$root)
