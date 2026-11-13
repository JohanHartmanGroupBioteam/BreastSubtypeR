# Measured PR #228 Linux allowances; every other value and all structure stay exact.
# Historical TCGA B uses the independently verified Windows-B == frozen-A mapping
# reference bridge. OSLO historical A/B is deliberately outside this policy.
issue199_portability_policy <- function() {
    data.frame(
        scope = c(rep("C_vs_approved_C", 18L), "B_mapping_vs_frozen_A"),
        dataset = c(rep("OSLO2EMIT0obj", 5L), rep("TCGABRCAobj", 14L)),
        path = c(
            "/data_input/se_SSP/assay",
            "/res/results/AIMS/all.probs/20",
            "/res/results/AIMS/data.used",
            "/res/results/sspbc/all.probs/24",
            "/res/results/sspbc/data.used",
            "/data_input/se_NC/assay",
            "/res/results/genefu.robust/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/genefu.robust/outList/testData",
            "/res/results/ssBC/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/ssBC/outList/testData",
            "/res/results/ssBC.v2/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/ssBC.v2/outList/testData",
            "/res/results/cIHC/outList/testData",
            "/res/results/cIHC.itr/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/cIHC.itr/outList/testData",
            "/res/results/PCAPAM50/outList/testData",
            "/res/results/AIMS/all.probs/20",
            "/res/results/sspbc/all.probs/24",
            "/data_input/se_NC/assay"),
        atol = c(0, 0, 0, 0, 0, 2^-48, 2^-46, 2^-48, 2^-46, 2^-48, 2^-46, 2^-48, 2^-48, 2^-46, 2^-48, 2^-48, 0, 0, 2^-48),
        rtol = c(2^-51, 2^-44, 2^-51, 2^-44, 2^-51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2^-44, 2^-45, 0),
        stringsAsFactors = FALSE)
}
