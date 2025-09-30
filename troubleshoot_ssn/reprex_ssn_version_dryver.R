library(SSN2)
set.seed(123)
sessionInfo()

#Import SSN
if (!file.exists("reprex_ssn.ssn")) {
  unzip("reprex_ssn.ssn.zip")
}
in_test_ssn <- SSN2::ssn_import("reprex_ssn.ssn", overwrite=TRUE)
SSN2::ssn_create_distmat(in_test_ssn)

# With partition factor 
ssn_from_scratch <- ssn_lm(
  formula = richness ~ log10(basin_area_km2),
  ssn.object = in_test_ssn, 
  tailup_type = 'linear',
  additive = "afv_qsqrt",
  partition_factor = ~ as.factor(campaign),
  # random = ~ country,
  estmethod = 'ml')

summary(ssn_from_scratch)

# Without partition factor 
ssn_from_scratch <- ssn_lm(
  formula = richness ~ log10(basin_area_km2),
  ssn.object = in_test_ssn, 
  tailup_type = 'linear',
  additive = "afv_qsqrt",
  # partition_factor = ~ as.factor(campaign),
  # random = ~ country,
  estmethod = 'ml')

summary(ssn_from_scratch)
