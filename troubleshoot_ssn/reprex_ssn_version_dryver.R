library(SSN2)

in_test_ssn <- SSN2::ssn_import("troubleshoot_ssn/reprex_ssn.ssn", overwrite=TRUE)

SSN2::ssn_create_distmat(in_test_ssn)

ssn_from_scratch <- ssn_lm(
  formula = richness ~ log10(basin_area_km2),
  ssn.object = in_test_ssn, 
  tailup_type = 'linear',
  additive = "afv_qsqrt",
  partition_factor = ~ as.factor(campaign),
  # random = ~ country,
  estmethod = 'ml')


summary(ssn_from_scratch)
