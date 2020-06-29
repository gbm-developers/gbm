test_that('project Stacked.SDM', {
  data(Env)
  data(Occurrences)
  SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 2,
                          Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                          Spcol = 'SPECIES', ensemble.thresh = 0,
                          verbose = FALSE, cores = 0)
  Env_new <- stack(Env[[1]]-0.3,Env[[2]],Env[[3]])
  SSDM_proj <- project(SSDM,Env_new)
  expect_is(SSDM_proj, 'Stacked.SDM')
  expect_false(all(is.na(values(SSDM_proj@diversity.map))))
})
