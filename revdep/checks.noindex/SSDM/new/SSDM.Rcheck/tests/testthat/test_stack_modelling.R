test_that('stack modelling function', {
  data(Env)
  data(Occurrences)
  SSDM <- stack_modelling(c('CTA', 'SVM'), Occurrences, Env, rep = 1,
                         Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                         Spcol = 'SPECIES', ensemble.thresh = 0,
                         verbose = FALSE, cores = 0)
  expect_is(SSDM, 'Stacked.SDM')
})
