test_that('ensemble modelling function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  ESDM <- ensemble_modelling(c('CTA', 'MARS'), Occurrences, Env, rep = 1,
                            Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                            ensemble.thresh = c(0), verbose = FALSE)
  expect_is(ESDM, 'Ensemble.SDM')
})
