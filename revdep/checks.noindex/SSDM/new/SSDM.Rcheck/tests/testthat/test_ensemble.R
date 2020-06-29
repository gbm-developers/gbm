test_that('ensemble function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  SDM1 <- modelling('GLM', Occurrences, Env, Xcol = 'LONGITUDE',
                    Ycol = 'LATITUDE', select.thresh = 0, verbose = FALSE)
  SDM2 <- modelling('SVM', Occurrences, Env, Xcol = 'LONGITUDE',
                    Ycol = 'LATITUDE', select.thresh = 0, verbose = FALSE)
  ESDM <- ensemble(SDM1, SDM2, ensemble.thresh = 0, verbose = FALSE)
  expect_is(ESDM, 'Ensemble.SDM')
})
