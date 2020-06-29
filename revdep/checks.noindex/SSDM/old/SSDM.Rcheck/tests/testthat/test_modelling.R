test_that('modelling function', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  SDM <- modelling('GLM', Occurrences, Env,
                   Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = F)
  expect_is(SDM, 'GLM.SDM')
})
