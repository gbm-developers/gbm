test_that('project Algorithm.SDM', {
  data(Env)
  data(Occurrences)
  Occurrences <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  SDM <- modelling('GLM', Occurrences, Env,
                   Xcol = 'LONGITUDE', Ycol = 'LATITUDE', verbose = F)
  Env_new <- stack(Env[[1]]-1,Env[[2]],Env[[3]])
  SDM_proj <- project(SDM,Env_new)
  expect_is(SDM_proj,'GLM.SDM')
  expect_false(all(is.na(values(SDM_proj@projection))))
})
