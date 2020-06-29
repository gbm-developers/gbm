#### Xpcol ####
test_that('Xcol argument checking works',{
  expect_error(.checkargs(Xcol = 1))
})

#### Ypcol ####
test_that('Ycol argument checking works',{
  expect_error(.checkargs(Ycol = NA))
})

#### Ppcol ####
test_that('Pcol argument checking works',{
  expect_error(.checkargs(Pcol = 1))
})

#### Spcol ####
test_that('Spcol argument checking works',{
  expect_error(.checkargs(Spcol = 2))
})

#### name ####
test_that('name argument checking works',{
  expect_error(.checkargs(name = 1))
})

#### Spname ####
test_that('Spname argument checking works',{
  expect_error(.checkargs(Spname = 1))
})

#### PA ####
test_that('PA argument checking works',{
  expect_error(.checkargs(PA = 1))
})

#### rep ####
test_that('rep argument checking works',{
  expect_error(.checkargs(rep = 'deux'))
})

#### cv ####
test_that('cv argument checking works',{
  expect_error(.checkargs(cv = 4))
})

#### cv.param ####
test_that('cv.param argument checking works',{
  expect_error(.checkargs(cv.param = 1))
})

#### thresh ####
test_that('thresh argument checking works',{
  expect_error(.checkargs(thresh = 'thousand'))
})

#### metric ####
test_that('metric argument checking works',{
  expect_error(.checkargs(metric = 2))
})

#### axes.metric ####
test_that('axes.metric argument checking works',{
  expect_error(.checkargs(axes.metric = 2))
})

#### select ####
test_that('select argument checking works',{
  expect_error(.checkargs(select = 'TRUE'))
})

#### select.metric ####
test_that('select.metric argument checking works',{
  expect_error(.checkargs(select.metric = 2))
})

#### select.thresh ####
test_that('select.thresh argument checking works',{
  expect_error(.checkargs(select.thresh = 'zero'))
})

#### verbose ####
test_that('verbose argument checking works',{
  expect_error(.checkargs(verbose = 'TRUE'))
})

#### GUI ####
test_that('GUI argument checking works',{
  expect_error(.checkargs(GUI = 'TRUE'))
})

#### uncertainty ####
test_that('uncertainty argument checking works',{
  expect_error(.checkargs(uncertainty = 'TRUE'))
})

#### tmp ####
test_that('tmp argument checking works',{
  expect_error(.checkargs(tmp = 'TRUE'))
})

#### ensemble.metric ####
test_that('ensemble.metric argument checking works',{
  expect_error(.checkargs(ensemble.metric = 2))
})

#### ensemble.thresh ####
test_that('ensemble.thresh argument checking works',{
  expect_error(.checkargs(ensemble.thresh = 'zero'))
})

#### weight ####
test_that('weight argument checking works',{
  expect_error(.checkargs(weight = 'TRUE'))
})

#### method ####
test_that('method argument checking works',{
  expect_error(.checkargs(method = 2))
})

#### rep.B ####
test_that('rep.B argument checking works',{
  expect_error(.checkargs(method = 'B', rep.B = 'thousand'))
})

#### GeoRes ####
test_that('GeoRes argument checking works',{
  expect_error(.checkargs(GeoRes = 'TRUE'))
})

#### reso ####
test_that('reso argument checking works',{
  expect_error(.checkargs(reso = 'zero'))
})

#### file ####
test_that('file argument checking works',{
  expect_error(.checkargs(file = 1))
})

#### files ####
test_that('files argument checking works',{
  expect_error(.checkargs(files = 1))
})

#### format ####
test_that('format argument checking works',{
  expect_error(.checkargs(format = 2))
})

#### categorical ####
test_that('categorical argument checking works',{
  expect_error(.checkargs(categorical = 1))
})

#### Norm ####
test_that('Norm argument checking works',{
  expect_error(.checkargs(Norm = 'TRUE'))
})

#### enm ####
test_that('enm argument checking works',{
  expect_error(.checkargs(enm = numeric()))
})

#### stack ####
test_that('stack argument checking works',{
  expect_error(.checkargs(stack = numeric()))
})

#### range ####
test_that('range argument checking works',{
  expect_error(.checkargs(range = 'one-two'))
})

#### endemism ####
test_that('endemism argument checking works',{
  expect_error(.checkargs(endemism = c(1,2)))
})

#### cores ####
test_that('cores argument checking works',{
  expect_error(.checkargs(cores = 'one'))
})

