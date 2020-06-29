skip_on_cran()
skip_on_appveyor()

# .get_fc_args
test_that("The correct fc arguments are created", {
  expect_equal(.get_fc_args("lqpht"),
               c("noautofeature", "threshold"))
  expect_error(.get_fc_args("lb"))
})

# .make_args
test_that("The correct arguments are created", {
  expect_equal(.make_args(1, "l", 500, "removeduplicates=false"),
               c("betamultiplier=1", "maximumiterations=500", "noautofeature",
                 "noquadratic", "noproduct", "nohinge",
                 "removeduplicates=false"))
})

# .get_lambdas
const <- c("linearPredictorNormalizer, 1.68",
           "densityNormalizer, 21.50",
           "numBackgroundPoints, 5000",
           "entropy, 5.80")

linear <- c("bio1, -0.44, 61.0, 42.0", const)
quadratic <- c("bio1^2, -0.44, 61.0, 42.0", const)
product <- c("bio1*bio16, -0.44, 61.0, 42.0", const)
hinge <- c("'bio1, -0.44, 61.0, 42.0", const)
rev_hinge <- c("`bio1, -0.44, 61.0, 42.0", const)
threshold <- c("(500<bio1), -0.44, 0.0, 1.0", const)
categ <- c("(biome=1.0), -0.44, 0.0, 1.0", const)

test_that("Constants are parsed correctly", {
  l <- .get_lambdas(linear)
  # Get the correct linear predictor normalizer
  expect_equal(l$lpn, 1.68)
  # Get the correct density normalizer
  expect_equal(l$dn, 21.50)
  # Get the correct entropy
  expect_equal(l$entropy, 5.80)
  # Lambdas has correct column names
  expect_named(l$lambdas, c("feature", "lambda", "min", "max"))
  # Min max has the correct column names
  expect_named(l$min_max, c("variable", "min", "max"))
  # Min max has the correct length
  expect_length(l$min_max$variable, 1)
  # Output has the correct column names
  expect_named(l, c("lambdas", "lpn", "dn", "entropy", "min_max"))
  l <- .get_lambdas(c("bio1, 0.0, 61.0, 42.0", linear))
  # Remove feature in which lambda is zero
  expect_length(l$lambdas$feature, 1)
})

test_that("The correct linear feature is created", {
  l <- .get_lambdas(linear)
  expect_equal(l$lambdas,
               data.frame(feature = "bio1", lambda = -0.44, min = 61.0,
                          max = 42.0, stringsAsFactors = FALSE))
  expect_equal(l$min_max,
               data.frame(variable = "bio1", min = 61.0, max = 42.0,
                          stringsAsFactors = FALSE))
})

test_that("The correct quadratic feature is created", {
  expect_equal(.get_lambdas(quadratic)$lambdas,
               data.frame(feature = "I(bio1^2)", lambda = -0.44, min = 61.0,
                          max = 42.0, stringsAsFactors = FALSE))
})

test_that("The correct product feature is created", {
  expect_equal(.get_lambdas(product)$lambdas,
               data.frame(feature = "I(bio1*bio16)", lambda = -0.44, min = 61.0,
                          max = 42.0, stringsAsFactors = FALSE))
})

test_that("The correct hinge feature is created", {
  expect_equal(.get_lambdas(hinge)$lambdas,
               data.frame(feature = ".hinge(bio1, 61, 42)", lambda = -0.44,
                          min = 61.0, max = 42.0, stringsAsFactors = FALSE))
  # Variable lower than minimum (hinge)
  expect_equal(.hinge(41.9, 42, 61), 0)
  # Variable equal to minimum (hinge)
  expect_equal(.hinge(42, 42, 61), 0)
  # Variable greater than minimum
  expect_equal(.hinge(350, 300, 400), 0.5)
  # Variable equal to maximum
  expect_equal(.hinge(61, 42, 61), 1)
})

test_that("The correct reverse hinge feature is created", {
  expect_equal(.get_lambdas(rev_hinge)$lambdas,
               data.frame(feature = ".rev_hinge(bio1, 61, 42)", lambda = -0.44,
                          min = 61.0, max = 42.0, stringsAsFactors = FALSE))
  # Variable equal to minimum
  expect_equal(.rev_hinge(42, 42, 61), 1)
  # Variable lower than maximum
  expect_equal(.rev_hinge(350, 300, 400), 0.5)
  # Variable equal to maximum (hinge)
  expect_equal(.rev_hinge(61, 42, 61), 0)
  # Variable greater than maximum (hinge)
  expect_equal(.rev_hinge(63, 42, 61), 0)
})

test_that("The correct threshold feature is created", {
  expect_equal(.get_lambdas(threshold)$lambdas,
               data.frame(feature = ".threshold(bio1, 500)", lambda = -0.44,
                          min = 0, max = 1, stringsAsFactors = FALSE))
  # Variable lower than threshold
  expect_equal(.threshold(300, 500), 0)
  # Variable equal to threshold
  expect_equal(.threshold(500, 500), 1)
  # Variable greater than threshold
  expect_equal(.threshold(600, 500), 1)
})

test_that("The correct categorical feature is created", {
  expect_equal(.get_lambdas(categ)$lambdas,
               data.frame(feature = ".categorical(biome, 1.0)", lambda = -0.44,
                          min = 0, max = 1, stringsAsFactors = FALSE))
  # Variable lower than category
  expect_equal(.categorical(0, 1), 0)
  # Variable equal to category
  expect_equal(.categorical(1, 1), 1)
  # Variable greater than category
  expect_equal(.categorical(2, 1), 0)
})

skip_on_travis()
skip_on_appveyor()
skip_on_covr()

#.train
test_that("The function trainMaxent produces the correct ouput", {
  m <- trainMaxent(data = SDMtune:::t, reg = 1.2, fc = "l")
  expect_s4_class(m, "SDMmodel")
  expect_s4_class(m@model, "Maxent")
  expect_s4_class(m@data, "SWD")
  expect_equal(m@model@reg, 1.2)
  expect_equal(m@model@fc, "l")
  expect_equal(m@data, SDMtune:::t)
})
