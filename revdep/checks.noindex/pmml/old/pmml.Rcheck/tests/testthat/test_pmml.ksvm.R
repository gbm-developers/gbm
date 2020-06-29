library(kernlab)
data(iris)

teardown({
  detach("package:kernlab", unload = TRUE)
})

test_that("pmml.ksvm error when a specified kernel is not supported", {
  expect_error(
    pmml(ksvm(Sepal.Length ~ ., data = iris, kernel = "laplacedot"), dataset = iris),
    "laplacedot kernel is not supported. Supported ksvm kernels: rbfdot, polydot, vanilladot, tanhdot."
  )

  invisible(capture.output(mod2 <- ksvm(Sepal.Length ~ ., data = iris, kernel = "besseldot")))
  expect_error(
    pmml(mod2, dataset = iris),
    "besseldot kernel is not supported. Supported ksvm kernels: rbfdot, polydot, vanilladot, tanhdot."
  )

  invisible(capture.output(mod3 <- ksvm(Sepal.Length ~ ., data = iris, kernel = "anovadot")))
  expect_error(
    pmml(mod3, dataset = iris),
    "anovadot kernel is not supported. Supported ksvm kernels: rbfdot, polydot, vanilladot, tanhdot."
  )

  invisible(capture.output(mod4 <- ksvm(Sepal.Length ~ ., data = iris, kernel = "splinedot")))
  expect_error(
    pmml(mod4, dataset = iris),
    "splinedot kernel is not supported. Supported ksvm kernels: rbfdot, polydot, vanilladot, tanhdot."
  )

  # The following test does not run due to a ksvm() error - "Error in .local(x, ...) :  kernel matrix not square!"
  # expect_error(pmml(ksvm(Sepal.Length ~ ., data=iris,kernel="matrix"),dataset=iris),
  #            "matrix kernel is not supported. Supported ksvm kernels: rbfdot, polydot, vanilladot, tanhdot.")
})
