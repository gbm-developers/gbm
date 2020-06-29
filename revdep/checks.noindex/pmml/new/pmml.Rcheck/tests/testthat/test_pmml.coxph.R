library(survival)
# most of the tests are from survival::coxph doc

teardown({
  detach("package:survival", unload = TRUE)
})

test_that("error when object is not coxph", {
  expect_error(pmml.coxph("foo"), "Not a legitimate coxph object")
})

test_that("no error for stratified model", {
  # stratified model
  test1 <- list(
    time = c(4, 3, 1, 1, 2, 2, 3),
    status = c(1, 1, 1, 0, 1, 1, 0),
    x = c(0, 2, 1, 1, 1, 0, 0),
    sex = c(0, 0, 0, 0, 1, 1, 1)
  )
  fit1 <- coxph(Surv(time, status) ~ x + strata(sex), test1)
  expect_error(pmml(fit1), NA) # expect no error

  # time-dependent model
  test2 <- list(
    start = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
    stop = c(2, 3, 6, 7, 8, 9, 9, 9, 14, 17),
    event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
    x = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
  )
  fit2 <- coxph(Surv(start, stop, event) ~ x, test2)
  expect_error(pmml(fit2), NA) # expect no error
})

test_that("Error for multiplicative strata variables", {
  # stratified model clustered on patients
  data(bladder)
  bladder1 <- bladder[bladder$enum < 5, ]
  fit3 <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) +
    cluster(id), bladder1)
  expect_error(pmml(fit3), "Multiplicative strata variables not yet supported in PMML")
})

test_that("Error for model with time-transform", {
  # time transform model using current age
  fit4 <- coxph(Surv(time, status) ~ ph.ecog + tt(age),
    data = lung,
    tt = function(x, t, ...) pspline(x + t / 365.25)
  )
  expect_error(
    pmml(fit4),
    "Special model equation terms 'cluster' and 'tt' not yet supported in PMML"
  )
})
