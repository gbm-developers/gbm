test_that("modelling works", {
  model <- Algorithm.SDM("GLM")
  model@name <- paste0("test", "GLM", ".SDM")
  model@parameters$data <- "presence-only data set"
  model@parameters$metric <- "SES"
  data("Occurrences")
  data <- data.frame(X = Occurrences$LONGITUDE, Y = Occurrences$LATITUDE)
  data$Presence <- 1
  model@data <- data
  data(Env)
  model = PA.select(model, Env, NULL, verbose = FALSE)
  model@parameters["PA"] = TRUE
  expect_equal(length(which(model@data$Presence == 0)), 1000)
  model = data.values(model, Env)
  expect_equal(dim(model@data), c(1057, 6))
  model = evaluate(model, cv = "holdout", cv.param = c(0.7, 2), thresh = 1001,
                   metric = "SES", Env)
  expect_equal(dim(model@evaluation), c(1, 7))
  model = project(model, Env)
  expect_equal(all(is.na(values(model@projection))), FALSE)
  model = evaluate.axes(model, cv.param = c(0.7, 2), thresh = 1001, metric = "SES",
                        axes.metric = "Pearson", Env)
  expect_equal(dim(model@variable.importance), c(1, 3))
  # expect_equal(NULL, print(model))
})
