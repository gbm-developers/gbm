skip_on_cran()
skip_on_appveyor()

data <- SDMtune:::t
mother <- SDMtune:::bm_maxnet
father <- train("Maxnet", data = data, fc = "l", reg = 2)
h <- list(fc = c("l", "lq", "lqph"), reg = c(1, 2))
metrics <- list(c(10, 11, 12), c(8, 10, 13))

test_that("Exception are raised", {
  expect_error(optimizeModel(mother, h, "auc", data, keep_best = 0.6,
                             keep_random = 0.6, pop = 3),
               "Sum of 'keep_best' and 'keep_random' cannot be more than 1!")
  expect_error(optimizeModel(mother, h, "auc", data, pop = 3),
               "Optimization algorithm interrupted at generation 0 because it overfits validation dataset!")
})

# TODO Rewrite this test
# test_that("The algorithm executes without errors", {
#   expect_s4_class(optimizeModel(mother, h, "auc", test = data, pop = 3, gen = 1),
#                   "SDMtune")
# })

test_that("Crossover is executed", {
  set.seed(30, kind = "Mersenne-Twister", sample.kind = "Rejection")
  x <- .breed(mother, father, h, mutation_chance = 0)
  # fc comes from father
  expect_equal(x@model@fc, father@model@fc)
  # reg comes from mother
  expect_equal(x@model@reg, mother@model@reg)
})

test_that("Mutation is executed", {
  # For an hyperparameter different from a
  set.seed(25, kind = "Mersenne-Twister", sample.kind = "Rejection")
  x <- .breed(mother, father, h, mutation_chance = 1)
  # fc comes from mutation
  expect_equal(x@model@fc, "lq")
})

test_that("The rank is correct", {
  # For AICc the most important is the one with the lowest metric
  expect_equal(.get_rank_index("aicc", metrics), c(1, 2, 3))
  # For AUC or TSS the most important is the one with the highest value not
  # overfitting
  expect_equal(.get_rank_index("auc", metrics), c(2, 1, 3))
  # All model are overfitting
  metrics <- list(c(10, 11, 12), c(11, 12, 13))
  expect_false(.get_rank_index("auc", metrics))
})
