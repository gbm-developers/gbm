context("Check hierarchical_importance() functions")

test_that("check hierarchical_importance function for aspects_importance",{
  library("DALEX")
  library("triplot")

  hi <- hierarchical_importance(x = apartments_num_lm_model,
                                data = apartments_num[,-1],
                                new_observation =
                                  apartments_num_new_observation[-1])

  expect_true("hclust" %in% class(hi[[1]]))
  expect_true("hierarchical_importance" %in% class(hi))
  expect_true("floor" %in% hi[[1]]$labels)
})

test_that("check hierarchical_importance function for feature_importance",{
  library("DALEX")
  library("triplot")

  hi <- hierarchical_importance(x = apartments_num_lm_model,
                                data = apartments_num[,-1],
                                y = apartments_num[,1],
                                type = "model")

  expect_true("hierarchical_importance" %in% class(hi))
  expect_true("floor" %in% hi[[1]]$labels)
})


test_that("check plot.hierarchical_importance  function",{
  library("DALEX")
  library("triplot")

  hi <- hierarchical_importance(x = apartments_num_lm_model,
                                data = apartments_num[,-1],
                                new_observation =
                                  apartments_num_new_observation[-1])

  p <- plot(hi, add_last_group = TRUE,
            absolute_value = TRUE)

  expect_true("ggplot" %in% class(p))
})

test_that("check for hierarchical_importance error",{
  library("DALEX")
  library("triplot")

  expect_error(hierarchical_importance(x = apartments_num_lm_model,
                                         data = apartments_num[,-1],
                                         type = "model"))
})
