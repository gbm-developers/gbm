context("Check triplot() functions")

test_that("check calculate_triplot.default function",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_num_lm_model,
                                      data = apartments_num[,-1],
                                      y = apartments_num[, 1],
                                      new_observation =
                                        apartments_num_new_observation[-1])
  expect_true("list" %in% class(apartments_tri))
})



test_that("check calculate_triplot.explainer function",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      new_observation =
                                        apartments_num_new_observation[-1])

  expect_true("list" %in% class(apartments_tri))
})

test_that("check print of triplot",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      new_observation =
                                        apartments_num_new_observation[-1])
  apartments_tri_fi <- calculate_triplot(x = apartments_num_lm_model,
                                      data = apartments_num[,-1],
                                      y = apartments_num[, 1],
                                      type = "model")
  expect_output(print(apartments_tri), "Triplot object")
  expect_output(print(apartments_tri_fi), "Triplot object")
})

test_that("check warning in calculate_triplot.explainer",{
  library("DALEX")
  library("triplot")

  apartments_num_explain_2 <- explain(model = apartments_num_lm_model,
                                      data = apartments_num,
                                      y = apartments_num[, 1],
                                      verbose = FALSE)

  expect_warning(calculate_triplot(apartments_num_explain_2,
                                   new_observation =
                                     apartments_num_new_observation))

})

test_that("check calculate_triplot.default function for FI",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_num_lm_model,
                                      data = apartments_num[,-1],
                                      y = apartments_num[, 1],
                                      type = "model")
  expect_true("list" %in% class(apartments_tri))
})

test_that("check calculate_triplot.explainer function for FI",{
  library("DALEX")
  library("triplot")
  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      type = "model")

  expect_true("list" %in% class(apartments_tri))
})


test_that("check plot.calculate_triplot function",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      new_observation =
                                        apartments_num_new_observation[-1])
  p <- plot(apartments_tri)
  
  p2 <- plot(apartments_tri,
             abbrev_labels = 5,
             text_size = 4,
             absolute_value = TRUE,
             add_importance_labels = TRUE,
             show_model_label = TRUE,
             add_last_group = FALSE,
             axis_lab_size = 4,
             bar_width = 3,
             margin_mid = 0.1)

  expect_true("patchwork" %in% class(p))
  expect_true("gg" %in% class(p))
  expect_true("ggplot" %in% class(p2))
})

test_that("check plot.calculate_triplot function for FI",{
  library("DALEX")
  library("triplot")

  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      type = "model")
  p <- plot(apartments_tri, 
            abbrev_labels = 5)
  
  p2 <- plot(apartments_tri,
            add_last_group = FALSE)


  expect_true("patchwork" %in% class(p))
  expect_true("gg" %in% class(p))
  expect_true("ggplot" %in% class(p2))
  
})

test_that("",{
  library("DALEX")
  library("triplot")
  
  apartments_tri <- calculate_triplot(x = apartments_explain,
                                      type = "model")

  expect_true(dim(apartments_tri[[1]])[1] == 66)
  expect_true(class(apartments_tri[[2]]) == "hierarchical_importance")
  expect_true("hclust" %in% class(apartments_tri[[3]]))
  expect_true(apartments_tri[[5]] == "model")
  
})


test_that("check triplot aliases",{
  library("DALEX")
  library("triplot")

  apartments_tri_model <- model_triplot(x = apartments_explain)

  apartments_tri_predict <- predict_triplot(x = apartments_explain,
                                    new_observation =
                                      apartments_num_new_observation[-1])

  expect_true("triplot" %in% class(apartments_tri_model))
  expect_true("triplot" %in% class(apartments_tri_predict))

})

test_that("check for triplot error",{
  library("DALEX")
  library("triplot")

  expect_error(predict_triplot(x = apartments_explain))


})

