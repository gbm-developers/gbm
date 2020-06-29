context("Check aspect_importance_single() functions")

test_that("check aspect_importance_single function",{
  library("DALEX")
  library("triplot")

  aspect_importance_titanic_single <-
    aspect_importance_single(x = apartments_lm_model,
                             data = 
                               apartments[, colnames(apartments) != "m2.price"],
                             new_observation = apartments_new_observation)
  
  expect_true("data.frame" %in% class(aspect_importance_titanic_single))
  expect_true(dim(aspect_importance_titanic_single)[1] == 5)
  expect_true(dim(aspect_importance_titanic_single)[2] == 5)

})

test_that("check aspect_importance_single.explainer function",{
  library("DALEX")
  library("triplot")

  titanic_without_target <- titanic_data[,colnames(titanic_data)!="survived"]

  titanic_explainer <- explain(model = titanic_glm_model,
                               data = titanic_without_target,
                               verbose = FALSE)

  aspect_importance_titanic_glm_single <-
    aspect_importance_single(titanic_explainer,
                             new_observation = titanic_new_observation)
  
  expect_true("data.frame" %in% class(aspect_importance_titanic_glm_single))

})
