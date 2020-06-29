context("Check group_variables() functions")

test_that("check cluster_variables function",{
  library("DALEX")
  library("triplot")

  cv <- cluster_variables(apartments_num, clust_method = "single")

  expect_true(length(cv$labels) == 5)
  expect_true("cluster_variables" %in% class(cv))
  expect_true("hclust" %in% class(cv))

})

test_that("check plot.cluster_variables function",{
  library("DALEX")
  library("triplot")

  cv <- cluster_variables(apartments_num, clust_method = "single")
  p <- plot(cv, p = 0.5)

  expect_true("ggplot" %in% class(p))
})

test_that("check list_variables function",{
  library("DALEX")
  library("triplot")

  cv <- cluster_variables(apartments_num, clust_method = "single")
  aspect_list <- list_variables(cv, 0.6)
  one_aspect <- list_variables(cv, 0)

  expect_true(class(aspect_list) == "list")
  expect_true(length(aspect_list) == 4)
  expect_true(length(one_aspect) == 1)
  expect_true("surface" %in% aspect_list$aspect.group3)

})

test_that("check group_variables function",{
  library("DALEX")
  library("triplot")

  aspect_list <- group_variables(x = apartments_num, h = 0.6,
                                 clust_method = "single")

  expect_true(class(aspect_list) == "list")

})
