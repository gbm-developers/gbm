context('corrgrapher working properly for explainers')

options(check.attributes = FALSE)
set.seed(2020)

custom_values <- data.frame(label = colnames(dragons)[-5],
                            value = rep(15, ncol(dragons) - 1),
                            stringsAsFactors = FALSE)
custom_values <- custom_values[order(custom_values$label),]
rownames(custom_values) <- NULL

model_pd_list <- list(numerical = model_pd)

test_that(
  'Function is working properly with just necessary arguments',{
    expect_is(corrgrapher(simple_model_exp),'corrgrapher')
    expect_is(corrgrapher(tit_model_exp),'corrgrapher')
  }
)

test_that('Values argument working', {
  expect_equal({
    df <- corrgrapher(model_exp,
                      values = custom_values,
                      partial_dependency = model_pd_list)[['nodes']][, c('label', 'value')]
    df[order(df$label),]
  },
  custom_values)})

test_that('Values argument overrides feature_importance_*',{
  expect_warning(
    cgr <- corrgrapher(
      model_exp,
      values = custom_values,
      feature_importance = model_fi,
      partial_dependency = model_pd_list
    )
  )
  expect_equal({
    df <- cgr[['nodes']][, c('label', 'value')]
    df <- df[order(df$label), ]
    rownames(df) <- NULL
    df
  },
  custom_values)
  expect_warning(
    cgr <- corrgrapher(
      model_exp,
      values = custom_values,
      feature_importance = list(),
      partial_dependency = model_pd_list
    )
  )
  expect_equal({
    df <- cgr[['nodes']][,c('label','value')]
    df <- df[order(df$label),]
    rownames(df) <- NULL
    df
  },
  custom_values)
})
  
test_that("Output type",{
  expect_is(cgr_exp, 'corrgrapher')
  expect_true(all(c("nodes", "edges", "pds") %in% names(cgr_exp)))
})


