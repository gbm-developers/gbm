context('checks working properly')

set.seed(2020)

expected_fi <- model_fi[model_fi$permutation == 0,
                        c('variable', 'dropout_loss')]
class(expected_fi) <- 'data.frame'
colnames(expected_fi) <- c('label', 'value')
model_pd_list <- list(numerical = model_pd)
pd_opts <- list(N = 100, grid_points = 81)

test_that('check_feature_importance working', {
  expect_error(check_feature_importance(feature_importance = 'ABC'),
               regexp = 'feature_importance')
  expect_silent(check_feature_importance(feature_importance = model_fi))
  expect_silent(check_feature_importance(feature_importance = list(loss_function = DALEX::loss_accuracy,
                                                                   type = 'raw')))
})

test_that('process_feature_importance working', {
  expect_equal(process_feature_importance(model_fi, model_exp),
               expected_fi,
               check.attributes = FALSE)
  expect_equal(process_feature_importance(
    list(loss_function = DALEX::loss_accuracy,
         type = 'raw'),
    model_exp
  ),
  expected_fi,
  check.attributes = FALSE)
})
test_that('check_partial_dependence working for numerical data', {
  expect_error(check_partial_dependence(partial_dependence = 'ABC',
                                        model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = list('ABC'),
                                        model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = model_pd,
                                        model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = list(model_pd),
                                        model_exp),
               regexp = 'partial_dependence')
  expect_silent(check_partial_dependence(partial_dependence = list(numerical = model_pd),
                                         model_exp))
  expect_silent(check_partial_dependence(partial_dependence = list(numerical = pd_opts),
                                         model_exp))
})  
test_that('check_partial_dependence working for mixed data', {
  expect_error(check_partial_dependence(partial_dependence = 'ABC',
                                        model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = list('ABC'),
                                        model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = list(numerical = 'ABC'),
                                        tit_model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = list(categorical = 'ABC'),
                                        tit_model_exp),
               regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = {
    l <- tit_model_pd
    l$categorical <- NULL
    l
  },
  tit_model_exp),
  regexp = 'partial_dependence')
  expect_error(check_partial_dependence(partial_dependence = {
    l <- tit_model_pd
    l$numerical <- NULL
    l
  },
  tit_model_exp),
  regexp = 'partial_dependence')
  expect_silent(check_partial_dependence(partial_dependence = tit_model_pd,
                                         tit_model_exp))
  expect_silent(check_partial_dependence(partial_dependence = list(numerical = tit_model_pd$numerical,
                                                                   categorical = pd_opts),
                                         tit_model_exp))
  expect_silent(check_partial_dependence(partial_dependence = list(categorical = tit_model_pd$categorical,
                                                                   numerical = pd_opts),
                                         tit_model_exp))
  expect_silent(check_partial_dependence(partial_dependence = list(numerical = pd_opts,
                                                                   categorical = pd_opts),
                                         tit_model_exp))
})

test_that('process_partial_dependence working', {
  expect_equal(process_partial_dependence(model_pd_list, model_exp),
               model_pd_list,
               check.attributes = FALSE)
  expect_equal(process_partial_dependence(list(numerical = pd_opts), 
                                          model_exp),
  model_pd_list,
  tolerance = 0.1,
  check.attributes = FALSE)
  expect_equal(process_partial_dependence(tit_model_pd, tit_model_exp),
               tit_model_pd,
               tolerance = 0.15,
               check.attributes = FALSE)
  expect_equal(process_partial_dependence(list(numerical = tit_model_pd$numerical,
                                               categorical = pd_opts), tit_model_exp),
               tit_model_pd,
               tolerance = 0.15,
               check.attributes = FALSE)
  expect_equal(process_partial_dependence(list(numerical = pd_opts,
                                               categorical = tit_model_pd$categorical), tit_model_exp),
               tit_model_pd,
               tolerance = 0.15,
               check.attributes = FALSE)
  expect_equal(process_partial_dependence(list(numerical = pd_opts,
                                               categorical = pd_opts), tit_model_exp),
               tit_model_pd,
               tolerance = 0.15,
               check.attributes = FALSE)
})



