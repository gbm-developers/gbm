context('corrgrapher working properly')

df <- as.data.frame(datasets::Seatbelts)[,1:7]
data('freeny')

test_that("Output type",{
  expect_is({
    cgr <- corrgrapher(df)
    cgr
    }, 'corrgrapher')
  expect_true(all(c("nodes", "edges") %in% names(cgr)))
  expect_is({
    cgr <- corrgrapher(freeny.x)
    cgr
  }, 'corrgrapher')
  expect_true(all(c("nodes", "edges") %in% names(cgr)))
})

test_that('Incorrect argument x caught',{
  expect_error(corrgrapher(1:5))
})

test_that('Incorrect argument cutoff caught',{
  expect_error(corrgrapher(df, cutoff = 'ABC'))
  expect_error(corrgrapher(df, cutoff = 1:5))
})

test_that('Unusual argument cutoff recognised', {
  expect_warning(corrgrapher(df, cutoff = 2))
  expect_warning(corrgrapher(df, cutoff = -1))
})

test_that('Incorrect value argument caught', {
  expect_error(corrgrapher(df, values = 'ABC'))
  expect_error(corrgrapher(df, values = data.frame(A = 1:5, B = 1:5)))
  expect_error(corrgrapher(df, values = data.frame(label = c('A', 'B', 'C'),
                                                   value = 1:3)))
  expect_error(corrgrapher(df, values = data.frame(label = colnames(df),
                                                  values = LETTERS[1:ncol(df)])))
})