context('calculate_cors working properly')

data("Seatbelts")
Seatbelts <- as.data.frame(Seatbelts)[-8]
data(mtcars)
data(titanic_imputed, package = 'DALEX')

titanic_imputed <- titanic_imputed[,-8]

mtcars$vs <- factor(mtcars$vs, labels = c('V-shaped', 'straight'))
mtcars$am <- factor(mtcars$am, labels = c('automatic', 'manual'))

data("HairEyeColor")

test_that('calculate_cors working properly for numeric data',{
  expect_is({cors <- calculate_cors(Seatbelts)
  cors}, 'matrix')
  expect_equal(ncol(cors), ncol(Seatbelts))
  expect_equal(nrow(cors), ncol(Seatbelts))
  expect_true(isSymmetric(cors))
  expect_true(all(abs(cors) <= 1))
  expect_equal(colnames(cors), colnames(Seatbelts))
  expect_equal(rownames(cors), colnames(Seatbelts))
})

test_that('calculate_cors working properly for mixed data with binary cats',{
  expect_is({cors <- calculate_cors(mtcars)
  cors}, 'matrix')
  expect_equal(ncol(cors), ncol(mtcars))
  expect_equal(nrow(cors), ncol(mtcars))
  expect_true(isSymmetric(cors))
  expect_equal(colnames(cors), colnames(mtcars))
  expect_equal(rownames(cors), colnames(mtcars))
})

test_that('calculate_cors working properly for mixed data with non-binary cats',{
  expect_is({cors <- calculate_cors(titanic_imputed)
  cors}, 'matrix')
  expect_equal(ncol(cors), ncol(titanic_imputed))
  expect_equal(nrow(cors), ncol(titanic_imputed))
  expect_true(isSymmetric(cors))
  expect_equal(colnames(cors), colnames(titanic_imputed))
  expect_equal(rownames(cors), colnames(titanic_imputed))
})

test_that('calculate_cors working properly for categorical data in contingency table',{
  expect_is({cors <- calculate_cors(HairEyeColor)
  cors}, 'matrix')
  expect_equal(ncol(cors), length(dim(HairEyeColor)))
  expect_equal(nrow(cors), length(dim(HairEyeColor)))
  expect_true(isSymmetric(cors))
  expect_equal(colnames(cors), names(dimnames(HairEyeColor)))
  expect_equal(rownames(cors), names(dimnames(HairEyeColor)))
})

test_that('numeric data requires only num_num_f',{
  expect_equal(calculate_cors(Seatbelts,
                              num_num_f = cor,
                              max_cor = 1),
               cor(Seatbelts))
  expect_error(calculate_cors(Seatbelts,
                              num_cat_f = wilcox.test))
})

test_that('mixed data requires all functions',{
  expect_error(calculate_cors(mtcars,
                              num_num_f = cor))
  expect_error(calculate_cors(mtcars,
                              num_num_f = cor,
                              num_cat_f = wilcox.test))
  expect_error(calculate_cors(mtcars,
                              num_num_f = cor,
                              cat_cat_f = chisq.test))
  expect_error(calculate_cors(mtcars,
                              num_cat_f = wilcox.test,
                              cat_cat_f = chisq.test))
})

test_that('mixed data with 1 cat column requires num_cat',{
  expect_error(calculate_cors(mtcars[,c('mpg', 'cyl', 'am')],
                              num_num_f = cor))
})

test_that('data with 1 numeric and 1 cat requires only num_cat',{
  expect_equal(calculate_cors(mtcars[,c('qsec', 'vs')],
                              num_cat = function(x,y) -log10(kruskal.test(x, y)[['p.value']]),
                              max_cor = 10),
               {
                 val <- -log10(kruskal.test(mtcars$qsec, mtcars$vs)$p.value) / 10
                 d <- matrix(c(1, val, val, 1), ncol = 2)
                 colnames(d) <- c('qsec', 'vs')
                 rownames(d) <- c('qsec', 'vs')
                 d
               }
               )
})
