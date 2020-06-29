test_that("function_to_pmml('1+2') outputs correct xml", {
  current <- function_to_pmml("1 + 2")

  node <- newXMLNode(name = "Apply", attrs = c("function" = "+"))
  c1 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "1")
  c2 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "2")
  expected <- addChildren(node, kids = c(c1, c2))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})

test_that("function_to_pmml('foo(bar(baz))') outputs correct xml", {
  current <- function_to_pmml("foo(bar(baz))")

  node <- newXMLNode(name = "Apply", attrs = c("function" = "foo"))
  c1 <- newXMLNode(name = "Apply", attrs = c("function" = "bar"))
  c2 <- newXMLNode(name = "FieldRef", attrs = c("field" = "baz"))
  expected <- addChildren(node, addChildren(c1, c2))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})

test_that("function_to_pmml('1(2)') throws unexpected end of input error", {
  expect_error(function_to_pmml("1(2"), regexp = "unexpected end of input")
})


test_that("function_to_pmml('-3') outputs correct xml", {
  current <- function_to_pmml("-3")

  node <- newXMLNode(name = "Apply", attrs = c("function" = "-"))
  c1 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "0")
  c2 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "3")
  expected <- addChildren(node, kids = c(c1, c2))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})


test_that("function_to_pmml('-(44*a)') outputs correct xml", {
  current <- function_to_pmml("-(44*a)")

  node <- newXMLNode(name = "Apply", attrs = c("function" = "-"))
  c1 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "0")
  c1node <- newXMLNode(name = "Apply", attrs = c("function" = "*"))
  c2 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "44")
  c3 <- newXMLNode(name = "FieldRef", attrs = c("field" = "a"))
  addChildren(c1node, kids = c(c2, c3))

  expected <- addChildren(node, kids = c(c1, c1node))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})

test_that("function_to_pmml('-a') outputs correct xml", {
  current <- function_to_pmml("-a")

  node <- newXMLNode(name = "Apply", attrs = c("function" = "-"))
  c1 <- newXMLNode(name = "Constant", attrs = c("dataType" = "double"), text = "0")
  c2 <- newXMLNode(name = "FieldRef", attrs = c("field" = "a"))
  expected <- addChildren(node, kids = c(c1, c2))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})


test_that("function_to_pmml('?3') throws error when ? is * or /", {
  expect_error(function_to_pmml("*3"), regexp = "<text>:1:1: unexpected '*'")
  expect_error(function_to_pmml("/3"), regexp = "<text>:1:1: unexpected '/'")
})

test_that("function_to_pmml outputs boolean TRUE/FALSE for if function", {
  current <- function_to_pmml("if(out < t){TRUE} else {FALSE}")

  c0node <- newXMLNode(name = "Apply", attrs = c("function" = "if"))
  c1node <- newXMLNode(name = "Apply", attrs = c("function" = "lessThan"))
  c2 <- newXMLNode(name = "FieldRef", attrs = c("field" = "out"))
  c3 <- newXMLNode(name = "FieldRef", attrs = c("field" = "t"))
  addChildren(c1node, kids = c(c2, c3))
  c3 <- newXMLNode(name = "Constant", attrs = c("dataType" = "boolean"), text = TRUE)
  c4 <- newXMLNode(name = "Constant", attrs = c("dataType" = "boolean"), text = FALSE)
  expected <- addChildren(c0node, kids = c(c1node, c3, c4))

  current_split <- strsplit(saveXML(current), split = "")[[1]]
  expected_split <- strsplit(saveXML(expected), split = "")[[1]]
  expect_equal(current_split, expected_split)
})
