folder <- tempfile("SDMtune")
settings <- list(update = FALSE)
data <- list()
.create_chart(folder = folder, script = "gridSearch.js", settings = settings,
              data = data)
.update_data(folder, data = settings)

test_that(".create_chart function works correctly", {
  # The folder lib is created and contains the correct files
  expect_true(file.exists(file.path(folder, "lib", "chart_script.js")))
  expect_true(file.exists(file.path(folder, "lib", "Chart.min.js")))
  expect_true(file.exists(file.path(folder, "lib", "jquery.min.js")))
  expect_true(file.exists(file.path(folder, "lib", "style.css")))
  # The template file is created
  expect_true(file.exists(file.path(folder, "chart_template.html")))
})

test_that(".render_script function render settings and data", {
  # Settings is rendered
  expect_equal(readLines(file.path(folder, "lib", "chart_script.js"),
                         encoding = "UTF-8")[2],
               "var settings = {\"update\":[false]};")
  # Data is rendered
  expect_equal(readLines(file.path(folder, "lib", "chart_script.js"),
                         encoding = "UTF-8")[3], "var data = [];")
})

test_that(".update_data function works corretly", {
  # The data.json file is created
  expect_true(file.exists(file.path(folder, "data.json")))
  # The data.json file contains the data
  expect_equal(readLines(file.path(folder, "data.json"), encoding = "UTF-8"),
               "{\"update\":[false]}")
})

teardown(unlink(folder, recursive = TRUE))
