# Basic tests for TiNDA package

test_that("TiNDA accepts valid input", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 100)

  expect_no_error({
    result <- TiNDA(test_df, sample_name = "test_1")
  })

  expect_s3_class(result, "TiNDA")
})

test_that("TiNDA validates input columns", {
  bad_df <- data.frame(
    CHR = c(1, 2),
    POS = c(100, 200)
  )

  expect_error(
    TiNDA(bad_df),
    "Missing required columns"
  )
})

test_that("TiNDA validates data_source parameter", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 100)

  expect_error(
    TiNDA(test_df, data_source = "invalid"),
    "data_source must be either"
  )
})

test_that("TiNDA validates input is data frame", {
  expect_error(
    TiNDA(list(a = 1)),
    "tbl must be a data frame"
  )
})

test_that("TiNDA returns expected structure", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 100)
  result <- TiNDA(test_df, sample_name = "test_1")

  expect_true("data" %in% names(result))
  expect_true("sample_name" %in% names(result))
  expect_true("classification_summary" %in% names(result))
  expect_true("parameters" %in% names(result))
})

test_that("TiNDA classification includes expected classes", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 200)
  result <- TiNDA(test_df, sample_name = "test_1")

  classes <- unique(result$data$TiN_Class)
  expect_true("Germline" %in% classes || "Somatic_Rescue" %in% classes)
})

test_that("print.TiNDA works without error", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 100)
  result <- TiNDA(test_df, sample_name = "test_1")

  expect_no_error(print(result))
})

test_that("summary.TiNDA returns expected fields", {
  data(hg19_length)
  test_df <- generate_test_data(hg19_length, num_variants = 100)
  result <- TiNDA(test_df, sample_name = "test_1")
  summ <- summary(result)

  expect_true("total_variants" %in% names(summ))
  expect_true("classification" %in% names(summ))
})

test_that("get_tinda_params returns correct defaults", {
  params <- get_tinda_params("WGS")

  expect_equal(params$max_control_af, 0.25)
  expect_equal(params$min_tumor_af, 0.01)
  expect_equal(params$min_clst_members, 0.85)
})

test_that("get_tinda_params validates data_source", {
  expect_error(
    get_tinda_params("invalid"),
    "data_source must be either"
  )
})
