context("RNA-seq analysis basic workfolw, including DE analysis")
test_that("RNA-seq analysis basic workfolw", {
  input = create_count_phe()
  count_data = input[[1]]
  col_data = input[[2]]
  expect_true(any(class(count_data)%in%c("matrix","array")))
  expect_equal(class(col_data),"data.frame")

  dds = create_DEseq(input[[1]],input[[2]],ref_level="untreated")
  expect_is(dds,"DESeqDataSet")

})
