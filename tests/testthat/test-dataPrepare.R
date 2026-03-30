

test_that("BSRDataModel", {
  data(sdc, package = "BulkSignalR")
  bsrdm <- BSRDataModel(sdc)
  expect_s4_class(bsrdm,"BSRDataModel")
  
})
