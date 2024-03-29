test_that("LengthAdded() errors", {
  trees <- inapplicable.trees[["Vinther2008"]]
  dataset <- inapplicable.phyData[["Vinther2008"]]
  
  expect_error(
    LengthAdded(trees, "dataset"),
    "`char` must be a character of class `phyDat`"
  )
  
  expect_error(
    LengthAdded(trees, dataset),
    "`char` must comprise a single character"
  )
  
  attr(dataset, "contrast")[6, ] <- 0
  expect_error(
    LengthAdded(trees, dataset[, 51]),
    "`char` contract matrix lacks levels for 6"
  )
})

test_that("LengthAdded()", {
  trees <- inapplicable.trees[["Vinther2008"]]
  dataset <- inapplicable.phyData[["Vinther2008"]]
  
  pe10 <- LengthAdded(trees, dataset[, 10])
  expect_equal(pe10["Neopilina"], c(Neopilina = 1))
  expect_equal(sum(pe10), 1)
  
  # Single tree
  expect_equal(LengthAdded(trees[[1]], dataset[, 10]), pe10)
  
  # No inapplicables
  appData <- dataset
  colnames(attr(appData, "contrast"))[1] <- "x"
  attr(appData, "levels")[1] <- "x"
  attr(appData, "allLevels")[4] <- "x"
  pe10 <- LengthAdded(trees, appData[, 10])
  expect_equal(pe10["Neopilina"], c(Neopilina = 1))
  expect_equal(sum(pe10), 1)
  
  # Implied weighting
  expect_equal(
    unname(PolEscapa(trees, dataset[, 11], concavity = 5)["Neopilina"]),
    as.numeric(TreeLength(trees[[1]], dataset[, 11], concavity = 5))
  )
  
  # minLength changes when only occurrence of 1 -> ?
  wiwaxia <- LengthAdded(trees, dataset[, 39], concavity = 10)
  expect_true(all(wiwaxia >= 0))
  
})
