context("pp_info_extra_step.R")

test_that("Information content of steps calculated correctly", {
  expect_equal(c(1, 2),
               as.double(ICS(2, 2, 10000, warn = FALSE) * NUnrooted(4)))
  expect_equal(c(3, 12), tolerance = 1e-5,
               as.double(ICS(2, 3, 10000, warn = FALSE) * NUnrooted(5)))
  
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(1)
  expect_equal(cumsum(as.double(ICS(3,3, 6000, warn = FALSE) * NUnrooted(6))),
               c(NUnrootedMult(c(3,3)),
                 NUnrootedMult(c(3,3)) + WithOneExtraStep(c(3,3)),
                 NUnrooted(6)),
               tolerance = 1)
    
  expect_equal(cumsum(as.double(ICS(3, 12, 60000, warn = FALSE))),
               c(NUnrootedMult(c(3,12)) / NUnrooted(3 + 12),
                 (NUnrootedMult(c(3,12)) + WithOneExtraStep(c(3,12)))
                 / NUnrooted(3 + 12),
                 1),
               tolerance = 5e-03)
})