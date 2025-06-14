# NB: RandomTreeScore uses C's RNG, so no point in setting seed.
MorphyAction <- function (Action) expect_equal("ERR_NO_ERROR", mpl_translate_error(Action))
MorphyWith <- function (char) {
  nTip <- nchar(char) - 1L
  morphyObj <- mpl_new_Morphy()
  MorphyAction(mpl_init_Morphy(nTip, 1, morphyObj)) 
  MorphyAction(mpl_attach_rawdata(char, morphyObj)) 
  MorphyAction(mpl_set_num_internal_nodes(nTip - 1L, morphyObj)) 
  MorphyAction(mpl_set_parsim_t(1, "FITCH", morphyObj))
  MorphyAction(mpl_set_charac_weight(1, 1, morphyObj)) 
  MorphyAction(mpl_apply_tipdata(morphyObj))
  class(morphyObj) <- "morphyPtr"
  morphyObj
}


context("pp: Tree randomness")
test_that("four-tip trees are randomly distributed", {
  nTrees <- 36000
  stringency <- 0.005 # low numbers mean you'll rarely fail by chance
  nTip <- 4
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees, 1/(nTip - 1))
  rTrees <- vapply(logical(nTrees), function (XX) 
    unlist(RandomMorphyTree(nTip)), integer((nTip * 4) - 3))
  expect_true(all(rTrees[1 + (seq_len(nTip - 1)), ] %in% nTip + seq_len(nTip - 2)))
  expect_lt(expectedBounds[1], sum(rTrees[2, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[2, ] == 5))
  expect_lt(expectedBounds[1], sum(rTrees[3, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[3, ] == 5))
  expect_lt(expectedBounds[1], sum(rTrees[4, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[4, ] == 5))

  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] > expectedBounds[1]))
  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] < expectedBounds[2]))

  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] < nTrees - expectedBounds[1]))
  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] > nTrees - expectedBounds[2]))
})

test_that("four-tip trees are randomly scored", {
  set.seed(0)
  
  nTrees <- 6000
  stringency <- 0.005
  nTip <- 4
  
  morphyObj <- MorphyWith("0011;")
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  
  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees,
                           NUnrooted(nTip - 1L) / NUnrooted(nTip))
  scores <- vapply(logical(nTrees), 
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  expect_lt(expectedBounds[1], sum(scores==1))
  expect_gt(expectedBounds[2], sum(scores==1))
})

test_that("five-tip trees are randomly scored", {
  set.seed(0)
  nTrees <- 6000
  stringency <- 0.005
  nTip <- 5
  morphyObj <- MorphyWith("00011;")
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           NUnrooted(nTip - 1) / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  expect_equal(2L, max(scores))
  expect_lt(expectedBounds[1], sum(scores == 1))
  expect_gt(expectedBounds[2], sum(scores == 1))
})


test_that("six-tip trees are randomly scored", {
  set.seed(0)
  
  nTrees <- 6000
  stringency <- 0.005
  nTip <- 6
  
  morphyObj <- MorphyWith("000011;")
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           NUnrooted(5) / NUnrooted(6))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  morphyObj <- UnloadMorphy(morphyObj)
  
  expect_true(max(scores) == 2)
  expect_lt(expectedBounds[1], sum(scores==1))
  expect_gt(expectedBounds[2], sum(scores==1))
  
  morphyObj <- MorphyWith("001122;")
  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees,
                           7 / NUnrooted(nTip))
  scores <- vapply(logical(nTrees), 
                   function (XX) RandomTreeScore(morphyObj),
                   integer(1))
  morphyObj <- UnloadMorphy(morphyObj)
  
  expect_true(all(scores %in% 2:4))
  expect_lt(expectedBounds[1], sum(scores == 2))
  expect_gt(expectedBounds[2], sum(scores == 2))
  
  morphyObj <- MorphyWith("000111;")
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           3 * 3 / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  # unloaded on exit; don't unload twice || morphyObj <- UnloadMorphy(morphyObj)
  
  expect_true(max(scores) == 3)
  expect_lt(expectedBounds[1], sum(scores == 1))
  expect_gt(expectedBounds[2], sum(scores == 1))
  
})

test_that("twelve-tip trees are randomly scored", {
  nTrees <- 12000 # 12000 seems to throw false +ve too often?
  stringency <- 0.01 #  increased from 0.005 to avoid false +ves
  nTip <- 12
  morphyObj <- MorphyWith("000000011111;")
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees, 
                           NUnrooted(7) * (2 * 7 - 3) *
                           NUnrooted(5) * (2 * 5 - 3) / NUnrooted(nTip))
  
  scores <- vapply(logical(nTrees), 
                   function (XX) RandomTreeScore(morphyObj),
                   integer(1L))
  # table(scores)
  
  expect_equal(5L, max(scores))
  nScoring1 <- sum(scores == 1)
  expect_lt(expectedBounds[1], nScoring1)
  expect_gt(expectedBounds[2], nScoring1)  
})
