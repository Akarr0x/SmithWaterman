test_that("error DNAString seq1", {
  seq1 <- 'GAATC'
  seq2 <- DNAString("CATACG")
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("error DNAString seq2", {
  seq1 <- DNAString('GAATC')
  seq2 <- 'CATACG'
  expect_error(globalSequenceAlignment(seq1, seq2))
})


test_that("error DNAString length 1", {
  seq1 <- DNAString('')
  seq2 <- DNAString('AACTTG')
  expect_error(globalSequenceAlignment(seq1, seq2))
})


test_that("error DNAString length 1", {
  seq1 <- DNAString('ATCCG')
  seq2 <- DNAString('')
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("exact result", {
    y <- c("*", " ", "*", " ", "*", "*")
    expect_equal(smithwaterman(DNAString("ACCTG"),
                 DNAString("ATCTG"))$result, y)
})

test_that("exact result", {
    y <- c("*", " ", "*", " ", "*", "*")
    expect_equal(smithwaterman(DNAString("ACCTG"),
                 DNAString("ATCTG"))$score, 10)
})

test_that("match score higher than zero", {
    expect_error(smithwaterman(DNAString("ACCTG"),
                 DNAString("ATCTG"), match = 0))
})


test_that("missmatch score lower than zero", {
    expect_error(smithwaterman(DNAString("ACCTG"),
                 DNAString("ATCTG"), missmatch = 2))
})


test_that("gap score lower than zero", {
    expect_error(smithwaterman(DNAString("ACCTG"),
                 DNAString("ATCTG"), gap = 1))
})