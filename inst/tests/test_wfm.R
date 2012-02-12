context('Check word frequency matrix functions')

load('mats.RData') # rmat (matrix) and mmat (Matrix)

## string split on whitespace with sane defaults
spl <- function(x){ s = strsplit(x, '[ ]+'); if (length(s)==1){ s[[1]] } else { s } }

# function(wfm, min.count=5, min.doc=5, sample.docs=NULL, sample.words=NULL, verbose=TRUE){

test_that("trim: min.count", {
  rmat.out <- trim(rmat, min.count=50, min.doc=1, verbose=FALSE) # W2 occurs just 49 times 	
  mmat.out <- trim(mmat, min.count=50, min.doc=1, verbose=FALSE)
  expect_equal(colnames(rmat.out), spl("W1 W3 W4"))
  expect_equal(colnames(mmat.out), spl("W1 W3 W4"))  
})

test_that("trim: min.doc", {
  rmat.out <- trim(rmat, min.count=1, min.doc=3, verbose=FALSE) # W4 does not appear in D1 	
  mmat.out <- trim(mmat, min.count=1, min.doc=3, verbose=FALSE)
  expect_equal(colnames(rmat.out), spl("W1 W2 W3"))
  expect_equal(colnames(mmat.out), spl("W1 W2 W3"))  
  rmat.out <- trim(rmat, min.count=1, min.doc=2, verbose=FALSE) # should keep them all
  mmat.out <- trim(mmat, min.count=1, min.doc=2, verbose=FALSE)
  expect_equal(rmat.out, rmat)
  expect_equal(mmat.out, mmat)  
})

test_that("trim: min.count and min.doc", {
  set.seed(12345)
  rmat.out <- trim(rmat, min.count=50, min.doc=3, verbose=FALSE) # toss W2 and W4
  mmat.out <- trim(mmat, min.count=50, min.doc=3, verbose=FALSE) 
  expect_equal(colnames(rmat.out), spl("W1 W3"))
  expect_equal(colnames(mmat.out), spl("W1 W3"))  
})

test_that("trim: sample.docs", {
  set.seed(12345)
  rmat.out <- trim(rmat, min.count=1, min.doc=1, sample.docs=2, verbose=FALSE) # keep everything
  mmat.out <- trim(mmat, min.count=1, min.doc=1, sample.docs=2, verbose=FALSE) # keep everything
  expect_equal(rownames(rmat.out), spl("D3 D2"))
  expect_equal(rownames(mmat.out), spl("D3 D2"))  
})

test_that("ldac: round tripping", {
  unlink('ldac-out', TRUE) ## prepare
  write_ldac(mmat, 'ldac-out')
  mmat.recon <- read_ldac('ldac-out')  
  expect_equal(mmat.recon, mmat)
  unlink('ldac-out', TRUE)
})

test_that("mtx: round tripping", {
  require(Matrix)
  unlink('mtx-out') ## prepare
  writeMM(mmat, 'mtx-out')
  mmat.recon <- read_mtx_format('mtx-out')
  expect_equal(mmat.recon, mmat)
  mmat.recon.2 <- read_mtx('mtx-out', jfreq=FALSE)  
  expect_equal(mmat.recon.2, mmat)
  unlink('mtx-out')
})

