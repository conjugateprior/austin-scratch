### --- Test setup ---

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library("RUnit")
    library("austin")
}

a <- 1
b <- 2

### --- Test functions ---

test.simple <- function(){
    checkTrue(a < b)
}

