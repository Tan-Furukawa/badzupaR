library(magrittr)

mapList <- function(list, mapFn = function(elem, index=NULL, array=NULL) {
  print(index)
  print(elem)
}) {
  seq_along(list) %>% lapply(function(i) {
    mapFn(list[[i]], i, list)
  })
}

group <- function(df, keyFn = function(df) df$key) {
  keys <- unique(keyFn(df))
  keys %>%
    lapply(function(key) {
      return (df[keyFn(df) == key,])
    })
}

set.seed(1234); testDat1 <- c(
  rnorm(100, mean = 100, sd = 10),
  rnorm(50, mean = 300, sd = 100)
)


