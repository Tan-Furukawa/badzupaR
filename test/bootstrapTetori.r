library(this.path)
cur_dir = dirname(this.path())
densFn <- function(x) IsoplotR::kde(x, plot=FALSE)

d <- read.csv(paste(cur_dir,'/data/tetori.csv', sep=''),header = T)
ageData <- d[,1]

boot <- Bootstrap$new(400, ageData, densFn, 0.1)
res <- boot$doAllProcess(show=FALSE)
