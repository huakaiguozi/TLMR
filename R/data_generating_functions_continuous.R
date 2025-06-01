logistic_function <- function(x) {
  return(1 / (1 + exp(-x)))
}


dimention_C <- function(p){
  list_1 <- rnorm(mean = p, sd = p/2 ,n = 2)
  list_2 <- round(list_1/sum(list_1) * p*(2/3))
  list_3 <- c(list_2,(p-sum(list_2)))
  return(list_3)
}
