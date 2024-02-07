#' Filter a dataframe to keep just half the most variable genes
#'@description
#' `Get_mostvar` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

Get_mostvar <- function(df, n){
  
  abdv <- apply(df, 1, function(x) {
    sum(
      abs(
        x - mean(x)
      )
    ) / length(x)
  })
  
  # to input either % or absolute n
  if ( n>0 && n<1){
    th_abdv.i <- nrow(df)*n
  } else{
    th_abdv.i <- n
  }
  
  abdv %>% sort(decreasing = T) %>% .[th_abdv.i] -> th_abdv
  df.f <- df[abdv >= th_abdv, , drop=F]
  return(df.f)
}