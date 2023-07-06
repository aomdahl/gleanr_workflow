
standardizeColNames <- function(cin)
{
    #tidyr doesn't like some names, so just trying to make sure things are consistent....
    gsub(x = cin, pattern = "\\.", replacement = "-") %>% 
          gsub(x = ., pattern = "_-", replacement = "_") %>%
          gsub(x = ., pattern = "-$", replacement = "") %>% gsub( x= ., pattern = "^X(\\w)", replacement = "\\1")
}