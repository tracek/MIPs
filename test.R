vect_lag <- function(v, n=1, forward=FALSE) {
    if (forward)
        c(v[(n+1):length(v)], rep(NA, n))
    else
        c(rep(NA, n), v[1:(length(v) - n)])
}

mydf <- data.frame(x=1:10, y=vect_lag(1:10, 1, forward = TRUE))

for (i in 1:(nrow(mydf)-2))
{
    print(mydf[i,]$x)
}