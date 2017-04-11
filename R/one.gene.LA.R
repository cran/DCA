one.gene.LA <-
function(array, x) # x is a row of the array, or a factor
{
    b<-x * t(array)
    d<-array %*% b
    d
}
