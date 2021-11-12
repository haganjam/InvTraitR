
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# this function takes two vectors and determines if they are nearly equal etc.

# https://stackoverflow.com/questions/25154930/less-or-equal-for-floats-in-r
near_equal <- function( x , y , tol = 1.5e-8 , mode = "ae" ){
  ae <- mapply( function(x,y) isTRUE( all.equal( x , y , tolerance = tol ) ) , x , y )    
  gt <- x > y
  lt <- x < y
  if( mode == "ae" )
    return( ae )
  if( mode == "gt" )
    return( gt )
  if( mode == "lt" )
    return( lt )
  if( mode == "ne.gt" )
    return( ae | gt )
  if( mode == "ne.lt" )
    return( ae | lt )
}

### END
