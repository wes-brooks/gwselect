epanechnikov = function(R, bw) {
  ifelse( R < bw, 1-(R/bw)**2, 0)
}

