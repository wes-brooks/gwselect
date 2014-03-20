epanechnikov = function(R, bw) {
  ifelse( R < bw, 0.75*(1 - (R/bw)**2), 0)
}

