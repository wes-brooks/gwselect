cubic = function(d, bw) {
    ifelse(d<bw, 1 - (7*(d/bw)**2 - 8.75*(d/bw)**3 + 3.5*(d/bw)**5 - 0.75*(d/bw)**7), 0)
}