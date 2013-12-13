spherical = function(d, bw) {
    ifelse(d<bw, 1 - 1.5*(d/bw) + 0.5*(d/bw)**3, 0)
}