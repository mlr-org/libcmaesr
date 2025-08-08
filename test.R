library(devtools)
library(roxygen2)
document()
roxygenize()
load_all()

f = function(x) {
  print(x)
  apply(x, 1, function(x) sum(x^2))
}
ctrl = cmaes_control(max_fevals = 10)
print(ctrl)
zz = cmaes(objective = f, x0 = 1, lower = 0, upper = 10, control = ctrl)
print(zz)








