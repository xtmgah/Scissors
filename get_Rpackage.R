library("devtools"); library("roxygen2");
package_path = "~/Dropbox/Research/hyochoi_GitHub/"

setwd(package_path)
# create("Scissors")

# Save functions
# Add documentation

setwd("./Scissors")
document()

setwd("..")
install("Scissors")

library(Scissors)
?annotate_pileup
