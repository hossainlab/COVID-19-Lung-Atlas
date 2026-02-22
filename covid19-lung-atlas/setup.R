# Install required packages for COVID-19 Lung Atlas Shiny App
required_packages <- c("shiny", "shinydashboard", "plotly", "DT",
                        "dplyr", "tidyr", "ggplot2", "viridis",
                        "data.table")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
  cat("Installed:", paste(new_packages, collapse = ", "), "\n")
} else {
  cat("All packages already installed.\n")
}

# Verify all packages load
invisible(lapply(required_packages, function(pkg) {
  library(pkg, character.only = TRUE)
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}))

cat("\nTo run the app:\n")
cat("  shiny::runApp('shiny_app.R')\n")
