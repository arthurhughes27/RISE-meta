# File to detect packages required for the project and installs anything missing on the user's machine
deps <- renv::dependencies()
pkgs <- unique(na.omit(deps$Package))

for(pkg in pkgs){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

rm(list = ls())