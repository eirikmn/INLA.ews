if(!rlang::is_installed("INLA")){
  warning(paste0("Package INLA not found. This can be installed by running: \n",
                 "install.packages(\"INLA\",repos=c(getOption(\"repos\"),INLA=\"https://inla.r-inla-download.org/R/stable\"), dep=TRUE)",
                 "\nfor the stable version or \n"),
          "install.packages(\"INLA\",repos=c(getOption(\"repos\"),INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)",
          "\nfor the testing version")
}
