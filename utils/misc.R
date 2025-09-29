check_and_install_package <- function(package_name) {
  # Check if the package is already installed.
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    # If not installed, try to install it.
    message(paste("Package '", package_name, "' not found. Attempting to install...", sep = ""))
    
    # Use tryCatch to handle potential installation errors.
    tryCatch({
      if (package_name!='TwoSampleMR'){
        install.packages(package_name, dependencies = TRUE)
      }else if(package_name=='TwoSampleMR'){
        library(remotes)
        remotes::install_github("MRCIEU/")
      }else if(package_name=='coloc'){
        library(remotes)
        remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
      }
      
      # Try to load the package after installation.
      if (require(package_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Package '", package_name, "' installed and loaded successfully.", sep = ""))
      } else {
        # This case is rare, but good to have a fallback.
        stop("Package installed, but could not be loaded. Please check your R environment.")
      }
    }, error = function(e) {
      message(paste("Error installing package '", package_name, "': ", e$message, sep = ""))
    })
  } else {
    # If the package is already installed, just load it.
    message(paste("Package '", package_name, "' is already installed and loaded.", sep = ""))
  }
}


set_env <- function(){
  my_token = readline('Copy and paste your OpenGWAS JWT token here: ')
  file.create('./.Renviron')
  cat(paste0("OPENGWAS_JWT=",my_token), file = ".Renviron", append = F)
  
  # R should be able to pick up the .Renviron file in the home directory. If not, load the environment manually:
  readRenviron("./.Renviron")
  # Check if recognised. A long list of token should be returned
  ieugwasr::get_opengwas_jwt()
  # Check if working. Should return no error
  ieugwasr::user()
}