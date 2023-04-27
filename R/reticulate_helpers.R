# global reference to scipy (will be initialized in .onLoad)
scipy <- NULL
numpy <- NULL
hdbscan <- NULL

#' onLoad
#'
#' @param libname
#' @param pkgname
#'
#' @return
#' @export
#'
#' @examples
.onLoad <- function(libname, pkgname) {

  reticulate::configure_environment(pkgname, force = TRUE)

  hdbscan <<- reticulate::import('hdbscan', delay_load = TRUE)

  }

#'  Install python packages to use `HDBSCAN`
#'
#' @param method
#' @param conda
#'
#' @return
#' @export
#'
#' @examples
install_python_packages <- function() {

  if(!reticulate::virtualenv_exists('hdbscanR')){
    if(!dir.exists(reticulate::miniconda_path())){
      reticulate::install_miniconda()
    }
    reticulate::install_python('3.9:latest')
    reticulate::virtualenv_create(envname = 'hdbscanR',
                                  version = '3.9:latest',
                                  packages = c('numpy','hdbscan'))

    reticulate::use_virtualenv('hdbscanR')
  } else {
    reticulate::install_python('3.9:latest')
    reticulate::use_virtualenv('hdbscanR')
  }
  reticulate::configure_environment('hdbscanR', force = TRUE)
}

