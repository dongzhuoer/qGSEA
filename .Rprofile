if (tolower(Sys.getenv('CI')) != 'true') {
    if (file.exists('~/.Rprofile')) source('~/.Rprofile')

    # for R-raw/*.Rmd
    options(knitr.package.root.dir = normalizePath('./'))

    # disable CRAN checks in pkgdown (unless the package is released to CRAN)
    options(pkgdown.internet = FALSE)  

    # portable, `%>%` used in R-raw/*.Rmd. If user want to develop the package, we assumed he has installed all dependencies
    options('defaultPackages' = unique(c(getOption('defaultPackages'), 'magrittr')));  
}
