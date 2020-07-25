qgsea_addin <- function() {
    app <- source(system.file('shiny/app.R', package = 'qGSEA'), local = T)$value;
    shiny::runGadget(app, viewer = shiny::browserViewer())
}
