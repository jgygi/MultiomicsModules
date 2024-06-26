# Define SPEAR class:
MOModules <- R6::R6Class("MOModules",
                           public = list(
                             modules = NULL,
                             ID_tables = NULL,
                             scores = list(),
                             params = list(),
                             reports = list(),
                             
                             # Functions:
                             initialize = function(
                               path.to.modules = NULL,
                               quiet = FALSE
                             ) {
                               # called by MOModules$new(...)
                               
                               # Params
                               self$params$quiet = quiet
                               self$params$color_scheme = c("#5558F7", "#00BB13", "#FF8000")
                               names(self$params$color_scheme) = c("Transcriptomics", "Proteomics", "Metabolomics")
                               
                               # Reports:
                               self$reports$map = list()
                               self$reports$percentage = list()
                               self$reports$fischer = list()
                               
                               # Initialize Object
                               if(!self$params$quiet){
                                 cat("----------------------------------------------------------------\n")
                                 cat("Multi-omics Modules R Tool version 1.0.0\nPlease direct all questions to Jeremy Gygi\n(jeremy.gygi@yale.edu) or Pramod Shinde (pshinde@lji.org)\n", sep = "")
                                 cat("----------------------------------------------------------------\n")
                                 cat("Generating Multi-omics Modules object...\n")
                               }
                               
                               # Data:
                               # Load RDS object and populate
                               # NOTE: Replace with 'data(...)' when making package
                               if(!is.null(path.to.modules)){load(path.to.modules)}
                               else{
                                 data("mo_modules_v1")
                               }
                               tmp = mo_modules_v1
                               self$modules = tmp$Modules
                               self$ID_tables = tmp$ID_Tables
                             },
    
    # print function:
    print = function(...){
      cat("Multi-omics Modules Object v 1.0\nTODO: Update print statement")
      return(invisible(self)) 
    },
    
    # Method functions: (each found in 1_helpers.R)
    add_analyte_scores = add_analyte_scores,  
    convert_id = convert_id,
    get_module_scores = get_module_scores,
    plot_fingerprint = plot_fingerprint,
    plot_mapping = plot_mapping,
    plot_report = plot_report
    
                           ), # end public
    private = list(
    )
)

#' Make a Multi-omics Modules object. Will return an R6 class MOModules used for the "MOModules" package.
#'@param quiet Mute print statements? Defaults to `FALSE`.
#'@export
make.multiomics.modules <- function(
    path.to.modules = NULL,
    quiet = FALSE
){
  return(MOModules$new(
    path.to.modules = path.to.modules,
    quiet = quiet
  ))
}

