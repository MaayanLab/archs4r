.onAttach <- function(libname, pkgname) {
    # Check and install BiocManager if not present
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        message("Installing BiocManager to manage Bioconductor dependencies...")
        install.packages("BiocManager")
    }
    
    # Check and install rhdf5 if not present
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        message("Installing required Bioconductor package: rhdf5")
        BiocManager::install("rhdf5", update = FALSE, ask = FALSE)
    }
    
    # Check and install preprocessCore if not present
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        message("Installing required Bioconductor package: preprocessCore")
        BiocManager::install("preprocessCore", update = FALSE, ask = FALSE)
    }
}