require(pls)
require(dplyr)

dfmsc <- function(x) {
  as.data.frame(msc(as.matrix(x[,sapply(x, is.numeric)])))
}


setGeneric("group_msc", function(obj, fun = mean, ...)
    standardGeneric("group_msc"))

setMethod("group_msc", "Spectra",
          function(obj, fun = mean, ...){
            
            # making up an id name from the aggregation function
            id_fun <- as.character(substitute(fun, env = parent.frame()))[1]
            id_obj <- as.character(substitute(obj, env = parent.frame()))
            id <- paste(id_fun, id_obj, sep = '.')
            
            # applying the function to the spectra
            nir <- aaply(.data = spectra(obj), .margins = 2, .fun = fun, ...)
            
            # Create and return Spectra object
            Spectra(wl = wl(obj), nir = nir, id = id, units = wl_units(obj))
          }
)

# In the case of a SDF, an id can be given to split the SDF and apply fun
#' @rdname group_msc
setMethod("group_msc", "SpectraDataFrame",
          function(obj, fun = mean, id = NULL, ...){
            
            # No split --> the whole data is aggregated together
            if (is.null(id)) {
              # making up an id name from the aggregation function
              id_fun <- as.character(substitute(fun, env = parent.frame()))[1]
              id_obj <- as.character(substitute(obj, env = parent.frame()))
              
              # Select and paste only alphanumeric chars
              id_obj <- paste(id_obj[grep(x = id_obj, pattern = '[[:alnum:]]')], collapse = '.')
              # Combine object name and function name into an id
              id <- paste(id_fun, id_obj, sep = '.')
              
              # applying the function to the spectra
              nir <- as.data.frame(msc(obj@nir))
              
              res <- Spectra(wl = wl(obj), nir = nir, id = id, units = wl_units(obj))
              
              data <- features(obj)
              
              res <- SpectraDataFrame(res, data = data.frame(matrix(data, nrow = 1, dimnames = list(id, names(obj)))))
            }
            
            # There is a variable against which the data will be aggregated
            else {
              if (id %in% names(features(obj))) {
                
                # Col index of the splitting variable
                idx <- paste("(names(features(obj)) == '", id, "')", sep = "")
                idx <- paste(idx, collapse="|")
                idx <- which(eval(parse(text=idx)))
                
                # Creating spectra splits
                s <- data.frame(id = features(obj)[, idx, drop = FALSE], spectra(obj))
                colnames(s) <- gsub("id.", "", colnames(s))
                s[] <- lapply(s, function(x) if(is.factor(x)) factor(x) else x)
                idx2 <- paste("((names(s)) == '", id, "')", sep = "")
                idx2 <- paste(idx2, collapse="|")
                p <- split( s[-which(eval(parse(text=idx2)))], s[,id] )
                # Apply msc
                N <- length(p)
                l <- list()
                for (i in 1:N) {
                  l[[i]] <- (dfmsc(p[[i]])) %>% 
                    tibble::rownames_to_column('ref')
                }
                
                z <- do.call(rbind, l)
                z$ref <- as.numeric(z$ref)
                z <- z[order(z$ref),]
                
                # new data slot
                d <- subset(features(obj), rownames(features(obj)) %in% rownames(z))
                d$ref <- as.numeric(rownames(d))
                d <- d[order(d$ref),]
                
                #Reassemble spectra
                res <- z
                spectra(res) <- ref ~ ...
                features(res, safe=TRUE, key="ref", exclude_id=TRUE) <- d
              }
              else
                stop('Bad aggregation identifier.')
            }
            
            res
          }
)
                