
# Get taxon ID function

# custom gbifid function to automatically choose the best order

get_gbifid2 <- function (sci, ask = TRUE, messages = TRUE, rows = NA, phylum = NULL, 
                         class = NULL, order = NULL, family = NULL, rank = NULL, method = "backbone", 
                         sciname = NULL, ...) 
{
  assert(sci, c("character", "taxon_state"))
  assert(ask, "logical")
  assert(messages, "logical")
  assert(phylum, "character")
  assert(class, "character")
  assert(order, "character")
  assert(family, "character")
  assert(rank, "character")
  assert(method, "character")
  assert_rows(rows)
  pchk(sciname, "sci")
  if (inherits(sci, "character")) {
    tstate <- taxon_state$new(class = "gbifid", names = sci)
    items <- sci
  }
  else {
    assert_state(sci, "gbifid")
    tstate <- sci
    sci <- tstate$taxa_remaining()
    items <- c(sci, tstate$taxa_completed())
  }
  prog <- progressor$new(items = items, suppress = !messages)
  done <- tstate$get()
  for (i in seq_along(done)) prog$completed(names(done)[i], 
                                            done[[i]]$att)
  prog$prog_start()
  for (i in seq_along(sci)) {
    direct <- FALSE
    mssg(messages, "\nRetrieving data for taxon '", 
         sci[i], "'\n")
    df <- switch(method, backbone = gbif_name_backbone(sci[i], 
                                                       ...), lookup = gbif_name_lookup(sci[i], ...))
    mm <- NROW(df) > 1
    if (is.null(df)) 
      df <- data.frame(NULL)
    if (nrow(df) == 0) {
      mssg(messages, m_not_found_sp_altclass)
      id <- NA_character_
      att <- "not found"
    }
    else {
      names(df)[1] <- "gbifid"
      id <- df$gbifid
      att <- "found"
    }
    if (length(id) == 0) {
      mssg(messages, m_not_found_sp_altclass)
      id <- NA_character_
      att <- "not found"
    }
    if (length(id) > 1) {
      matchtmp <- df[as.character(df$canonicalname) %in% 
                       sci[i], "gbifid"]
      if (length(matchtmp) == 1) {
        id <- as.character(matchtmp)
        direct <- TRUE
      }
      else {
        if (!is.null(phylum) || !is.null(class) || !is.null(order) || 
            !is.null(family) || !is.null(rank)) {
          df <- filt(df, "phylum", phylum)
          df <- filt(df, "class", class)
          df <- filt(df, "order", order)
          df <- filt(df, "family", family)
          df <- filt(df, "rank", rank)
        }
        df <- sub_rows(df, rows)
        if (NROW(df) == 0) {
          id <- NA_character_
          att <- "not found"
        }
        else {
          id <- df$gbifid
          if (length(id) == 1) {
            rank_taken <- as.character(df$rank)
            att <- "found"
          }
        }
        if (length(id) > 1) {
          if (ask) {
            df <- df[, switch(method, backbone = gbif_cols_show_backbone, 
                              lookup = gbif_cols_show_lookup)]
            message("\n\n")
            message("\nMore than one GBIF ID found for taxon '", 
                    sci[i], "'!\n\n            Enter rownumber of taxon (other inputs will return 'NA'):\n")
            rownames(df) <- 1:nrow(df)
            # print(df)
            # take <- scan(n = 1, quiet = TRUE, what = "raw")
            take <- (df$rank == "order" & df$status == "ACCEPTED" & df$matchtype == "EXACT")
            
            if (length(take) == 0) {
              take <- "notake"
              att <- "nothing chosen"
            }
            if (take %in% seq_len(nrow(df))) {
              take <- as.numeric(take)
              message("Input accepted, took gbifid '", 
                      as.character(df$gbifid[take]), "'.\n")
              id <- as.character(df$gbifid[take])
              att <- "found"
            }
            else {
              id <- NA_character_
              att <- "not found"
              mssg(messages, "\nReturned 'NA'!\n\n")
            }
          }
          else {
            if (length(id) != 1) {
              warning(sprintf(m_more_than_one_found, 
                              "gbifid", sci[i]), call. = FALSE)
              id <- NA_character_
              att <- m_na_ask_false
            }
          }
        }
      }
    }
    res <- list(id = id, att = att, multiple = mm, direct = direct)
    prog$completed(sci[i], att)
    prog$prog(att)
    tstate$add(sci[i], res)
  }
  out <- tstate$get()
  ids <- structure(as.character(unlist(pluck(out, "id"))), 
                   class = "gbifid", match = pluck_un(out, "att", 
                                                      ""), multiple_matches = pluck_un(out, "multiple", 
                                                                                       logical(1)), pattern_match = pluck_un(out, "direct", 
                                                                                                                             logical(1)))
  on.exit(prog$prog_summary(), add = TRUE)
  on.exit(tstate$exit, add = TRUE)
  add_uri(ids, get_url_templates$gbif)
}

environment(get_gbifid2) <- asNamespace('taxize')


# args:

# database_function: name of the function - "itis" (get_tsn), "bold" - (get_boldid), "gbif" - (get_gbifid)
# taxon_name: name that you want to query
# ask_or_not: accept imperfect matches (default is FALSE)
# tries: how many times to try the function if it keeps getting an error (default = 5)

get_taxon_id <- function(database_function = "itis", taxon_name, ask_or_not = FALSE, tries = 5) {
  
  # choose the correct function based on the database
  if(database_function == "itis") {
    
    func_string <- "get_tsn"
    
  } else if (database_function == "bold") {
    
    func_string <- "get_boldid"
    
  } else if (database_function == "gbif") {
    
    func_string <- "get_gbifid2"
    
  } else {
    
    stop("error, choose a supported database: itis, bold, gbif")
    
  }
  
  library(taxize)
  x <- try(stop("!"), silent = TRUE)
  i <- 1
  while( class(x) == "try-error" ) {
    
    if (i > tries){
      break
    }
    
    i <- i + 1
    x <- try( do.call(func_string, list(sci = taxon_name, ask = ask_or_not)) )
    
  }
  
  if ( class(x) == "try-error" ) {
    
    attr(x, "match") <- "could not access database"
    warning("could not access database, use more tries or run function at a later stage")
    
  }
  
  return(x)
  
}


# define a function to extract the genus from a binomial

# args
# binomial - binomial character string separated by a space (e.g. "Loxodonta loxodonta")

extract_genus <- function(binomial) {
  z <- unlist( strsplit(x = binomial, split = " ", fixed = TRUE) )
  if (length(z) > 1) {
    search.name <- z[1]
  } else {search.name <- binomial}
  return(search.name)
}

### END
