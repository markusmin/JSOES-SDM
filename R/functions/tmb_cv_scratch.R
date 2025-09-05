# cross-validation for bespoke TMB models. Adapted from the sdmTMB_cv() function from
# the sdmTMB package (Anderson et al. 2022)

# for testing
formula = n_per_km ~ 0 + as.factor(year)
data = csyif
mesh = inla_mesh_cutoff15_sdmTMB
time = "year"
spatial = "on"
spatiotemporal = "off"
family = tweedie(link = "log")
k_folds = 4

sdmTMB_cv <- function(
    formula, data, mesh_args, mesh = NULL, time = NULL,
    k_folds = 8, fold_ids = NULL,
    lfo = FALSE,
    lfo_forecast = 1,
    lfo_validations = 5,
    parallel = TRUE,
    use_initial_fit = FALSE,
    future_globals = NULL,
    ...) {
  if (k_folds < 1) cli_abort("`k_folds` must be >= 1.")
  
  spde <- mesh
  data[["_sdm_order_"]] <- seq_len(nrow(data))
  constant_mesh <- missing(mesh_args)
  if (missing(mesh_args)) mesh_args <- NULL
  if (missing(spde)) spde <- NULL
  if (lfo) fold_ids <- NULL
  # add column of fold_ids stratified across time steps
  if (is.null(time)) {
    time <- "_sdmTMB_time"
    data[[time]] <- 0L
  }
  if (is.null(fold_ids)) {
    if (lfo) {
      if (length(unique(data[[time]])) < (lfo_validations + lfo_forecast)) {
        cli_abort("Not enough time steps for the desired validation period. Either decrease `lfo_validations` or add more data.")
      }
      # Create lfo_validations + 1 folds, ordered sequentially
      data$cv_fold <- 1
      t_validate <- sort(unique(data[[time]]), decreasing = TRUE)
      for (t in seq(1, lfo_validations+lfo_forecast)) {
        # fold id increasing order + forecast
        data$cv_fold[data[[time]] == t_validate[t]] <- lfo_validations - t + 1 + lfo_forecast
      }
    } else {
      dd <- lapply(split(data, data[[time]]), function(x) {
        x$cv_fold <- sample(rep(seq(1L, k_folds), nrow(x)), size = nrow(x))
        x
      })
      data <- do.call(rbind, dd)
    }
    fold_ids <- "cv_fold"
  } else {
    # fold_ids passed in; can be numeric, or a named column in `data`
    data$cv_fold <- NA
    if (length(fold_ids) == nrow(data)) {
      data$cv_fold <- fold_ids
    }
    if (length(fold_ids) == 1L && is.character(fold_ids)) {
      if (!fold_ids %in% names(data)) {
        cli_abort("Name of fold identifier not found in data.")
      }
      data$cv_fold <- data[[fold_ids]]
    }
    if (length(fold_ids) > 1 && length(fold_ids) < nrow(data)) {
      cli_abort("Dimension of `fold_ids` doesn't match data and is not a named variable.")
    }
    if (length(which(is.na(data$cv_fold))) > 0) {
      cli_abort("NAs found in `fold_ids`; please check `fold_ids` are specified correctly.")
    }
    k_folds <- length(unique(data$cv_fold))
  }
  if (time == "_sdmTMB_time") { # undo changes above, make time NULL
    data[["_sdmTMB_time"]] <- NULL
    time <- NULL
  }
  
  dot_args <- as.list(substitute(list(...)))[-1L]
  if ("weights" %in% names(dot_args)) {
    cli_abort("`weights` cannot be specified within sdmTMB_cv().")
  }
  if ("offset" %in% names(dot_args)) {
    if (!is.character(dot_args$offset)) {
      cli_abort("Please use a character value for 'offset' (indicating the column name) for cross validation.")
    }
    .offset <- eval(dot_args$offset)
  } else {
    .offset <- NULL
  }
  
  if (k_folds > 1) {
    # data in kth fold get weight of 0:
    weights <- ifelse(data$cv_fold == 1L, 0, 1)
  } else {
    weights <- rep(1, nrow(data))
  }
  if (lfo) weights <- ifelse(data$cv_fold == 1L, 1, 0)
  
  if (use_initial_fit) {
    # run model on first fold to get starting values:
    
    if (!constant_mesh) {
      if (lfo) {
        dat_fit <- data[data$cv_fold == 1L, , drop = FALSE]
      } else {
        dat_fit <- data[data$cv_fold != 1L, , drop = FALSE]
      }
      
      mesh_args[["data"]] <- dat_fit
      mesh <- do.call(make_mesh, mesh_args)
    } else {
      mesh <- spde
      dat_fit <- data
    }
    dot_args <- list(dot_args)[[1]]
    dot_args$offset <- NULL
    .args <- c(list(
      data = dat_fit, formula = formula, time = time, mesh = mesh,
      weights = weights, offset = .offset
    ), dot_args)
    fit1 <- do.call(sdmTMB, .args)
  }
  
  fit_func <- function(k) {
    if (lfo) {
      weights <- ifelse(data$cv_fold <= k, 1, 0)
    } else {
      # data in kth fold get weight of 0:
      weights <- ifelse(data$cv_fold == k, 0, 1)
    }
    
    if (k == 1L && use_initial_fit) {
      object <- fit1
    } else {
      if (!constant_mesh) {
        if (lfo) {
          dat_fit <- data[data$cv_fold <= k, , drop = FALSE]
        } else {
          dat_fit <- data[data$cv_fold != k, , drop = FALSE]
        }
        mesh_args[["data"]] <- dat_fit
        mesh <- do.call(make_mesh, mesh_args)
      } else {
        mesh <- spde
        dat_fit <- data
      }
      dot_args <- as.list(substitute(list(...)))[-1L] # re-evaluate here! issue #54
      dot_args <- list(...)
      dot_args$offset <- NULL
      args <- c(list(
        data = dat_fit, formula = formula, time = time, mesh = mesh, offset = .offset,
        weights = weights, previous_fit = if (use_initial_fit) fit1 else NULL
      ), dot_args)
      object <- do.call(sdmTMB, args)
    }
    
    if (lfo) {
      cv_data <- data[data$cv_fold == (k + lfo_forecast), , drop = FALSE]
    } else {
      cv_data <- data[data$cv_fold == k, , drop = FALSE]
    }
    
    # FIXME: only use TMB report() below to be faster!
    # predict for withheld data:
    predicted <- predict(object, newdata = cv_data, type = "response",
                         offset = if (!is.null(.offset)) cv_data[[.offset]] else rep(0, nrow(cv_data)))
    
    cv_data$cv_predicted <- predicted$est
    response <- get_response(object$formula[[1]])
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted
    
    # FIXME: get LFO working with the TMB report() option below!
    # calculate log likelihood for each withheld observation:
    # trickery to get the log likelihood of the withheld data directly
    # from the TMB report():
    if (!lfo) {
      tmb_data <- object$tmb_data
      tmb_data$weights_i <- ifelse(tmb_data$weights_i == 1, 0, 1) # reversed
      new_tmb_obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = get_pars(object),
        map = object$tmb_map,
        random = object$tmb_random,
        DLL = "sdmTMB",
        silent = TRUE
      )
      lp <- object$tmb_obj$env$last.par.best
      r <- new_tmb_obj$report(lp)
      cv_loglik <- -1 * r$jnll_obs
      cv_data$cv_loglik <- cv_loglik[tmb_data$weights_i == 1]
    } else { # old method; doesn't work with delta models!
      cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    }
    
    ## test
    # x2 <- ll_sdmTMB(object, withheld_y, withheld_mu)
    # identical(round(cv_data$cv_loglik, 6), round(x2, 6))
    # cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    
    list(
      data = cv_data,
      model = object,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  }
  
  if (requireNamespace("future.apply", quietly = TRUE) && parallel) {
    message(
      "Running fits with `future.apply()`.\n",
      "Set a parallel `future::plan()` to use parallel processing."
    )
    if (!is.null(future_globals)) {
      fg <- structure(TRUE, add = future_globals)
    } else {
      fg <- TRUE
    }
    if (lfo) {
      out <- future.apply::future_lapply(seq_len(lfo_validations), fit_func, future.seed = TRUE, future.globals = fg)
    } else {
      out <- future.apply::future_lapply(seq_len(k_folds), fit_func, future.seed = TRUE, future.globals = fg)
    }
  } else {
    message(
      "Running fits sequentially.\n",
      "Install the future and future.apply packages,\n",
      "set a parallel `future::plan()`, and set `parallel = TRUE` to use parallel processing."
    )
    if (lfo) {
      out <- lapply(seq_len(lfo_validations), fit_func)
    } else {
      out <- lapply(seq_len(k_folds), fit_func)
    }
  }
  
  models <- lapply(out, `[[`, "model")
  data <- lapply(out, `[[`, "data")
  fold_cv_ll <- vapply(data, function(.x) sum(.x$cv_loglik), FUN.VALUE = numeric(1L))
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  data[["_sdmTMB_time"]] <- NULL
  row.names(data) <- NULL
  pdHess <- vapply(out, `[[`, "pdHess", FUN.VALUE = logical(1L))
  max_grad <- vapply(out, `[[`, "max_gradient", FUN.VALUE = numeric(1L))
  converged <- all(pdHess)
  out <- list(
    data = data,
    models = models,
    fold_loglik = fold_cv_ll,
    sum_loglik = sum(data$cv_loglik),
    converged = converged,
    pdHess = pdHess,
    max_gradients = max_grad
  )
  `class<-`(out, "sdmTMB_cv")
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

#' @export
#' @import methods
print.sdmTMB_cv <- function(x, ...) {
  nmods <- length(x$models)
  nconverged <- sum(x$pdHess)
  cat(paste0("Cross validation of sdmTMB models with ", nmods, " folds.\n"))
  cat("\n")
  cat("Summary of the first fold model fit:\n")
  cat("\n")
  print(x$models[[1]])
  cat("\n")
  cat("Access the rest of the models in a list element named `models`.\n")
  cat("E.g. `object$models[[2]]` for the 2nd fold model fit.\n")
  cat("\n")
  cat(paste0(nconverged, " out of ", nmods, " models are consistent with convergence.\n"))
  cat("Figure out which folds these are in the `converged` list element.\n")
  cat("\n")
  cat(paste0("Out-of-sample log likelihood for each fold: ", paste(round(x$fold_loglik, 2), collapse = ", "), ".\n"))
  cat("Access these values in the `fold_loglik` list element.\n")
  cat("\n")
  cat("Sum of out-of-sample log likelihoods:", round(x$sum_loglik, 2), "\n")
  cat("More positive values imply better out-of-sample prediction.\n")
  cat("Access this value in the `sum_loglik` list element.\n")
}