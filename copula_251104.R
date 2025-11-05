# 2025-11-04 
# attempt at spatial copula-based stochastic simulation approach


# --- Helpers ----------------------------------------------------------------

# Robust clamp away from {0,1} to avoid ±Inf in qnorm
.clamp01 <- function(u, eps = 1e-6) pmin(pmax(u, eps), 1-eps)

# Build empirical CDFs across *space* for one layer and return uniforms for valid cells
ecdf_uniforms <- function(x) {
  ok <- is.finite(x)
  u  <- rep(NA_real_, length(x))
  Fx <- ecdf(x[ok])
  u[ok] <- Fx(x[ok])
  .clamp01(u)
}

# Normal-score transform of uniforms
u_to_z <- function(u) qnorm(.clamp01(u))

# Fit bivariate copula (collocated) from U(0,1) pairs; returns family, par, par2 + Kendall tau
fit_bicop <- function(u1, u2) {
  ok <- is.finite(u1) & is.finite(u2)
  fit <- VineCopula::BiCopSelect(u1[ok], u2[ok],
                                 familyset = c(1,2,3,4,5,6,7,8,9,10)) # broad set incl. t
  list(family = fit$family, par = fit$par, par2 = fit$par2,
       tau = suppressWarnings(cor(u1[ok], u2[ok], method = "kendall")))
}

# Empirical variogram on Z-scores, fit simple model (exponential default)
fit_vario <- function(z, xy, model = "Exp") {
  ok <- is.finite(z)
  if (sum(ok) < 50) stop("Too few valid locations for variogram.")
  df  <- data.frame(x = xy$x[ok], y = xy$y[ok], z = z[ok])
  sp::coordinates(df) <- ~ x + y
  # empirical variogram (direct)
  vg  <- gstat::variogram(z ~ 1, data = df, cutoff = max(diff(range(df@coords[,1])),
                                                         diff(range(df@coords[,2]))) * 0.7)
  # initial guess
  psill0 <- var(df$z, na.rm = TRUE)
  range0 <- max(vg$dist, na.rm = TRUE) / 3
  nug0   <- max(min(vg$gamma, na.rm = TRUE) * 0.1, 1e-4)
  mod0   <- gstat::vgm(psill = psill0, model = model, range = range0, nugget = nug0)
  fit    <- gstat::fit.variogram(vg, model = mod0, fit.method = 6) # WLS
  fit
}

# Convert a single fitted VineCopula "family+par(+par2)" to an approximate
# collocated Pearson rho
copula_rho0 <- function(fam, par, par2) {
  # For Gaussian/t, par is rho directly (t uses par=rho, par2=df)
  if (fam %in% c(1, 2)) {    # 1=Gaussian, 2=Student t in VineCopula encoding
    return(par)
  }
  # For Archimedean, convert Kendall's tau to rho via an approximation:
  #   First get tau from (fam, par, par2), then rho ≈ sin(pi*tau/2) * 1.097 
  #   (rough heuristic).
  #   Better would be Monte Carlo, but this is fast and serviceable as a 
  #   collocated mapping.
#  tau <- VineCopula::BiCopTau2Par(family = fam, tau = NA) # we can't invert like this; fallback MC
  # Fallback: Monte Carlo approximate rho from the fitted copula
  set.seed(1492)
  sim <- VineCopula::BiCopSim(50000, family = fam, par = par, par2 = par2)
  return(cor(sim[,1], sim[,2]))
}

# Given variogram, simulate an unconditional Gaussian field on this grid (sp::SpatialPixels)
simulate_gaussian_field <- function(vgm_fit, rast_template) {
  library(sp)
  library(gstat)
  
  # grid for prediction
  grd <- as.data.frame(terra::crds(rast_template, df = TRUE))
  colnames(grd) <- c("x","y")
  SP <- SpatialPixels(SpatialPoints(grd))
  
  # minimal non-empty dummy data (one arbitrary point)
  dummy <- data.frame(x = 0, y = 0, z = 0)
  coordinates(dummy) <- ~ x + y
  
  # gstat object
  g <- gstat(id = "z", formula = z ~ 1, data = dummy,
             model = vgm_fit, beta = 0)
  
  # unconditional simulation
  sim <- predict(g, newdata = SP, nsim = 1)
  
  # convert to raster
  r <- rast_template[[1]]
  r[] <- sim$sim1
  r
}

range_km_to_cells <- function(r, range_km) {
  # If CRS is in meters, just divide by pixel size
  if (!is.lonlat(r)) {
    return(max(1, round(range_km*1000 / mean(res(r)))))
  }
  # lon/lat: degrees to km ≈ 111 km * cos(lat) for lon; 111 km for lat
  # use the raster’s median latitude
  ymid <- mean(ext(r)[c(3,4)])
  km_per_deg_y <- 111.32
  km_per_deg_x <- 111.32 * cos(ymid * pi/180)
  # convert the larger of the two pixel spacings to km, then to cells
  px_km <- max(res(r)[1]*km_per_deg_x, res(r)[2]*km_per_deg_y)
  max(1, round(range_km / px_km))
}


# --- Fast unconditional Gaussian field via FFT -------------------------------

library(terra)

simulate_gaussian_field_safe <- function(rast_template,
                                         range = 5,
                                         variance = 1,
                                         mean = 0,
                                         seed = NULL,
                                         mask_na = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Save NA pattern
  na_mask <- is.na(rast_template)
  r <- rast_template[[1]]
  r[] <- rnorm(ncell(r))
  if (mask_na) r[na_mask] <- NA
  r
  # 
  # # Convert range to valid kernel size (cells)
  # nrow_r <- nrow(r)
  # ncol_r <- ncol(r)
  # max_radius <- floor(min(nrow_r, ncol_r) / 4)
  # k <- min(range, max_radius)
  # if (k < 1) k <- 1
  # cat("Effective smoothing radius (cells):", k, "\n")
  # 
  # w <- matrix(1, nrow = 2 * k + 1, ncol = 2 * k + 1)
  # w <- w / sum(w)
  # 
  # # ---- use expand=TRUE so the output covers full extent ----
  # smooth <- focal(r, w = w, fun = "mean", na.policy = "omit", expand = TRUE)
  # 
  # # Normalize and rescale
  # vals <- values(smooth, mat = FALSE)
  # sdv <- sd(vals, na.rm = TRUE)
  # if (is.na(sdv) || sdv == 0) {
  #   warning("Standard deviation of smoothed field is 0; returning NA raster")
  #   r[] <- NA
  #   return(r)
  # }
  # 
  # vals <- (vals - mean(vals, na.rm = TRUE)) / sdv
  # vals <- vals * sqrt(variance) + mean
  # smooth[] <- vals
  # 
  # # Reapply NA mask so extent is original but holes remain NA
  # if (mask_na) smooth[na_mask] <- NA
  # 
  # smooth
}


# # test
# rtest <- rast(ncol = 128, nrow = 128, ext = c(0, 1, 0, 1))
# rtest[] <- 0
# r_sim <- simulate_gaussian_field_safe(rtest, range = 5, variance = 1, seed = 123)
# plot(r_sim)



# Per-cell inverse CDF from climatology (10-year history). Uses quantile interpolation.
# Inputs: matrix nCell x nYear; u-vector length nCell -> returns quantile per row.
inv_ecdf_cellwise <- function(mat, u) {
  stopifnot(nrow(mat) == length(u))
  # Use type=7 quantiles at probs u for each row
  # Apply is memory heavy; do in chunks for big rasters
  n  <- nrow(mat)
  ans <- numeric(n)
  chunk <- 5e5
  idxs <- split(seq_len(n), ceiling(seq_len(n)/chunk))
  for (I in idxs) {
    ans[I] <- vapply(I, function(i) {
      x <- mat[i, ]
      x <- x[is.finite(x)]
      if (length(x) < 3) return(NA_real_)
      stats::quantile(x, probs = u[i], type = 7, names = FALSE, na.rm = TRUE)
    }, numeric(1))
  }
  ans
}

##### --- (A) Fit within-year structure for all years --------------------------------

fits <- vector("list", yrs)
for (t in seq_len(yrs)) {
  cat(sprintf("Fitting year %d/%d ...\n", t, yrs))
  r_t <- values(rainsum[[t]], mat = FALSE)
  p_t <- values(mprices[[t]], mat = FALSE)
  ok  <- is.finite(r_t) & is.finite(p_t)
  
  # Uniforms across space (margins per year)
  u_r <- ecdf_uniforms(r_t)
  u_p <- ecdf_uniforms(p_t)
  
  # Collocated copula (at same pixel)
  bic <- fit_bicop(u_r, u_p)
  
  # Spatial correlation via variograms on normal scores
  z_r <- u_to_z(u_r)
  z_p <- u_to_z(u_p)
  vgr <- fit_vario(z_r, grid_xy)
  vgp <- fit_vario(z_p, grid_xy)
  
  fits[[t]] <- list(bic = bic, vgr = vgr, vgp = vgp)
}
fits
saveRDS(fits, file = "fits.rds")
# use if necessary to save time: fits <- loadRDS(file = "fits.rds")

# --- Pool parameters across years ---------------------------------------------

# Collocated dependence: choose the most frequent family; average rho0 (or median)
fam_vec <- sapply(fits, \(f) f$bic$family)
fam_star <- as.numeric(names(sort(table(fam_vec), decreasing = TRUE))[1])

rho0_vec <- numeric(yrs)
# this takes several minutes to run:
for (t in seq_len(yrs)) {
  cat(sprintf("Fitting year %d/%d ...\n", t, yrs))
  b <- fits[[t]]$bic
  rho0_vec[t] <- copula_rho0(b$family, b$par, b$par2)
}
rho0_star <- median(rho0_vec, na.rm = TRUE)

# Variograms: take median of key params from fitted models (nugget, range, psill)
extract_vg <- function(vgfit) {
  # Exponential model commonly returns a structure with model "Exp" + a nugget row.
  # Build tidy named list: nugget, range, psill
  nug <- sum(vgfit$psill[vgfit$model == "Nug"])
  exp_row <- vgfit[vgfit$model != "Nug", ][1, ]  # first structural component
  list(nugget = nug,
       range  = exp_row$range,
       psill  = exp_row$psill)
}

pars_r <- do.call(rbind, lapply(fits, \(f) as.data.frame(extract_vg(f$vgr))))
pars_p <- do.call(rbind, lapply(fits, \(f) as.data.frame(extract_vg(f$vgp))))
pars_r; pars_p
# !!!!
# price parameters do not look right: 
# see outlying values in years 6 & 7, and low overall nugget values
plot(mprices[[seq(5,8)]]) # but the input data for these years looks okay
# !!!!

med_r <- apply(pars_r, 2, median, na.rm = TRUE)
med_p <- apply(pars_p, 2, median, na.rm = TRUE)
med_r; med_p

# variograms
vg_r_star <- gstat::vgm(psill = med_r["psill"], model = "Exp",
                        range = med_r["range"], nugget = med_r["nugget"])
vg_p_star <- gstat::vgm(psill = med_p["psill"], model = "Exp",
                        range = med_p["range"], nugget = med_p["nugget"])
vg_r_star; vg_p_star
plot(vg_r_star,1000)
plot(vg_p_star,100)

cat(sprintf("Pooled: fam=%s (VineCopula code), rho0≈%.3f\n", fam_star, rho0_star))
print(vg_r_star); print(vg_p_star)

# --- (B) Build per-cell climatology margins (10-year) -------------------------

# Matrices: nCell x nYear (WARNING: may be memory heavy for very large rasters)
rmat <- as.matrix(values(rainsum))
pmat <- as.matrix(values(mprices))
ok_rows <- rowSums(is.finite(rmat) & is.finite(pmat)) == ncol(rmat)
if (!all(ok_rows)) {
  message(sprintf("Dropping %d cells with incomplete history.", sum(!ok_rows)))
  rmat <- rmat[ok_rows, , drop = FALSE]
  pmat <- pmat[ok_rows, , drop = FALSE]
}
dim(rmat) # 30551    10

# Keep an index to put simulated values back into full rasters
all_cells <- seq_len(ncell(rainsum))
kept_cells <- which(ok_rows)
  length(all_cells)
  length(kept_cells)
  # 30551 out of 46513 kept for analysis

# --- (C) Simulator: one future joint outcome ----------------------------------


simulate_multiscale <- function(r, range_km_LF = 150, range_km_HF = 30,
                                w_LF = 0.8, w_HF = 0.2, seed = 1) {
  kL <- range_km_to_cells(r, range_km_LF)
  kH <- range_km_to_cells(r, range_km_HF)
  ZL <- simulate_gaussian_field_safe(r, range = kL, variance = 1, seed = seed)
  ZH <- simulate_gaussian_field_safe(r, range = kH, variance = 1, seed = seed+1)
  # standardize both (numeric, not across NA)
  zL <- values(ZL); zL <- (zL - mean(zL, na.rm=TRUE))/sd(zL, na.rm=TRUE)
  zH <- values(ZH); zH <- (zH - mean(zH, na.rm=TRUE))/sd(zH, na.rm=TRUE)
  out <- r[[1]]; out[] <- w_LF*zL + w_HF*zH
  out
}

simulate_future_once <- function(rainsum, mprices,
                                 range_rain_LF = 200, range_rain_HF = 60,
                                 range_price_LF = 250, range_price_HF = 80,
                                 rho = 0.4,
                                 seed = 123) {
  n_years <- min(nlyr(rainsum), nlyr(mprices))
  rain_sims <- vector("list", n_years)
  price_sims <- vector("list", n_years)
  
  for (i in seq_len(n_years)) {
    cat("  Year", i, "...\n")
    
    # Simulate smooth rainfall and price fields
    Zr <- simulate_multiscale(rainsum[[i]], 
                              range_km_LF = range_rain_LF, 
                              range_km_HF = range_rain_HF, 
                              seed = seed + i)
    
    Zp <- simulate_multiscale(mprices[[i]], 
                              range_km_LF = range_price_LF, 
                              range_km_HF = range_price_HF, 
                              seed = seed + 1000 + i)
    
    # Impose cross-field correlation rho
    zr <- values(Zr)
    zp <- values(Zp)
    ok <- is.finite(zr) & is.finite(zp)
    zp[ok] <- rho * zr[ok] + sqrt(1 - rho^2) * zp[ok]
    
    Zp[] <- zp
    rain_sims[[i]]  <- Zr
    price_sims[[i]] <- Zp
  }
  
  list(
    rain_sim = rast(rain_sims),
    price_sim = rast(price_sims)
  )
}



# --- (D) Run simulations ------------------------------------------------------

# result <- simulate_future_once(rainsum, mprices, 
#                                range_rain = 6, range_price = 8,
#                                variance_rain = 1, variance_price = 1,
#                                rho = 0.4, seed = 42)
# plot(result$rain_sim[[1]], main = "Simulated rainfall (year 1)")
# plot(result$price_sim[[1]], main = "Simulated price (year 1)")

# Example: generate 3 stochastic futures
set.seed(42)
nsim <- 3
sim_list <- vector("list", nsim)

for (s in seq_len(nsim)) {
  cat(sprintf("Simulating future scenario %d/%d ...\n", s, nsim))
  sim_list[[s]] <- simulate_future_once(
    rainsum, mprices,
    rho = 0.4,
    seed = 1000 + s
  )
}

# ---- stack all simulated futures ----
rain_future_stack  <- rast(lapply(sim_list, \(x) x$rain_sim))
price_future_stack <- rast(lapply(sim_list, \(x) x$price_sim))

names(rain_future_stack)  <- paste0("rain_future_",  seq_len(dim(rain_future_stack)[3]))
names(price_future_stack) <- paste0("price_future_", seq_len(dim(rain_future_stack)[3]))


# Save (optional)
writeRaster(rain_future_stack,  "rain_future_stack.tif", overwrite=TRUE)
writeRaster(price_future_stack, "price_future_stack.tif", overwrite=TRUE)

# --- (E) Quick diagnostics ----------------------------------------------------

# visualize
plot(mprices[[seq(1,4)]])
plot(price_future_stack[[seq(1,4)]])

# Compute Kendall's tau (rank correlation) across all pixels and all years
sim_u_check <- function(r_sim, p_sim) {
  # Convert to matrices (one column per year)
  rv <- values(r_sim, mat = TRUE)
  pv <- values(p_sim, mat = TRUE)
  
  # Vectorize across years, keeping finite cells only
  ok <- is.finite(rv) & is.finite(pv)
  cor(rank(rv[ok]), rank(pv[ok]), method = "kendall")
}

kendall_s <- vapply(sim_list, \(s)
                    sim_u_check(s$rain_sim, s$price_sim), numeric(1))

cat("Simulated Kendall's tau (collocated), per draw:\n")
print(kendall_s)
