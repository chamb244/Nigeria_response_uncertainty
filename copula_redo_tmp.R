
# --- Packages ---------------------------------------------------------------
library(terra)
library(VineCopula)
#library(gstat)
#library(sp)
library(stats)


# --- Inputs (provide these) -------------------------------------------------
setwd("C:/Github/chamb244/EiA2030-ex-ante/Nigeria_response_uncertainty")
rainsum <- rast("rainsum.tif")
mprices <- rast("mprices.tif")

stopifnot(nlyr(rainsum) == nlyr(mprices))
stopifnot(all.equal(ext(rainsum), ext(mprices)))
stopifnot(all.equal(res(rainsum), res(mprices)))

# Align masks
m <- !app(rainsum, is.na) & !app(mprices, is.na)
rain  <- mask(rainsum,  m)
price <- mask(mprices, m)

cat("rainsum layers:", nlyr(rainsum), " | mprices layers:", nlyr(mprices), "\n")
cat("rainsum valid share (layer 1):", mean(is.finite(values(rainsum[[1]]))), "\n")
cat("mprices valid share (layer 1):", mean(is.finite(values(mprices[[1]]))), "\n")

# 1) Single-layer field on real data (should be visible)
Zr_test <- simulate_gaussian_field_safe(rainsum[[1]], range = 8, seed = 123)
cat("Zr_test valid share:", mean(is.finite(values(Zr_test))), "\n")
plot(Zr_test, main = "Zr_test (single-layer)")

# 2) Multiscale field (should show large gradients)
Zr_ms <- simulate_multiscale(rainsum[[1]], range_km_LF = 200, range_km_HF = 60, seed = 123)
cat("Zr_ms valid share:", mean(is.finite(values(Zr_ms))), "\n")
plot(Zr_ms, main = "Zr_ms (multiscale)")

# 3) One full future (10 layers)
one <- simulate_future_once(rainsum, mprices, rho = 0.4, seed = 1001)
cat("one$rain_sim layers:", nlyr(one$rain_sim), " valid share (y1):", 
    mean(is.finite(values(one$rain_sim[[1]]))), "\n")
plot(one$rain_sim[[1]], main = "one$rain_sim (year 1)")


# Convenience
yrs <- nlyr(rainsum)
#grid_xy <- as.data.frame(geom(rainsum, cells = TRUE)[, c("x","y","cell")]) # coords (+ cell id)
xy_coords <- xyFromCell(rainsum, 1:ncell(rainsum))
grid_xy <- as.data.frame(cbind(xy_coords, 1:ncell(rainsum)))
colnames(grid_xy) <- c("x","y","cell")
head(grid_xy)


#---------------------------- 
# Helper: convert range in km to cells
#----------------------------
# Note: assumes roughly constant pixel size across the map in lon/lat
# --> may want to redo everything in meter based projection
#     if doing this at SSA level?
range_km_to_cells <- function(r, range_km) {
  if (!is.lonlat(r)) {
    return(max(1, round(range_km * 1000 / mean(res(r)))))
  }
  ymid <- mean(ext(r)[c(3,4)])
  km_per_deg_y <- 111.32
  km_per_deg_x <- 111.32 * cos(ymid * pi/180)
  px_km <- max(res(r)[1] * km_per_deg_x, res(r)[2] * km_per_deg_y)
  max(1, round(range_km / px_km))
}



#----------------------------
# Simulate Gaussian field
#----------------------------

simulate_gaussian_field_safe <- function(rast_template,
                                         range = 5,      # radius in cells
                                         variance = 1,
                                         mean = 0,
                                         seed = NULL,
                                         mask_na = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- 1. Base mask and white noise
  na_mask <- is.na(rast_template)
  M <- rast_template[[1]]
  M[] <- ifelse(is.finite(values(rast_template)), 1, 0)
  
  Z <- rast_template[[1]]
  Z[] <- rnorm(ncell(Z))
  
  # --- 2. Circular kernel
  k <- max(1, round(range))
  d <- 2 * k + 1
  W <- matrix(0, nrow = d, ncol = d)
  for (i in 1:d) for (j in 1:d)
    if ((i - k - 1)^2 + (j - k - 1)^2 <= k^2) W[i, j] <- 1
  W <- W / sum(W)
  
  # --- 3. Edge-corrected smoothing
  # use fill=0 so edge cells get partial support
  num <- focal(Z * M, w = W, fun = "sum", fill = 0, expand = TRUE)
  den <- focal(M,     w = W, fun = "sum", fill = 0, expand = TRUE)
  
  num_vals <- values(num)
  den_vals <- values(den)
  eps <- 1e-6
  
  # divide safely; keep everywhere except full NA regions
  out_vals <- rep(NA_real_, length(num_vals))
  ok <- den_vals > eps
  out_vals[ok] <- num_vals[ok] / pmax(den_vals[ok], eps)
  
  # --- 4. Normalize and rescale
  S <- rast_template[[1]]
  S[] <- out_vals
  
  vals <- values(S)
  if (!any(is.finite(vals))) {
    warning("All values NA; check mask or range.")
    return(S)
  }
  
  mu <- mean(vals, na.rm = TRUE)
  sdv <- sd(vals, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) {
    warning("Zero variance field; returning flat raster.")
    S[] <- mean
  } else {
    vals <- (vals - mu) / sdv
    vals <- vals * sqrt(variance) + mean
    S[] <- vals
  }
  
  # --- 5. Reapply NA mask only at very end
  if (mask_na) S[na_mask] <- NA
  
  return(S)
}



#----------------------------
# Multiscale simulator (with gradients)
#----------------------------
simulate_multiscale <- function(r, range_km_LF = 150, range_km_HF = 30,
                                w_LF = 0.8, w_HF = 0.2, seed = 1) {
  kL <- range_km_to_cells(r, range_km_LF)
  kH <- range_km_to_cells(r, range_km_HF)
  
  ZL <- simulate_gaussian_field_safe(r, range = kL, variance = 1, seed = seed)
  ZH <- simulate_gaussian_field_safe(r, range = kH, variance = 1, seed = seed + 1)
  
  zL <- values(ZL); zL <- (zL - mean(zL, na.rm = TRUE)) / sd(zL, na.rm = TRUE)
  zH <- values(ZH); zH <- (zH - mean(zH, na.rm = TRUE)) / sd(zH, na.rm = TRUE)
  
  out <- r[[1]]
  out[] <- w_LF * zL + w_HF * zH
  out
}


#### Future simulator (multi-year + joint correlation) ####
# --------------------------------------------------------

simulate_future_once <- function(rainsum, mprices,
                                 range_rain_LF = 200, range_rain_HF = 60,
                                 range_price_LF = 250, range_price_HF = 80,
                                 rho = 0.4, seed = 123) {
  
  n_years <- min(nlyr(rainsum), nlyr(mprices))
  rain_sims  <- vector("list", n_years)
  price_sims <- vector("list", n_years)
  
  for (i in seq_len(n_years)) {
    cat(sprintf("  Year %d/%d ...\n", i, n_years))
    
    # Multiscale Gaussian fields (each smooth + realistic)
    Zr <- simulate_multiscale(rainsum[[i]],
                              range_km_LF = range_rain_LF,
                              range_km_HF = range_rain_HF,
                              seed = seed + i)
    
    Zp <- simulate_multiscale(mprices[[i]],
                              range_km_LF = range_price_LF,
                              range_km_HF = range_price_HF,
                              seed = seed + 1000 + i)
    
    # Impose correlation ρ
    zr <- values(Zr)
    zp <- values(Zp)
    ok <- is.finite(zr) & is.finite(zp)
    zp_new <- zp
    zp_new[ok] <- rho * zr[ok] + sqrt(1 - rho^2) * zp[ok]
    Zp[] <- zp_new
    
    rain_sims[[i]]  <- Zr
    price_sims[[i]] <- Zp
  }
  
  list(
    rain_sim  = rast(rain_sims),
    price_sim = rast(price_sims)
  )
}



#----------------------------
# Main multi-simulation driver
#----------------------------
simulate_all_futures <- function(rainsum, mprices, nsim = 3, rho = 0.4,
                                 range_rain_LF = 200, range_rain_HF = 60,
                                 range_price_LF = 250, range_price_HF = 80,
                                 seed_base = 1000) {
  
  template <- rainsum[[1]]  # enforce consistent geometry
  sim_list <- vector("list", nsim)
  
  for (s in seq_len(nsim)) {
    cat(sprintf("\nSimulating future scenario %d/%d ...\n", s, nsim))
    sim_raw <- simulate_future_once(
      rainsum, mprices,
      range_rain_LF  = range_rain_LF,
      range_rain_HF  = range_rain_HF,
      range_price_LF = range_price_LF,
      range_price_HF = range_price_HF,
      rho   = rho,
      seed  = seed_base + s * 100
    )
    
    # 🔧 ensure exact same extent/resolution/mask
    sim_list[[s]] <- list(
      rain_sim  = resample(sim_raw$rain_sim,  template, method = "bilinear"),
      price_sim = resample(sim_raw$price_sim, template, method = "bilinear")
    )
  }
  
  # Stack all aligned layers
  rain_future_stack  <- rast(do.call(c, lapply(sim_list, \(x) x$rain_sim)))
  price_future_stack <- rast(do.call(c, lapply(sim_list, \(x) x$price_sim)))
  
  # Automatic naming: sim#, year#
  n_years <- nlyr(rainsum)
  names(rain_future_stack) <-
    unlist(lapply(seq_len(nsim),
                  \(s) paste0("rain_future_sim", s, "_year", seq_len(n_years))))
  names(price_future_stack) <-
    unlist(lapply(seq_len(nsim),
                  \(s) paste0("price_future_sim", s, "_year", seq_len(n_years))))
  
  list(
    rain_future_stack  = rain_future_stack,
    price_future_stack = price_future_stack,
    sim_list           = sim_list
  )
}



#----------------------------
# Run full workflow
#----------------------------
res <- simulate_all_futures(rainsum, mprices,
                            nsim = 3, rho = 0.4,
                            range_rain_LF = 200, range_rain_HF = 60,
                            range_price_LF = 250, range_price_HF = 80)

# rain_future_stack <- c(sim_list[[1]]$rain_sim, sim_list[[2]]$rain_sim, sim_list[[3]]$rain_sim)
# price_future_stack <- c(sim_list[[1]]$price_sim, sim_list[[2]]$price_sim, sim_list[[3]]$price_sim)
rain_future_stack  <- do.call(c, lapply(sim_list, \(x) x$rain_sim))
price_future_stack  <- do.call(c, lapply(sim_list, \(x) x$price_sim))

#visualize just to check 
plot(rain_future_stack[[c(1,11,21)]])
plot(price_future_stack[[c(1,11,21)]])

cat("Rain stack layers:", nlyr(rain_future_stack),
    " | Price stack layers:", nlyr(price_future_stack), "\n")


# Label layers
n_years <- nlyr(rainsum)
names(rain_future_stack) <-
  unlist(lapply(seq_len(nsim),
                \(s) paste0("rain_future_sim", s, "_year", seq_len(n_years))))
names(price_future_stack) <-
  unlist(lapply(seq_len(nsim),
                \(s) paste0("price_future_sim", s, "_year", seq_len(n_years))))

names(rain_future_stack)
names(price_future_stack)

