# 1) Multiscale anomaly
Zr <- simulate_multiscale(rainsum[[1]], range_km_LF=200, range_km_HF=60, seed=101)

# 2) Back-transform to rainfall units with per-cell mean/SD
mu_r <- app(rainsum, mean, na.rm=TRUE)
sd_r <- app(rainsum, sd, na.rm=TRUE)
rain_future <- mu_r + sd_r * Zr

# 3) Plot
par(mfrow=c(1,2))
plot(rainsum[[1]], main="Historical (example year)")
plot(rain_future, main="Simulated future (multiscale + climatology)")