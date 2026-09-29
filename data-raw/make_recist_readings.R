## Generates the synthetic reader-study dataset shipped as `recist_readings`.
## Design follows Yoon, Yoon & Hahn (2019, BMC Med Res Methodol 19:90):
## 249 target lesions, 6 readers, 2 replicate sessions at baseline; readings
## in whole millimetres. Error is additive on the square-root scale.
set.seed(2019)
n_les <- 249; n_rd <- 6; n_ses <- 2
truth <- pmax(rlnorm(n_les, meanlog = log(25), sdlog = 0.55), 8)
v_rd <- 0.010; v_int <- 0.050; v_e <- 0.0087          # sqrt-scale variance components
dat <- expand.grid(lesion = seq_len(n_les), reader = seq_len(n_rd), session = seq_len(n_ses))
rd_eff  <- rnorm(n_rd, 0, sqrt(v_rd))
int_eff <- matrix(rnorm(n_les * n_rd, 0, sqrt(v_int)), n_les, n_rd)
dat$diameter <- round((sqrt(truth[dat$lesion]) + rd_eff[dat$reader] +
                       int_eff[cbind(dat$lesion, dat$reader)] +
                       rnorm(nrow(dat), 0, sqrt(v_e)))^2)
dat$reader  <- factor(paste0("R", dat$reader))
dat$session <- as.integer(dat$session)
dat$lesion  <- as.integer(dat$lesion)
recist_readings <- dat[order(dat$lesion, dat$reader, dat$session), c("lesion", "reader", "session", "diameter")]
rownames(recist_readings) <- NULL
save(recist_readings, file = "data/recist_readings.rda", compress = "xz")
