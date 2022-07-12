library(binbayesrct)
library(parallel)
library(tidyverse)

y <- expand_grid(y0 = 0:50, y1 = 0:50)
grid_supr <- function(y, n) {
  apply(y, 1, function(x) beta_prob_supr_approx(x + 1, n - x + 1, reverse = T))
}
y <- cbind(y, p = grid_supr(y, 50))

tibble_res <- function(res, ...) {
  dplyr::bind_rows(parallel::mclapply(res, function(i) {
    lapply(1:length(i),
           function(x) {
             tidyr::gather(tidyr::as_tibble(i[[x]], rownames = "interim"), "arm", !!names(i)[x], -interim)
           }) %>%
      purrr::reduce(dplyr::full_join, by = c("interim", "arm")) %>%
      dplyr::mutate(arm = forcats::fct_inorder(arm)) %>%
      dplyr::arrange(interim, arm)}, ...), .id = "trial") %>%
    dplyr::mutate(trial = as.numeric(trial), interim = as.numeric(interim)) %>%
    dplyr::arrange(trial, interim)
}

trigger_superior <- function(res) {
  res %>%
    filter(arm == 2, interim > 0) %>%
    select(trial, interim, p_supr) %>%
    mutate(
      `0.950` = p_supr > 0.95,
      `0.960` = p_supr > 0.96,
      `0.970` = p_supr > 0.97,
      `0.980` = p_supr > 0.98,
      `0.990` = p_supr > 0.99,
      `0.995` = p_supr > 0.995
    ) %>%
    group_by(trial) %>%
    summarise_if(is.logical, findfirst) %>%
    gather(x, interim, -trial) %>%
    count(x, interim, .drop = F) %>%
    ungroup() %>%
    complete(x, interim = 1:15, fill = list(n = 0)) %>%
    group_by(x) %>%
    mutate(p = n / sum(n),
           cp = cumsum(p)) %>%
    ungroup() %>%
    mutate(x = as.numeric(x))
}

sims <- 1e4

dat_null <- mclapply(1:sims, function(i) gen_potential_outcomes(3000, c(0.2, 0.2)), mc.cores = 14)
res_null <- mclapply(dat_null, function(x) beta_brar_trial(seq(200, 3000, 200), x, approx = T, mwu = T), mc.cores = 14)
res_null <- tibble_res(res_null, mc.cores = 14)

dat_moderate <- mclapply(1:sims, function(i) gen_potential_outcomes(3000, c(0.2, 0.15)), mc.cores = 14)
res_moderate <- mclapply(dat_moderate, function(x) beta_brar_trial(seq(200, 3000, 200), x, approx = T, mwu = T), mc.cores = 14)
res_moderate <- tibble_res(res_moderate, mc.cores = 14)

dat_large <- mclapply(1:sims, function(i) gen_potential_outcomes(3000, c(0.2, 0.1)), mc.cores = 14)
res_large <- mclapply(dat_large, function(x) beta_brar_trial(seq(200, 3000, 200), x, approx = T, mwu = T), mc.cores = 14)
res_large <- tibble_res(res_large, mc.cores = 14)

p <- matrix(rbeta(2*sims, 1, 1), sims, 2)
dat_prior <- mclapply(1:sims, function(i) gen_potential_outcomes(3000, p[i, ]), mc.cores = 14)
res_prior <- mclapply(dat_prior, function(x) beta_brar_trial(seq(200, 3000, 200), x, approx = T, mwu = T), mc.cores = 14)
res_prior <- tibble_res(res_prior, mc.cores = 14)

trigger_res <- trigger_superior(res_null) %>%
  mutate(Scenario = "Null") %>%
  bind_rows(
    trigger_superior(res_moderate) %>%
      mutate(Scenario = "Moderate")
  ) %>%
  bind_rows(
    trigger_superior(res_large) %>%
      mutate(Scenario = "Large")
  ) %>%
  bind_rows(
    trigger_superior(res_prior) %>%
      mutate(Scenario = "Prior (all)")
  ) %>%
  bind_rows(
    trigger_superior(res_prior %>% filter(trial %in% which(p[, 1] > p[, 2]))) %>%
      mutate(Scenario = "Prior (superior)")
  ) %>%
  bind_rows(
    trigger_superior(res_prior %>% filter(trial %in% which(p[, 1] <= p[, 2]))) %>%
      mutate(Scenario = "Prior (inferior)")
  ) %>%
  mutate(Scenario = factor(Scenario,
                           levels = c("Null", "Moderate", "Large",
                                      "Prior (all)", "Prior (inferior)", "Prior (superior)")))
trigger_res %>%
  filter(!is.na(interim)) %>%
  ggplot(., aes(interim, cp)) +
  facet_wrap( ~ Scenario, scales = "free_y") +
  geom_point(aes(colour = x)) +
  geom_line(aes(group = x, colour = x))
