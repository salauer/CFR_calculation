library(tidyverse)
zmeanHDT <- 13
zmedianHDT <- 9.1
muHDT <- log(zmedianHDT)
sigmaHDT <- sqrt(2*(log(zmeanHDT) - muHDT))
cCFRBaseline <- 1.38
cCFREstimateRange <- c(1.23, 1.53)
#cCFRIQRRange <- c(1.3, 1.4)

# Hospitalisation to death distribution
hospitalisation_to_death_truncated <- function(x) {
    dlnorm(x, muHDT, sigmaHDT)
}
# Function to work out correction CFR
scale_cfr <- function(data_1_in, delay_fun){
    case_incidence <- data_1_in$new_cases
    death_incidence <- data_1_in$new_deaths
    cumulative_known_t <- 0 # cumulative cases with known outcome at time tt
    # Sum over cases up to time tt
    for(ii in 1:nrow(data_1_in)){
        known_i <- 0 # number of cases with known outcome at time ii
        for(jj in 0:(ii - 1)){
            known_jj <- (case_incidence[ii - jj]*delay_fun(jj))
            known_i <- known_i + known_jj
        }
        cumulative_known_t <- cumulative_known_t + known_i # Tally cumulative known
    }
    # naive CFR value
    b_tt <- sum(death_incidence)/sum(case_incidence)
    # corrected CFR estimator
    p_tt <- sum(death_incidence)/cumulative_known_t
    data.frame(nCFR = b_tt,
               cCFR = p_tt,
               total_deaths = sum(death_incidence),
               cum_known_t = round(cumulative_known_t),
               total_cases = sum(case_incidence))
}
# Get data
allDat <- NCoVUtils::get_ecdc_cases()
reportDataFinal <- allDat %>%
    dplyr::arrange(country, date) %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>%
    dplyr::select(date, country, new_cases = cases, new_deaths = deaths) %>%
    dplyr::filter(country != "CANADA",
                  country != "Cases_on_an_international_conveyance_Japan") %>%
    dplyr::group_by(country) %>%
    padr::pad() %>%
    dplyr::mutate(new_cases = tidyr::replace_na(new_cases, 0),
                  new_deaths = tidyr::replace_na(new_deaths, 0),
                  cum_deaths = sum(new_deaths)) %>%
    dplyr::filter(cum_deaths > 0) %>%
    dplyr::select(-cum_deaths) %>%
    dplyr::do(scale_cfr(., delay_fun = hospitalisation_to_death_truncated)) %>%
    dplyr::filter(cum_known_t > 0,
                  total_deaths<=total_cases,
                  total_deaths<=cum_known_t) %>%
    dplyr::mutate(nCFR_UQ = binom.test(total_deaths, total_cases)$conf.int[2],
                  nCFR_LQ = binom.test(total_deaths, total_cases)$conf.int[1],
                  cCFR_UQ = binom.test(total_deaths, cum_known_t)$conf.int[2],
                  cCFR_LQ = binom.test(total_deaths, cum_known_t)$conf.int[1],
                  underreporting_estimate = cCFRBaseline / (100*cCFR),
                  lower = cCFREstimateRange[1] / (100 * cCFR_UQ),
                  upper = cCFREstimateRange[2] / (100 * cCFR_LQ),
                  quantile25 = binom.test(total_deaths, cum_known_t, conf.level = 0.5)$conf.int[1],
                  quantile75 = binom.test(total_deaths, cum_known_t, conf.level = 0.5)$conf.int[2]) %>%
    # dplyr::filter(total_deaths > 10) %>%
    dplyr::select(country, total_cases, total_deaths,
                  underreporting_estimate, lower,
                  upper) %>%
    ungroup() %>%
    dplyr::mutate(underreporting_estimate = ifelse(underreporting_estimate <= 1,
                                                   underreporting_estimate, 1) %>%
                      signif(2),
                  upper = ifelse(upper <= 1, upper, 1),
                  lower = signif(lower, 2),
                  upper = signif(upper, 2),
                  country = country %>% stringr::str_replace_all("_", " "),
                  underreporting_estimate_clean = paste0(underreporting_estimate*100,
                                                         "% (",lower*100,"% - ",upper*100,"%)"))

# saveRDS(reportDataFinal, "data/all_together_clean.rds")
