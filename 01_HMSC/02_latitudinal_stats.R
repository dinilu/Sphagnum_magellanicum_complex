library(tidyverse)
data_cc
cm_cc

data <- cm_cc |>
  bind_cols(data_cc) |>
  group_by(Loc) |>
  summarize(diabolicum = max(diabolicum),
            divinum = max(divinum),
            magniae = max(magniae),
            medium = max(medium),
            sympatry = diabolicum + divinum + magniae + medium,
            latitude = mean(Lat))

sdiab_test <- glm(diabolicum ~ latitude, data = data, family = binomial)
sdivi_test <- glm(divinum ~ latitude, data = data, family = binomial)
smagn_test <- glm(magniae ~ latitude, data = data, family = binomial)
smedi_test <- glm(medium ~ latitude, data = data, family = binomial)
sympa_test <- glm(sympatry ~ latitude, data = data, family = poisson)

sdiab_test |> summary()
sdivi_test |> summary()
smagn_test |> summary()
smedi_test |> summary()
sympa_test |> summary()

library(mgcv)
sympa_gam <- mgcv::gam(sympatry ~ s(latitude, sp = -1), data = data, family = poisson, method="REML")
sympa_gam |> plot()
