library(tidyverse)
library(Hmsc)
library(pals)

# Set defualt ggplot theme for publication figures
source(
  "https://raw.githubusercontent.com/koundy/ggplot_theme_Publication/master/ggplot_theme_Publication-2.R"
)
source("code/00_Functions.R")

theme_set(theme_Publication(base_size = 8))

model = readRDS("results/fit/model_chains_4_samples_1000_thin_1000.rds")

m = model

postBeta = getPostEstimate(m, parName = "Beta")

spNames = m$spNames %>% str_to_sentence()
# covNames<-m$covNames
covNames = c(
  "Intercept",
  "Width",
  "Watertemp",
  "Turbidity",
  "Salinity",
  "pH",
  bquote(Height[Bank]),
  bquote(Height[Emerg]),
  bquote(Cover[Float]),
  bquote(Cover[Bank]),
  bquote(Cover[Emerg]),
  "Shaded",
  bquote(DO[2])
)

mbeta = postBeta$mean
betaP = postBeta$support

supportLevel = 0.95

toPlot = mbeta
toSig = toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
rownames(toSig) = covNames
colnames(toSig) = spNames
toSig =
  ifelse(toSig != 0, "*", "") %>%
  as.data.frame() %>%
  slice(-1) %>%
  rownames_to_column(var = "coef") %>%
  pivot_longer(-coef, names_to = "species", values_to = "sig")

# replace zeros with NA
toPlot[toPlot == 0] <- NA

# format the data as matrix and add column and row names
betaMat = matrix(toPlot, nrow = m$nc, ncol = ncol(m$Y))
colnames(betaMat) <- spNames
rownames(betaMat) <- covNames
betaMat = betaMat[-1, ]
betaMat = betaMat[, rev(order(colnames(betaMat)))]
betaMat = betaMat %>%
  as.data.frame() %>%
  rownames_to_column(var = "coef") %>%
  pivot_longer(-coef, names_to = "species") %>%
  left_join(toSig)

betaMat = betaMat %>% mutate(taxa = ifelse(species %in% spNames[1:4], "Mosquito", "Predator"))

# Manipulate Species Names
betaMat = betaMat %>%
  mutate(species = ifelse(
    species %in% spNames[1:4],
    str_replace(species, "_", ". "),
    str_replace_all(species, "_", " ")
  ))

#  Plot -------------------------------------------------------------------

cols_to_use = pals::brewer.puor(3)
cols_to_use = c("#fe7d33", "white", "#48227c")

beta = ggplot(betaMat,
              aes(x = coef, y = species, fill = value)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = sig), cex = 3, nudge_y = -0.1) +
  # scale_fill_gradientn(colours = pals::brewer.rdylbu(3)) +
  scale_fill_steps2(
    # n.breaks = 7,
    n.breaks = 7,
    high = cols_to_use[1],
    mid = 'white',
    # mid = cols_to_use[2],
    low = cols_to_use[3],
    nice.breaks = TRUE,
    na.value = "white",
  ) +
  theme(
    # plot.background = element_rect(fill = "white"),
    # legend.background = element_rect(fill = "white"),
    # panel.border = element_rect(fill = NA, color = NA),
    # legend.margin = margin(l = 1, unit = 'cm'),
    # legend.title = element_text(vjust = 0.1),
    # legend.key.width = unit(1, 'cm'),
    # legend.key.height = unit(2.25, 'cm'),
    # legend.position = "top",
    legend.text = element_text(vjust = 0.5, size = 5, angle = 90),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  facet_grid(taxa ~ ., scales = "free", space = "free") +
  scale_x_discrete(expand = c(0, 0), labels = ~ parse(text = .x)) +
  scale_y_discrete(expand = c(0, 0)) +
  # labs(x = bquote(beta ~ paramaters), y = "Species", fill = "Std. Eff", caption = "Support = 90% CI")
  labs(x = bquote(beta ~ -paramaters),
       y = "Species",
       fill = "Std. Eff")
beta
# ggsave(beta, "../Desktop/beta.pdf", width = 85, height = 85, dpi = 600, units = "mm")
ggsave(
  plot = beta,
  "results/plots/beta_90.png",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)
ggsave(
  "results/plots/beta_90.pdf",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)

# -------------------------------------------------------------------------

# Omega Plot
OmegaCor = computeAssociations(model)

## These are the rhyne level residual correlations
toplotrhyne =
  ((OmegaCor[[1]]$support > supportLevel) +
     (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean


colnames(toplotrhyne) = str_to_sentence(colnames(toplotrhyne))
colnames(toplotrhyne) = ifelse(
  colnames(toplotrhyne) %in% spNames[1:4],
  str_replace(colnames(toplotrhyne), "_", ". "),
  str_replace_all(colnames(toplotrhyne), "_", " ")
)

rownames(toplotrhyne) = str_to_sentence(rownames(toplotrhyne))
rownames(toplotrhyne) = ifelse(
  rownames(toplotrhyne) %in% spNames[1:4],
  str_replace(rownames(toplotrhyne), "_", ". "),
  str_replace_all(rownames(toplotrhyne), "_", " ")
)


omega = ggcorrplot::ggcorrplot(toplotrhyne,
                               colors = c(cols_to_use[3], "white", cols_to_use[1])) +
  # theme_light(base_size = 15) +
  theme_light(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = NA),
    # legend.margin = margin(l = 1, unit = 'cm'),
    legend.title = element_text(hjust = 0.08),
    # legend.key.width = unit(1, 'cm'),
    # legend.key.height = unit(2.25, 'cm'),
    # legend.text = element_text(size = 8),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  # bin into + and -
  labs(x = NULL, y = NULL)
# Plot
omega
# ggsave("../Desktop/omega.pdf", width = 85, height = 85, dpi = 600, units = "mm")
ggsave(
  "results/plots/omega.png",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)
ggsave(
  "results/plots/omega_90.pdf",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)

# plotgrid = cowplot::plot_grid(beta, omega, labels = "AUTO")
# ggsave("../Desktop/plotgrid.pdf", width = 170, height = 85, dpi = 600, units = "mm")

omega2 = toplotrhyne |>
  reshape::melt() |>
  mutate(
    value = ifelse(value == 0, NA, value),
    sign = sign(value),
    correlation = ifelse(sign == 1, "+", "-")
  ) |>
  ggplot(aes(x = X1, y = X2, fill = correlation)) +
  geom_tile(col = "white") +
  scale_fill_manual(values = c(cols_to_use[3], cols_to_use[1]),
                    na.translate = FALSE) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.95,
      vjust = 0.2,
      face = "italic"
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 0.95,
      vjust = 0.2,
      face = "italic"
    )
  ) +
  coord_fixed() +
  labs(fill = "Correlation",
       x = NULL,
       y = NULL)

omega2

ggsave(
  "results/plots/omega2.png",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)
ggsave(
  "results/plots/omega2_90.pdf",
  width = 85,
  height = 85,
  dpi = 600,
  units = "mm"
)


# VP ----------------------------------------------------------------------

# Variance Partitioning B
preds = computePredictedValues(model)
MF = evaluateModelFit(hM = model, predY = preds)

# GGPLOT DATA PREP
vp_to_gg = function(VP) {
  rownames = rownames(VP$vals)
  
  tibble = data.frame(VP$vals, Group = rownames)
  
  df = pivot_longer(
    tibble,
    cols = -Group,
    # names_repair = "minimal",
    names_to = "Species",
    values_to = "Variance"
  )
  return(df)
}

# group = c(3, 1, 1, 1, 1, 2, 2, 2, 3, 1)
group = c(3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 1)

groupnames = c("Chemical", "Vegetation", "Structural")

VP = computeVariancePartitioning(model, group = group, groupnames = groupnames)
vals = VP$vals
rownames(vals)

# table
VP$vals |>
  t() |> as.data.frame() |>
  rownames_to_column(var = "Taxa") |>
  select(Taxa, everything()) |>
  gt::gt() |>
  gt::fmt_number() |>
  gt::gtsave("results/tables/VP_SuppTable.docx")

foo = VP$vals |> t() |> as.data.frame() |>
  rowwise() |>
  mutate(RE = sum(`Random: year`,  `Random: season`, `Random: plot_id`))

sd(foo$RE)

### Total variance explained assign to R2
R2 = NULL
if (!is.null(MF$TjurR2)) {
  TjurR2 = MF$TjurR2
  vals = rbind(vals, TjurR2)
  R2 = TjurR2
}
if (!is.null(MF$SR2)) {
  R2 = MF$SR2
  vals = rbind(vals, R2)
}

if (!is.null(R2)) {
  VPr = VP
  for (k in 1:model$ns) {
    VPr$vals[, k] = R2[k] * VPr$vals[, k]
  }
  
  VPr$vals = VPr$vals[, order(-R2)]
}

VPr$vals |>
  t() |> as.data.frame() |>
  rownames_to_column(var = "Taxa") |>
  select(Taxa, everything()) |>
  arrange(Taxa) |>
  gt::gt() |>
  gt::fmt_number() |>
  gt::gtsave("VPr_SuppTable.docx")

VPr$vals |>
  t() |>
  as.data.frame() %>%
  rowSums()

lvl = vp_to_gg(VP) %>%
  mutate(Group =  fct_reorder(Group, Variance, .fun = 'median')) %>%
  pull(Group) %>%
  levels()

spplvl = levels(with(vp_to_gg(VPr), reorder(Species,-Variance)))

tmp_df = vp_to_gg(VPr) %>%
  mutate(Group = factor(Group, levels = lvl)) %>%
  mutate(Group = case_when(# Detect "Radom:"
    str_detect(Group, "Random") ~ "Random",
    .default = Group)) %>%
  mutate(
    Species = Species |> str_to_sentence() |> str_replace_all("_", " "),
    Species = case_when(
      Species == "Cs annulata" ~ "Cs. annulata",
      Species == "Cx pipiens" ~ "Cx. pipiens",
      Species == "An claviger" ~ "An. claviger",
      Species == "An maculipennis" ~ "An. maculipennis",
      .default = Species
    ),
    Taxa =  case_when(
      str_detect(Species, "annulata|pipiens|clav|maculi") ~ "Mosquito",
      .default = "Predator"
    )
  )


mypal = c("#c1272d", "#0000a7", "#eecc16", "#008176", "#b3b3b3")

F2 = tmp_df %>%
  mutate(Group = fct_relevel(Group, "Random"),) %>%
  ggplot(aes(
    x = reorder(Species, -Variance),
    y = Variance,
    fill = Group,
    col = Group
  )) +
  geom_col() +
  scale_fill_manual(values = mypal) +
  scale_color_manual(values = mypal) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 0.95,
    vjust = 0.2,
    face = "italic"
  )) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  facet_grid(~ Taxa, scales = "free_x", space = "free") +
  labs(y = "Total Variance Explained", x = "Species")
F2
ggsave(
  plot = F2,
  "results/plots/VP.png",
  width = 175,
  height = 125,
  units = "mm",
  dpi = 1200
)
ggsave(
  plot = F2,
  "results/plots/VP.pdf",
  width = 175,
  height = 125,
  units = "mm",
  dpi = 1200
)

# Table for VP
library(gtsummary)
VPr_todocx = VPr$vals |> t() |> as.data.frame()
VPr_todocx$Total = rowSums(VPr_todocx)
VPr_todocx = VPr_todocx |> rownames_to_column("Species")
gt::gt(VPr_todocx) |>
  gt::fmt_number() |>
  gt::gtsave("results/tables/Supplmentary_VP.docx")

gt::gt(VPr_todocx) |>
  gt::fmt_number() |>
  gt::gtsave("results/tables/Supplmentary_VP.tex")

# CV ----------------------------------------------------------------------

ModelCV = read_rds("results/fit/ModelCV.rds")
ModelCV_df = bind_rows(ModelCV, .id = "Type")
ModelCV_df$Species = rep(spNames, 2)
ModelCV_df |> write_csv("results/tables/CV.csv")


# ESS ---------------------------------------------------------------------
# Assuming you have already fitted the HMSC model and stored it in the variable 'm'
convergence = ess(model)

pdf("results/plots/ESS_PSRF.pdf")
ess(model)
dev.off()


