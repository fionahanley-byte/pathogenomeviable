---
title: "Lab 13: Environmental Predictors of Biosynthetic Gene Cluster Potential in Soil Microbial Communities"
author: "Fiona Hanley"
format:
  html:
    toc: true
    toc_float: true
    embed-resources: true
    theme: cerulean
execute: 
  warning: false
  message: false
---

## This lab was generated entirely through the use of AI.

# Introduction

Microbial secondary metabolites play a central role in human health and medicine, forming the basis of many antibiotics, anticancer agents, immunosuppressants, and signaling molecules. The genes responsible for producing these compounds are typically organized into **biosynthetic gene clusters (BGCs)**—co-localized sets of genes that collectively encode complex metabolic pathways. While BGCs have been extensively studied in cultured microbes and clinical isolates, far less is known about how biosynthetic potential is distributed across natural environments or how environmental conditions shape microbial investment in secondary metabolism. Soils represent one of the most microbially diverse ecosystems on Earth and are a major reservoir of natural products with pharmaceutical relevance. In this study, we leverage metagenome-assembled genomes (MAGs) from the National Ecological Observatory Network (NEON), combined with associated soil chemistry and environmental metadata, to investigate how environmental factors predict community-level biosynthetic potential. Because direct BGC annotation requires complete genome sequences that are not always available, we employ a multi-proxy framework that integrates taxonomic identity, genome architectural properties, and hybrid metrics to infer biosynthetic capacity across environmental gradients.

# Load packages and data

```{r}
library(tidyverse)
library(lme4)
library(broom.mixed)
library(vegan)

```

```{r}
soilMAGs <- read_csv("NEON_soilMAGs_soilChem.csv")

```

# Renaming columns

```{r}
bgc_data <- soilMAGs %>%
  select(
    site = site,
    site_ID = site_ID,
    subplot = subplot,
    horizon = layer,
    phylum = phylum,
    genus = genus,
    genome_size_bp = total_number_of_bases,
    completeness = bin_completeness,
    contamination = bin_contamination,
    soil_pH = soilInWaterpH,
    soil_pH_CaCl = soilInCaClpH
  )

```

# Proxy 1: taxonomy-based BGC potential

To estimate biosynthetic gene cluster (BGC) potential based on evolutionary lineage, we assigned each MAG an expected BGC richness value using broad, literature-supported trends across bacterial phyla. Lineages known to be enriched in secondary metabolism, such as Actinobacteria and Firmicutes, were assigned higher proxy values, while lineages typically associated with streamlined genomes were assigned lower values. This taxonomy-based proxy captures lineage-level constraints on biosynthetic capacity without requiring direct genome sequence analysis.

```{r}
bgc_data <- bgc_data %>%
  mutate(
    bgc_taxonomy_proxy = case_when(
      phylum == "Actinobacteria" ~ 10,
      phylum == "Firmicutes" ~ 6,
      phylum == "Proteobacteria" ~ 3,
      phylum == "Bacteroidota" ~ 2,
      TRUE ~ 1
    )
  )

```

# Proxy 2: Genome-property-based BGC proxy

Because biosynthetic gene clusters are large, energetically costly genomic features, we used genome size as a proxy for overall metabolic and biosynthetic capacity. Genome completeness was included to account for variation in MAG assembly quality and the likelihood of capturing full metabolic pathways. These variables were standardized and combined into a continuous genome-property–based proxy representing the relative capacity of each MAG to encode complex secondary metabolism independent of taxonomic identity.

```{r}
bgc_data <- bgc_data %>%
  mutate(
    genome_size_mb = genome_size_bp / 1e6,
    bgc_genome_proxy =
      scale(genome_size_mb)[,1] +
      scale(completeness)[,1]
  )

```

# Proxy 3: Hybrid proxy (taxonomy x genome size)

To integrate lineage-level expectations with genome-specific investment, we constructed a hybrid BGC proxy by scaling the taxonomy-based BGC estimates by genome size. This hybrid metric captures within-lineage variation in genome expansion and provides a more nuanced estimate of biosynthetic potential, reflecting both evolutionary history and genome architecture.

```{r}
bgc_data <- bgc_data %>%
  mutate(
    bgc_hybrid_proxy = bgc_taxonomy_proxy * genome_size_mb
  )

bgc_data
```

# Aggregate proxies by environmental sample

Because biosynthetic gene clusters are large, energetically costly genomic features, we used genome size as a proxy for overall metabolic and biosynthetic capacity. Genome completeness was included to account for variation in MAG assembly quality and the likelihood of capturing full metabolic pathways. These variables were standardized and combined into a continuous genome-property–based proxy representing the relative capacity of each MAG to encode complex secondary metabolism independent of taxonomic identity.

```{r}
bgc_env <- bgc_data %>%
  group_by(site, subplot, horizon) %>%
  summarise(
    taxonomic_BGC = sum(bgc_taxonomy_proxy, na.rm = TRUE),
    genome_BGC = sum(bgc_genome_proxy, na.rm = TRUE),
    hybrid_BGC = sum(bgc_hybrid_proxy, na.rm = TRUE),
    mean_pH = mean(soil_pH, na.rm = TRUE),
    .groups = "drop"
  )

bgc_env
```

# Environmental models

To assess environmental predictors of biosynthetic potential, we fit linear mixed-effects models with each BGC proxy as a response variable. Soil pH and soil layer were included as fixed effects, while site was included as a random effect to account for spatial non-independence among samples. Identical model structures were applied across all proxy definitions to allow direct comparison of environmental effect strength and direction.

```{r}
library(lme4)

mod_tax <- lmer(taxonomic_BGC ~ mean_pH + horizon + (1 | site), data = bgc_env)
mod_gen <- lmer(genome_BGC ~ mean_pH + horizon + (1 | site), data = bgc_env)
mod_hyb <- lmer(hybrid_BGC ~ mean_pH + horizon + (1 | site), data = bgc_env)

mod_tax
mod_gen
mod_hyb
```

# Compare results across proxies

Model coefficients from each proxy-specific analysis were extracted and compared to evaluate the consistency of environmental effects across proxy definitions. Concordant effect directions and magnitudes across independent proxies were interpreted as robust signals of environmental control on biosynthetic potential, while discrepancies among proxies provided insight into the sensitivity of inference to underlying biological assumptions.

```{r}
library(broom.mixed)

bind_rows(
  tidy(mod_tax) %>% mutate(proxy = "Taxonomy"),
  tidy(mod_gen) %>% mutate(proxy = "Genome"),
  tidy(mod_hyb) %>% mutate(proxy = "Hybrid")
)

```

# Cross proxy agreement

To quantify agreement among BGC proxy metrics, we calculated correlations between taxonomy-based, genome-based, and hybrid estimates of biosynthetic potential at the environmental sample level. This analysis assesses whether different proxy approaches identify similar biosynthetic “hotspots” and helps distinguish generalizable ecological patterns from proxy-specific artifacts.

```{r}
bgc_env %>%
  select(taxonomic_BGC, genome_BGC, hybrid_BGC) %>%
  cor(use = "pairwise.complete.obs")

```

# Visualization of proxy relationships

Relationships among BGC proxies were visualized using scatterplots to illustrate agreement and divergence across methods. Visualization provides an intuitive assessment of methodological robustness and highlights samples or environments where biosynthetic potential estimates differ depending on proxy definition.

### Scatterplot 1: Taxonomy-based vs genome-based proxy

#### Do BGC-rich lineages also tend to have genome architectures consistent with high biosynthetic capacity?

-   Lineages known to be BGC-rich also tend to have large, complex genomes.

-   Environmental selection is acting on both **lineage composition** and **genome architecture** in the same direction.

```{r}
ggplot(bgc_env, aes(x = taxonomic_BGC, y = genome_BGC)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Taxonomy-based BGC proxy",
    y = "Genome-based BGC proxy",
    title = "Agreement between taxonomy- and genome-based BGC proxies"
  )

```

### Scatterplot 2: Taxonomy-based vs hybrid proxy

#### Within the same evolutionary groups, does genome expansion matter?

-   Genome size does not vary much within lineages

-   Hybrid proxy closely tracks taxonomy alone

```{r}
ggplot(bgc_env, aes(x = taxonomic_BGC, y = hybrid_BGC)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Taxonomy-based BGC proxy",
    y = "Hybrid BGC proxy",
    title = "Taxonomy × genome hybrid BGC potential"
  )

```

### Scatterplot 3: Genome-based vs hybrid proxy

#### Does genome size mean the same thing across all lineages?

-   Genome size contributes similarly to biosynthetic potential across lineages

```{r}
ggplot(bgc_env, aes(x = genome_BGC, y = hybrid_BGC)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Genome-based BGC proxy",
    y = "Hybrid BGC proxy",
    title = "Relationship between genome-based and hybrid BGC proxies"
  )

```

### Scatterplot 4: Soil horizons

#### Do **O vs M horizons** cluster differently?

-   Vertical soil stratification in biosynthetic strategies

```{r}
ggplot(bgc_env, aes(x = taxonomic_BGC, y = genome_BGC, color = horizon)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Taxonomy-based BGC proxy",
    y = "Genome-based BGC proxy",
    color = "Soil layer"
  )

```

# Conclusion

Using a multi-proxy approach, this study demonstrates that environmental conditions are strongly associated with variation in inferred biosynthetic potential across NEON soil microbial communities. Concordance among taxonomy-based, genome-property–based, and hybrid BGC proxies suggests that environmental drivers act consistently across evolutionary and genomic scales, while divergence among proxies highlights the importance of within-lineage genome variation. Soil pH and soil layer emerged as key predictors, indicating that both chemical conditions and vertical stratification structure microbial strategies for secondary metabolism. Although proxy-based, these results provide ecologically grounded evidence that natural soil environments selectively favor microbes with greater biosynthetic capacity, reinforcing soils as critical reservoirs of medically relevant natural products. More broadly, this work illustrates how large-scale environmental genomics datasets can be used to generate novel, testable hypotheses about the ecological controls on natural product diversity, even in the absence of complete genome sequences, and motivates future integration with direct BGC annotation and functional validation.
