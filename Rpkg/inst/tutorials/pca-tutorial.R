library(mbl2018)
library(ggplot2)
theme_set(theme_bw())

mouse_may <- mbl_load_rnaseq("mouse", "may")

# PCA on all the samples
pca_may <- mbl_pca(mouse_may)
ggplot(pca_may$data, aes(PC1, PC2, color = genotype)) +
  geom_point()

ggplot(pca_may$data, aes(PC1, PC2, color = genotype)) +
  geom_point()

ggplot(pca_may$data, aes(PC1, PC2, color = source)) +
  geom_point()

ggplot(pca_may$data, aes(PC1, PC2, color = source, shape = genotype)) +
  geom_point(size = 2)

ggplot(pca_may$data, aes(PC2, PC3, color = source, shape = genotype)) +
  geom_point(size = 2)

# PCA on just the samples from the palate
pca_may <- mbl_pca(mouse_may[, mouse_may$samples$source == "palate"])
ggplot(pca_may$data, aes(PC1, PC2, color = genotype)) +
  geom_point()


mouse_mbl <- mbl_load_rnaseq("mouse", "mbl")
pca_mbl <- mbl_pca(mouse_mbl[, !mouse_mbl$samples$genotype %in% c("Old", "Young")])

ggplot(pca_mbl$data, aes(PC2, PC3, color = genotype)) +
  geom_point()


# Use plotly for interactive ggplot2 graphics
library(plotly)

gg <- ggplot(pca_mbl$data, aes(PC2, PC3, color = genotype,
                               text = paste0("sampleid: ", sample_id))) +
  geom_point()
ggplotly(gg)

# 3D PLOTS!
# Note that plotly uses the pipe "%>%" instead of ggplot2's "+"
plot_ly(pca_mbl$data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~genotype) %>%
  add_markers()

