################################################################################
## DESCRIPTIVE ANALYSIS: BIRDS IN WET-GRASSLAND PATCHES  (base R + vegan)
## One row = one patch. Nothing is transformed in the data: logs are applied
## inside plots or models only.
################################################################################
# install.packages("vegan")      # once, if needed
library(vegan)
options(scipen = 999)            # avoid 1e+05 style axis labels

## ---- 0. LOAD DATA AND BUILD BASIC VARIABLES ---------------------------------------
birds <- read.csv("dataset_Influence_of_habitat_structure_and_spatial_configuration_on_bird_abundance_in_wet_grassland_patches.csv",
                  encoding = "latin1")
birds <- birds[!is.na(birds$Area), ]                         # drop the empty last row

species <- birds[, (which(names(birds) == "Tall") + 1):ncol(birds)]   # abundance matrix (patches x species)
pa      <- species > 0                                       # presence/absence matrix
cover   <- birds[, c("Bare", "Water", "Tree", "Bog", "Short", "Tall")] # habitat cover inside the patch (m2)

birds$S       <- rowSums(pa)                                 # species richness per patch
birds$N       <- rowSums(species)                            # total abundance per patch
birds$Shannon <- diversity(species, index = "shannon")       # Shannon diversity of birds
birds$Simpson <- diversity(species, index = "simpson")       # Simpson diversity of birds
birds$Evenness <- birds$Shannon / log(birds$S)               # Pielou's evenness (0-1)
birds$HabDiv  <- diversity(cover, index = "shannon")         # habitat diversity inside the patch

cat("Patches:", nrow(birds), "| Species:", ncol(species), "| Sites:", length(unique(birds$Site)), "\n")
summary(birds[, c("Area", "Proximity", "Distance", "Shape", "S", "N", "HabDiv")])

################################################################################
## PART A. DESCRIBING THE PATCHES
################################################################################

## ---- A1. Histograms of the main variables ---------------------------------------------
## QUESTION: what does the sampling design look like? (area is strongly right-skewed)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), las = 1)
hist(birds$Area, breaks = 15, col = "steelblue", border = "white", main = "Patch area", xlab = "Area (m2)")
hist(log10(birds$Area), breaks = 15, col = "steelblue", border = "white", main = "Patch area (log10)", xlab = "log10 Area")
hist(birds$Proximity, breaks = 15, col = "darkgreen", border = "white", main = "Proximity index", xlab = "Proximity (high = isolated)")
hist(birds$S, breaks = 10, col = "orange3", border = "white", main = "Species richness", xlab = "Species per patch")
hist(birds$N, breaks = 15, col = "orange3", border = "white", main = "Total abundance", xlab = "Birds per patch")
hist(birds$HabDiv, breaks = 10, col = "seagreen4", border = "white", main = "Habitat diversity", xlab = "Shannon index of cover")

## ---- A2. Map of the patches ----------------------------------------------------------------
## QUESTION: how are the patches arranged in space, and where are the rich patches?
rich_col <- colorRampPalette(c("lightyellow", "darkblue"))(100)[cut(birds$S, 100, labels = FALSE)]
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
plot(birds$Lon, birds$Lat, pch = 21, bg = rich_col, col = "grey20",
     cex = 1 + 3 * sqrt(birds$Area / max(birds$Area)),
     asp = 1 / cos(mean(birds$Lat) * pi / 180), xlab = "Longitude", ylab = "Latitude",
     main = "Patches: size = area, colour = richness (dark = more species)")

## ---- A3. Richness and abundance among sites and years ---------------------------------------
## QUESTION: how much do patches differ among sites and sampling years?
par(mfrow = c(1, 2))
boxplot(S ~ factor(Site), data = birds, col = "orange3", xlab = "Site", ylab = "Species richness", main = "Richness by site")
boxplot(S ~ factor(Year), data = birds, col = "orange3", xlab = "Year", ylab = "Species richness", main = "Richness by year")

## ---- A4. Relationships among variables ---------------------------------------------------------
## QUESTION: which variables are related before any modelling? (logs only inside the plot)
pairs(log10(birds[, c("Area", "Proximity", "Distance", "Shape", "S", "N")]),
      pch = 19, col = rgb(0, 0, 0, 0.4), main = "Pairs plot (all variables log10)")

## ---- A5. Habitat composition of each patch -------------------------------------------------------
## QUESTION: does the mix of habitat types change from small to large patches?
prop <- as.matrix(cover / rowSums(cover))                    # proportion of each cover type
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 7), xpd = TRUE)
barplot(t(prop)[, order(birds$Area)], col = hcl.colors(6, "Set 2"), border = NA, names.arg = rep("", nrow(prop)),
        xlab = "Patches (smallest to largest)", ylab = "Proportion of patch area",
        main = "Habitat composition")
legend("right", inset = c(-0.17, 0), legend = colnames(prop), fill = hcl.colors(6, "Set 2"), cex = 0.8, bty = "n")
par(xpd = FALSE, mar = c(4.5, 4.5, 3, 1))

################################################################################
## PART B. DESCRIBING THE BIRD COMMUNITY
################################################################################

## ---- B1. Rank-abundance curve and most abundant species --------------------------------------------
## QUESTION: is the community dominated by a few species, with many rare ones?
tot <- sort(colSums(species), decreasing = TRUE)             # total abundance of each species
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(tot, log = "y", type = "b", pch = 19, col = "steelblue", xlab = "Species rank",
     ylab = "Total abundance (log scale)", main = "Rank-abundance curve")
par(mar = c(4.5, 9, 3, 1))
barplot(rev(head(tot, 10)), horiz = TRUE, las = 1, col = "steelblue", cex.names = 0.7,
        names.arg = gsub("_", " ", rev(names(head(tot, 10)))), xlab = "Total abundance",
        main = "10 most abundant species")

## ---- B2. Occupancy-frequency and abundance-occupancy relationship -------------------------------------
## QUESTION: are most species rare (few patches) or widespread? Are widespread species also abundant?
occ <- colSums(pa)                                           # number of patches occupied by each species
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
hist(occ, breaks = 15, col = "seagreen4", border = "white", xlab = "Number of patches occupied",
     ylab = "Number of species", main = "Occupancy-frequency")
plot(occ, colSums(species) / occ, log = "xy", pch = 19, col = "seagreen4",
     xlab = "Patches occupied (log)", ylab = "Mean abundance where present (log)",
     main = "Abundance-occupancy")
abline(lm(log10(colSums(species) / occ) ~ log10(occ)), lwd = 2)

## ---- B3. Species accumulation curve ---------------------------------------------------------------------
## QUESTION: is the sampling complete? Does the number of species level off as patches are added?
acc <- sapply(1:100, function(i) rowSums(apply(pa[sample(nrow(pa)), ], 2, cumsum) > 0))   # 100 random orders
par(mfrow = c(1, 1))
matplot(acc, type = "l", lty = 1, col = "grey80", xlab = "Number of patches sampled",
        ylab = "Cumulative number of species", main = "Species accumulation curve")
lines(rowMeans(acc), lwd = 3, col = "firebrick")

## ---- B4. Diversity indices per patch ------------------------------------------------------------------------
## QUESTION: how do richness, diversity and evenness vary among patches?
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
hist(birds$Shannon,  breaks = 10, col = "orange3", border = "white", main = "Shannon diversity", xlab = "H'")
hist(birds$Simpson,  breaks = 10, col = "orange3", border = "white", main = "Simpson diversity", xlab = "1 - D")
hist(birds$Evenness, breaks = 10, col = "orange3", border = "white", main = "Evenness", xlab = "Pielou's J")

## ---- B5. Richness vs abundance ----------------------------------------------------------------------------------
## QUESTION: how strongly are number of species and number of birds linked?
par(mfrow = c(1, 1))
plot(birds$N, birds$S, log = "x", pch = 19, col = "purple4", xlab = "Total abundance (log scale)",
     ylab = "Species richness", main = "More birds, more species")
lines(lowess(log10(birds$N), birds$S) |> (\(l) list(x = 10^l$x, y = l$y))(), lwd = 2)

################################################################################
## PART C. BASIC COMMUNITY PATTERNS
################################################################################

## ---- C1. Heatmap of the abundance matrix ------------------------------------------------------------------------------
## QUESTION: what does the whole dataset look like at a glance?
## Patches ordered by area (largest at the top); species ordered by total abundance.
ordP <- order(birds$Area, decreasing = TRUE)
ordS <- order(colSums(species), decreasing = TRUE)
heat <- log10(as.matrix(species)[rev(ordP), ordS] + 1)       # log10(abundance + 1) only for the picture
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
image(1:ncol(heat), 1:nrow(heat), t(heat), col = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
      xlab = "Species (most to least abundant)", ylab = "Patches (small to large)",
      main = "Abundance matrix (log10(N + 1))")

## ---- C2. Cluster dendrogram of patches ---------------------------------------------------------------------------------------
## QUESTION: which patches have similar species lists? Do they group by site?
jac <- vegdist(pa * 1, method = "jaccard", binary = TRUE)    # Jaccard dissimilarity between patches
clus <- hclust(jac, method = "average")
par(mar = c(2, 4.5, 3, 1))
plot(clus, labels = birds$Patch, cex = 0.6, hang = -1, main = "Cluster of patches (Jaccard, average linkage)", xlab = "", sub = "")

## ---- C3. Beta diversity summary -------------------------------------------------------------------------------------------------
## QUESTION: how different are patches from each other, on average?
cat("\nMean Jaccard dissimilarity between patches:", round(mean(jac), 2),
    "| range:", round(min(jac), 2), "-", round(max(jac), 2), "\n")
par(mar = c(4.5, 4.5, 3, 1))
hist(jac, breaks = 20, col = "grey60", border = "white", xlab = "Jaccard dissimilarity (0 = identical, 1 = no shared species)",
     main = "How different are patches?")

## ---- C4. Core and satellite species ---------------------------------------------------------------------------------------------
## QUESTION: how many species are widespread (core) and how many are rare (satellite)?
frac <- occ / nrow(pa)                                       # proportion of patches occupied
cat("Core species (> 80% of patches):", sum(frac > 0.8), "\n")
cat("Satellite species (< 20% of patches):", sum(frac < 0.2), "\n")
cat("Singletons (1 patch only):", sum(occ == 1), "\n")
