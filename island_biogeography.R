## Bird Communities in Patchy Grasslands: An Island Biogeography Perspective 
##
## Data source: "Influence of habitat structure and spatial configuration on bird abundance in wet grassland patches" 
#(Silva et al. 2026, Biodiversity and Conservation 519: 111674)
#https://zenodo.org/records/22853918

## data structure:
## 60 sites (Lines) x 138 variables (Columns; environmental variables and bird species)
## 
## Bird counts (114 species) in the patch (number of individuals per species per patch)
## Independent variables:
##  Area       focal patch area (m2)                  -> ISLAND SIZE
##  Proximity  inverse of sum(neighbor area / squared distance) in the buffer
##             HIGHER = MORE ISOLATED (fewer, smaller or more distant neighbors)  -> ISOLATION
##  Distance   edge-to-edge distance to nearest neighboring patch (m)     -> ISOLATION (simple version)
##  Patch shape      shape index (higher = more elongated / complex)
##  Bare, Water, Tree, Bog, Short, Tall   cover inside the patch (m2)      -> HABITAT STRUCTURE
##  Wet, Dry, Natural.forest, Agriculture, Afforestation   area in the buffer (m2)   -> LANDSCAPE
###################################### 



rm(list = ls()) # clear the workspace
graphics.off() # clear the graphics
options(scipen = 999) # avoid  axis labels with scientific notation

## ---- 1. LOADING DATA -----------------------------------------------------------
birds <- read.csv("dataset Influence of habitat structure and spatial configuration on bird abundance in wet grassland patches.csv")

birds <- birds[!is.na(birds$Area), ]              # drop the empty last row
names(birds)                                      # check the column names

# Creating new objects for analysis
species <- birds[, (which(names(birds) == "Tall") + 1):ncol(birds)]   # species columns
pa      <- species > 0                            # presence/absence (patches x species)
birds$S <- rowSums(pa)                            # species richness per patch
birds$N <- rowSums(species)                       # total abundance per patch

## ---- 2. AREA EFFECT: species-area relationship -----------------------------------
## TESTING: do larger patches hold more species?  (H0: z = 0; H1: z > 0)
## Model: log10(S) = log10(c) + z * log10(Area)
sar <- lm(log10(S) ~ log10(Area), data = birds)
summary(sar)                                      # slope = z
confint(sar)

##plot
plot(birds$Area, birds$S, log = "xy", pch = 19, col = "steelblue",
     xlab = "Patch area (m2, log scale)", ylab = "Species richness (log scale)", main = "Area effect")
abline(sar, lwd = 2, col = "firebrick")
legend("topleft", bty = "n", legend = paste("z =", round(coef(sar)[2], 2)))



## ---- 2.1. Compare alternative models of the species-area relationship --------------------------------

#Linear model
linear <- lm(S ~ Area, data = birds)

# Power law:   S = c * Area^z
power <- nls(S ~ c * Area^z, data = birds,
             start = list(c = exp(coef(sar)[1]), z = coef(sar)[2]))


# Asymptotic:  S = a + b / Area   (b < 0: rises and levels off at a)
asymptotic <- nls(S ~ a + b / Area, data = birds,
                  start = list(a = coef(sar)[1], b = coef(sar)[2]))

# ---- Compare models using AICc ----
# install.packages("AICcmodavg") if not installed

library(bbmle)
AICctab( linear, power, asymptotic, 
        nobs = length(birds$Area), weights = TRUE, delta = TRUE, base = TRUE)

#plot
par(las = 1, bty = "l")
plot(birds$Area, birds$S, pch = 16, col = "grey30",
     xlab = "Patch area (m2)", ylab = "Species richness (S)")

# Generate ordered prediction values
area_seq <- seq(
  min(birds$Area),
  max(birds$Area),
  length.out = 200
)

# Linear model
abline(linear, col = "red",lwd = 2)

# Power model
lines(
  area_seq,
  predict(power, newdata = data.frame(Area = area_seq)),
  col = "purple",
  lwd = 2
)

# Asymptotic model
lines(
  area_seq,
  predict(asymptotic, newdata = data.frame(Area = area_seq)),
  col = "blue",
  lwd = 2
)

#legend
legend("topleft", legend = c("Linear", "Power", "Asymptotic"),
       col = c("red", "purple", "blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)

## ---- 3. ISOLATION EFFECT ---------------------------------------------------------------
## TESTING: are more isolated patches (high Proximity index) poorer in species?
## (H0: slope = 0; H1: slope < 0, because HIGH Proximity = ISOLATED)
iso <- lm(log10(S) ~ log10(Proximity), data = birds)
summary(iso)

##Plot
plot(birds$Proximity, birds$S, log = "xy", pch = 19, col = "darkgreen",
     xlab = "Proximity index (log scale; high = isolated)", ylab = "Species richness (log scale)",
     main = "Isolation effect")
abline(iso, lwd = 2, col = "firebrick")


## ---- 4. COMBINED EFFECT OF AREA AND ISOLATION ---------------------------------------------
## TESTING: do area and isolation both matter when analysed together?
## Richness is a count -> GLM (quasipoisson allows for overdispersion).
combined.effect <- glm(S ~ log10(Area) + log10(Proximity), data = birds, family = quasipoisson)
summary(combined.effect)

## TESTING: does isolation matter more in small patches? (interaction area x isolation)
inter <- glm(S ~ log10(Area) * log10(Proximity), data = birds, family = quasipoisson)
anova(combined.effect, inter, test = "F")                    # p < 0.05 = interaction present

##PLOT 

##Create some values for prediction
xA <- 10^seq(log10(min(birds$Area)),      log10(max(birds$Area)),      length.out = 100)
xP <- 10^seq(log10(min(birds$Proximity)), log10(max(birds$Proximity)), length.out = 100)


#plotting

iso_col <- colorRampPalette(c("gold", "darkred"))(100)[cut(log10(birds$Proximity), 100, labels = FALSE)]
plot(birds$Area, birds$S, log = "x", pch = 19, col = iso_col,
     xlab = "Patch area (m2, log scale)", ylab = "Species richness", main = "Comparing isolation levels")
connected <- unname(quantile(birds$Proximity, 0.1)); isolated <- unname(quantile(birds$Proximity, 0.9))
lines(xA, predict(combined.effect, data.frame(Area = xA, Proximity = connected), type = "response"), lwd = 2, col = "gold3")
lines(xA, predict(combined.effect, data.frame(Area = xA, Proximity = isolated),  type = "response"), lwd = 2, col = "darkred")
legend("topleft", bty = "n", lwd = 2, col = c("gold3", "darkred"), legend = c("well connected", "isolated"))


## ---- 6. PASSIVE SAMPLING HYPOTHESIS ---------------------
## TESTING: big patches may have more species simply because they hold more birds.
## If so, richness depends on number of individuals, not on area itself.
summary(lm(log10(N) ~ log10(Area), data = birds))                              # more area -> more birds?
summary(glm(S ~ log10(N), data = birds, family = quasipoisson))                # more birds -> more species?
summary(glm(S ~ log10(Area) + log10(N), data = birds, family = quasipoisson))  # area effect left once N is known?

par(mfrow = c(1, 2))
fitN <- lm(log10(N) ~ log10(Area), data = birds)
plot(birds$Area, birds$N, log = "xy", pch = 19, col = "orange3", xlab = "Patch area (m2)",
     ylab = "Total number of birds", main = "Bigger patch, more birds")
lines(xA, 10^predict(fitN, data.frame(Area = xA)), lwd = 2)
fitSN <- glm(S ~ log10(N), data = birds, family = quasipoisson)
xN <- 10^seq(log10(min(birds$N)), log10(max(birds$N)), length.out = 100)
plot(birds$N, birds$S, log = "x", pch = 19, col = "purple4", xlab = "Total number of birds",
     ylab = "Species richness", main = "More birds, more species")
lines(xN, predict(fitSN, data.frame(N = xN), type = "response"), lwd = 2)

## ---- 7. HABITAT STRUCTURE INSIDE THE PATCH ----------------------------------------------
## TESTING: does habitat variety (bare, water, trees, bog, short grass, tall grass) explain
## richness in addition to area?  Habitat diversity = Shannon index of the cover classes.
cover <- birds[, c("Bare", "Water", "Tree", "Bog", "Short", "Tall")]
birds$HabDiv <- diversity(cover, index = "shannon")

#creating our own shannon-wienner function
birds$HabDiv <- apply(cover, 1, function(x) { p <- x[x > 0] / sum(x); -sum(p * log(p)) })

summary(lm(HabDiv ~ log10(Area), data = birds))                                      # big patches more varied?
summary(glm(S ~ log10(Area) + HabDiv, data = birds, family = quasipoisson))          # both matter?

par(mfrow = c(1, 2))
fitH <- lm(HabDiv ~ log10(Area), data = birds)
plot(birds$Area, birds$HabDiv, log = "x", pch = 19, col = "seagreen4", xlab = "Patch area (m2)",
     ylab = "Habitat diversity", main = "(a) Are big patches more varied?")
lines(xA, predict(fitH, data.frame(Area = xA)), lwd = 2)
plot(birds$HabDiv, birds$S, pch = 19, col = "seagreen4", xlab = "Habitat diversity",
     ylab = "Species richness", main = "(b) Varied patches, more species")
abline(lm(S ~ HabDiv, birds), lwd = 2)

## ---- 8. ISLAND vs HABITAT vs LANDSCAPE: what matters most? ------------------------------------
## TESTING: effects of area, isolation, shape, habitat diversity and the surrounding
## landscape (forest, agriculture, silviculture in the buffer) on richness.
## scale() puts all predictors on the same scale (effect per 1 SD), so they can be compared.
full <- glm(S ~ scale(log10(Area)) + scale(log10(Proximity)) + scale(Shape) + scale(HabDiv) +
              scale(Natural.forest) + scale(Agriculture) + scale(Silviculture),
            data = birds, family = quasipoisson)
summary(full)

est <- coef(full)[-1]; ci <- confint.default(full)[-1, ]
labs <- c("Area", "Proximity (isolation)", "Shape", "Habitat diversity", "Forest", "Agriculture", "Silviculture")
par(mfrow = c(1, 1), mar = c(4.5, 11, 3, 1))
plot(est, 7:1, xlim = range(ci), pch = 19, yaxt = "n", ylab = "", cex = 1.5,
     xlab = "Effect on log richness per 1 SD", main = "What explains richness?")
segments(ci[, 1], 7:1, ci[, 2], 7:1, lwd = 2); abline(v = 0, lty = 2)
axis(2, 7:1, labs, las = 1)


## ---- 9. SINGLE LARGE OR SEVERAL SMALL? (the SLOSS debate) --------------------------------------------
## TESTING: for the same total area, is it better to protect the largest patches first,
## or many small ones? We add patches in different orders and count the species gained.
## ---- SLOSS: single large or several small? ----------------------------------------
## Question: for the same total area, do the largest patches or the smallest patches
## contain more species when added one by one?

# Vector with the area (m2) of each patch
area <- birds$Area

# Function that counts how many species have been found after adding patches one by one.
# 'ord' = the order in which the patches are added (row numbers of the matrix 'pa')
accum <- function(ord) {
  pa_ordered <- pa[ord, ]                      # presence/absence matrix with patches in that order
  running    <- apply(pa_ordered, 2, cumsum)   # running total per species: how many patches so far hold it
  found      <- running > 0                    # TRUE once a species has appeared in at least one patch
  rowSums(found)                               # number of species found after 1, 2, 3, ... patches
}

# The two deliberate orders
small_first <- order(area)                     # patch numbers sorted from smallest to largest
large_first <- order(-area)                    # patch numbers sorted from largest to smallest

# Empty plot with limits that fit all curves:
# x = from 0 to total area of all patches, y = from 0 to total number of species
par(mar = c(4.5, 4.5, 3, 1))
plot(NA, xlim = c(0, sum(area)), ylim = c(0, ncol(pa)),
     xlab = "Cumulative area (m2)", ylab = "Cumulative number of species",
     main = "Single large vs several small")

# Reference: 100 random orders (grey lines) = what you get without any strategy
for (i in 1:100) {
  r <- sample(nrow(pa))                        # random order of the patches
  lines(cumsum(area[r]), accum(r), col = "grey85")   # x = area added so far, y = species found so far
}

# Smallest patches first (blue) and largest patches first (red)
lines(cumsum(area[small_first]), accum(small_first), lwd = 3, col = "steelblue")
lines(cumsum(area[large_first]), accum(large_first), lwd = 3, col = "firebrick")

# Legend
legend("bottomright", bty = "n", lwd = c(3, 3, 1), col = c("firebrick", "steelblue", "grey70"),
       legend = c("largest patches first", "smallest patches first", "random order"))

## How to read it: all lines end at the same point (all patches, all species).
## Blue above red  -> several small patches accumulate species faster per m2.
## Red above blue  -> the largest patches are more efficient.
## Grey cloud      -> variation expected by chance; a thick line above it beats random order.
## If the blue curve is ABOVE the red one, several small patches hold more species per m2.


## ---- Beta-diversity Analysis --------------------------------------------------------------------------------


## ---- 10. NESTEDNESS: are poor patches subsets of rich patches? ------------------------------------------
## Using function NODF from package VEGAN (nestedness metric based on overlap and decreasing fill)
require(vegan)
nest=nestednodf(species, order = TRUE)
nest
par(mar = c(1, 3, 10, 1))   # bottom, left, top, right (in lines)
plot(nest,
     col   = "steelblue4",  # color of occupied cells
     names = TRUE,          # show site/species labels
     main  = "",
     cex.axis = 0.7,
     las = 2)


## ---- 11. DO NEARBY PATCHES SHARE MORE SPECIES? (distance decay + ordination) -----------------------------
## TESTING: similarity in species composition vs geographic distance between patches.
km_x <- birds$Lon * cos(mean(birds$Lat) * pi / 180) * 111      # rough conversion to km
km_y <- birds$Lat * 111
geo  <- dist(cbind(km_x, km_y))
jac  <- dist(pa * 1, method = "binary")                        # Jaccard dissimilarity
cor.test(as.vector(geo), as.vector(jac), method = "spearman", exact = FALSE)

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), las = 1)
plot(as.vector(geo), 1 - as.vector(jac), pch = 19, col = rgb(0, 0, 0, 0.2),
     xlab = "Distance between patches (km)", ylab = "Shared species (Jaccard similarity)",
     main = "(a) Distance decay")
lines(lowess(as.vector(geo), 1 - as.vector(jac)), lwd = 3, col = "firebrick")

pcoa <- cmdscale(jac, k = 2)                                   # ordination of communities
pcoa
cols <- colorRampPalette(c("gold", "darkgreen"))(100)[cut(log10(birds$Area), 100, labels = FALSE)]
plot(pcoa, pch = 19, cex = 1.6, col = cols, xlab = "Axis 1", ylab = "Axis 2",
     main = "(b) Community composition (green = larger patch)")

###Bonus: WHICH SPECIES NEED BIG PATCHES? --------------------------------------------------------------
## TESTING: for each species, does the chance of being present rise with patch area?
## A50 = patch size (m2) at which a species has a 50% chance of being present.
occ <- colSums(pa)
sel <- names(occ)[occ >= 8 & occ <= nrow(pa) - 8]            # skip very rare / very common species
fits <- lapply(sel, function(s) suppressWarnings(glm(pa[, s] ~ log10(birds$Area), family = binomial)))
names(fits) <- sel
b1  <- sapply(fits, function(f) coef(f)[2])
A50 <- 10^sapply(fits, function(f) -coef(f)[1] / coef(f)[2])
sig <- sapply(fits, function(f) summary(f)$coefficients[2, 4]) < 0.05
cat(sum(b1 > 0), "of", length(b1), "species are more likely in larger patches;", sum(sig & b1 > 0), "significantly\n")

ok  <- b1 > 0 & A50 < 10 * max(birds$Area)                   # keep sensible A50 values only
ord <- order(A50[ok]); nm <- gsub("_", " ", names(A50)[ok])[ord]
par(mfrow = c(1, 1), mar = c(4.5, 11, 3, 1))
dotchart(A50[ok][ord], labels = nm, log = "x", pch = 19, cex = 0.7,
         color = ifelse(sig[ok][ord], "firebrick", "grey60"),
         xlab = "Patch area needed for 50% chance of presence (m2, log)",
         main = "Area-demanding species (red = significant)")

