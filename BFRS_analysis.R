# BFRS_analysis

# set up ----

## libraries ----
packageVersion("Rmisc")
citation("lme4")

# usage confirmed
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library("phyloseq") # V.1.46.0







## import data ----

### feature table ----
# full otu frequency data (exported FeatureTable[frequency] from qiime2)
# note, otu = asv , they are used interchangeably in naming

OTU_df <- read.csv("otu_table.csv", row.names = 1)
# 6843 obs (OTUs) of 36 vars (samples)

# make phyloseq object
OTU <- otu_table(OTU_df, taxa_are_rows = TRUE)


### sample data ----
sample_df <- read.csv("sample_data.csv")
# 36 obs (samples) of 3 vars 

# set row names & remove unnecessary column
rownames(sample_df) <- sample_df$samplename
sample_df$samplename <- NULL
str(sample_df)

# set all as factors
sample_df$site <- as.factor(sample_df$site)
sample_df$plot <- as.factor(sample_df$plot)
sample_df$treatment <- as.factor(sample_df$treatment)
sample_df$treatment <- factor(sample_df$treatment, levels = c("Control", "Fire", "Mech", "Mech+Fire"))
str(sample_df)

# make phyloseq object
SAMPLE <- phyloseq::sample_data(sample_df)


### tax data ----
tax_df <- read.csv("tax_table.csv")
rownames(tax_df) <- tax_df$id
tax_df$id <- NULL
tax_matrix <- as.matrix(tax_df)

# make phyloseq object
TAX <- tax_table(tax_matrix)


### phylosequ-class object ----
physeq_object <- phyloseq(OTU, SAMPLE, TAX)
# note, leaving out the phylogenetic tree in this version


# data exploration ----
# option to filter out taxa with low abundance (e.g. fewer than 5 reads in less than 20% of samples) - could be artifacts
physeq_object_filtered <- filter_taxa(physeq_object, function(x) sum(x > 5) > (0.2*length(x)), prune = TRUE)

## reads ----
# check total number of reads per sample
# unfiltered & filtered look pretty different!

### unfiltered ----
sample_sums(physeq_object)

# viz library size distribution
hist(sample_sums(physeq_object), main="Distribution of sample read counts", xlab="Read Counts")
min_lib_size <- min(sample_sums(physeq_object)) # 108130
max_lib_size <- max(sample_sums(physeq_object)) # 185656
185656-108130 # = 77526

# total sequence reads across samples
sum(OTU) # = 5604530 

# viz differences in taxa abundance among samples 
plot_bar(physeq_object)

### filtered ----
sample_read_counts <- sample_sums(physeq_object_filtered)

# viz library size distribution
hist(sample_read_counts, main="Distribution of (filtered) sample read counts", xlab="Read Counts")
min_lib_size <- min(sample_read_counts) # 16358
108130 - 16358 # 91772 less than before filter
max_lib_size <- max(sample_read_counts) # 185656
185656-108130 # = 77526

# total sequence reads across samples
sum(sample_read_counts) # 3123265
3123265/5604530 # 56 % of original read counts

# viz differences in taxa abundance among samples 
plot_bar(physeq_object_filtered)



## alpha ----
plot_richness(physeq_object, x="treatment", measures=c("Shannon", "Observed"))

## beta ----
ord_nmds <- ordinate(physeq_object, method="NMDS")
plot_ord_samples <- plot_ordination(physeq_object, ord_nmds, color = "treatment")

# left off here ----
## normalize data ----

# rarefy full otu with rarefaction
set.seed(1782) # set seed for analysis reproducibility
OTU_rar = rarefy_even_depth(otu_table(data_phylo), rngseed = TRUE, replace = FALSE)
# 110 OTUs were removed because they are no longer present in any sample after random subsampling
# second time I did this (including PT), only 86 were removed - not sure why

data_otu_rar = data.frame(otu_table(OTU_rar)) # create otu-class phyloseq object
sum(data_otu_rar) # 1180188, total num squ(reads)
sum_seq_rar <-rowSums(data_otu_rar) 
sum_seq_rar # now each sample has 32783 reads/counts (randomly selected)
plot(sum_seq_rar, main=c("Counts / Sample - Rarefied ASVs"), xlab=c("Samples"))     
# create a new phyloseq-class object 
data_phylo_rar <- phyloseq(OTU_rar, TAX, SAM, PT) 
data_phylo_rar # has 2483 taxa and 36 samples (filtered OTU had 351 taxa after rarefaction)
# second time ran this (with PT), has 2507 taxa - not sure why diff
# (2483/6843 = 36.28526% of original taxa remain, 63.71% removed)
plot_bar(data_phylo_rar) # see frequency distribution (abundance vs samples)

# creating a second rarefied phyloseq object, without PT bc having issues with PERMANOVA results
set.seed(1782) # set seed for analysis reproducibility
OTU_rar2 = rarefy_even_depth(otu_table(data_phylo2), rngseed = TRUE, replace = FALSE)
# 110 OTUs removed
data_phylo_rar2 <- phyloseq(OTU_rar2, TAX, SAM)
data_phylo_rar2

###### a-div: indix calcs & graphs
dim(data_otu_rar)
# all indices 
data_richness_rar <- estimateR(data_otu_rar) # calculate richness and Chao1 using vegan package

data_evenness_rar <- diversity(data_otu_rar) / log(specnumber(data_otu_rar)) # calculate evenness index using vegan package
data_evenness_rar

data_shannon_rar <- diversity(data_otu_rar, index = "shannon") # calculate Shannon index using vegan package

data_alphadiv_rar <- cbind(data_grp2, t(data_richness_rar), data_shannon_rar, data_evenness_rar) # combine all indices in one data table
data_alphadiv_rar
# use tidy format (longer table)(alphadiv_index as 1 column)(obs_values as 1 column)
data_alphadiv_rar_tidy <- 
  data_alphadiv_rar %>%
  mutate(sample_id = rownames(data_alphadiv_rar)) %>%
  gather(key   = alphadiv_index,
         value = obs_values,
         -plot, -site, -treatment)
data_alphadiv_rar_tidy

# plots (original, not accounting for nested design, see bottom of lmer section for adjusted graphs)
P1_rar <-ggplot(data_alphadiv_rar, aes(x=treatment, y=S.obs)) +
  geom_boxplot(fill=c("#247202","#a60000","#ea9a05","#012476")) +
  labs(title= 'Richness', x= ' ', y= '', tag = "A") +
  geom_point()
P1_rar
# adjust plot
P1_rar + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"))

P2_rar <- ggplot(data_alphadiv_rar, aes(x=treatment, y=S.chao1)) +
  geom_boxplot(fill=c("#247202","#a60000","#ea9a05","#012476")) +
  labs(title= 'Chao1', x= ' ', y= '', tag = "B") +
  geom_point()
P2_rar

P3_rar <- ggplot(data_alphadiv_rar, aes(x=treatment, y=data_evenness)) +
  geom_boxplot(fill=c("#247202","#a60000","#ea9a05","#012476")) +
  labs(title= 'Eveness', x= ' ', y= '', tag = "C") +
  geom_point()
P3_rar

P4_rar <- ggplot(data_alphadiv_rar, aes(x=treatment, y=data_shannon)) +
  geom_boxplot(fill=c("#247202","#a60000","#ea9a05","#012476")) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "D") +
  geom_point()
P4_rar

###### a-div: lmer & eemeans - richness
Adiv_rar_lmer_richness
Adiv_rar_lmer_richness <- lmer(S.obs~treatment+(1|site), data = data_alphadiv_rar)
# site as random factor (plot not included bc smallest unit = DNA pooled by plot, each "sample" = "plot")
summary(Adiv_rar_lmer_richness)
anova(Adiv_rar_lmer_richness)
ranova(Adiv_rar_lmer_richness)

# test assumptions
qqnorm(residuals(Adiv_rar_lmer_richness))
qqline(residuals(Adiv_rar_lmer_richness))
plot(Adiv_rar_lmer_richness) # looks pretty good

boxplot(residuals(Adiv_rar_lmer_richness)~data_alphadiv_rar$treatment)
boxplot(residuals(Adiv_rar_lmer_richness)~data_alphadiv_rar$site) # bit of hetero(variance) among sites, less among treat&plot

shapiro.test(data_alphadiv_rar$S.obs) #  p=0.3641 (want it >0.05).. this is good, right?

### Just to test, did not use transformed data in rest of analyses ###
# try transformations to make more normal   ???   Ex code from filtered, rarefied data & shannon
# adding log transformed shannon as new column to data_alphadiv
data_alphadiv_rar$logT_rar_richness <- log(data_alphadiv_rar$data_)
str(data_alphadiv) # see new columna dded
shapiro.test(data_alphadiv$logT_shannon) # even less normal, p=0002664
# adding sqrt trans shannon as new column
data_alphadiv$sqrt_shannon <- sqrt(data_alphadiv$data_shannon)
str(data_alphadiv) # see new column dded
shapiro.test(data_alphadiv$sqrt_shannon) # o=0.001124
# maybe best to stick with original shannon for now

# DHARMa for richness
# just looking at dispersion
testDispersion(Adiv_rar_lmer_richness)
# recommended step 1 (calc residuals)
DharmaOutput_rich <- simulateResiduals(fittedModel = Adiv_rar_lmer_richness, plot = F)
# access the residuals
residuals(DharmaOutput_rich)
# see the QQ plot
plot(DharmaOutput_rich)
# see the predictor plot
plotResiduals(DharmaOutput_rich, form = data_alphadiv_rar$treatment)
# homogeneity should range from 0.25-0.75, all my treat go beyond

#run post-hoc test on Adiv_lmer_shannon with emmeans
Adiv_rich_rar_tukey <- emmeans(Adiv_rar_lmer_richness, pairwise~treatment, adjusted="tukey")
Adiv_rich_rar_tukey
# no diff detected btwn groups

###### a-div: lmer & eemeans - evenness
data_evenness_rar 
# is column name for evenness values in data_alphadiv_rar
Adiv_rar_lmer_even <- lmer(data_evenness_rar~treatment+(1|site), data = data_alphadiv_rar)     
# site as random factor (plot not included bc smallest unit = DNA pooled by plot, each "sample" = "plot")

summary(Adiv_rar_lmer_even)
anova(Adiv_rar_lmer_even)
ranova(Adiv_rar_lmer_even)

# test assumptions
qqnorm(residuals(Adiv_rar_lmer_even))
qqline(residuals(Adiv_rar_lmer_even))
plot(Adiv_rar_lmer_even) # looks odd

boxplot(residuals(Adiv_rar_lmer_even)~data_alphadiv_rar$treatment) # bit of hetero(variance)
boxplot(residuals(Adiv_rar_lmer_even)~data_alphadiv_rar$site) # bit of hetero(variance) among sites

shapiro.test(data_alphadiv_rar$data_evenness_rar) #  p=0.0003521 (want it >0.05).. this is not good, right?

# try transformations to make more normal   ??? 
# adding log transformed evenness as new column to data_alphadiv_rar
data_alphadiv_rar$logT_rar_even <- log(data_alphadiv_rar$data_evenness_rar)
str(data_alphadiv_rar) # see new columna dded
shapiro.test(data_alphadiv_rar$logT_rar_even) # even less normal, p=8.081e-06
# adding sqrt trans even as new column
data_alphadiv_rar$sqrt_rar_sh <- sqrt(data_alphadiv_rar$data_shannon_rar)
str(data_alphadiv_rar) # see new column dded
shapiro.test(data_alphadiv_rar$sqrt_rar_sh) # p-value = 0.0001035
# maybe best to stick with original shannon for now

# DHARMa for evenness
# just looking at dispersion
testDispersion(Adiv_rar_lmer_even)
# recommended step 1 (calc residuals)
DharmaOutput <- simulateResiduals(fittedModel = Adiv_rar_lmer_even, plot = F)
# access the residuals
residuals(DharmaOutput)
# see the QQ plot
plot(DharmaOutput)
# see the predictor plot
plotResiduals(DharmaOutput, form = data_alphadiv_rar$treatment)
# homogeneity should range from 0.25-0.75, my control treat goes beyond

# run post-hoc test on Adiv_lmer_shannon with emmeans
emmeans(Adiv_rar_lmer_even, pairwise~treatment, adjusted="tukey")
# no diff detected btwn groups

# generate tables with basic summary stats for richness & evenness
data_alphadiv_rar
# richness
r.mean <- tapply(data_alphadiv_rar$S.obs, data_alphadiv_rar$treatment, mean)
r.sd <- tapply(data_alphadiv_rar$S.obs, data_alphadiv_rar$treatment, sd)
r.n <- tapply(data_alphadiv_rar$S.obs, data_alphadiv_rar$treatment, length)
r.se <- r.sd/sqrt(r.n)
r2 <- data.frame(mean=r.mean,sd=r.sd,n=r.n,se=r.se)
r2
# evenness
e.mean <- tapply(data_alphadiv_rar$data_evenness_rar, data_alphadiv_rar$treatment, mean)
e.sd <- tapply(data_alphadiv_rar$data_evenness_rar, data_alphadiv_rar$treatment, sd)
e.n <- tapply(data_alphadiv_rar$data_evenness_rar, data_alphadiv_rar$treatment, length)
e.se <- e.sd/sqrt(e.n)
e2 <- data.frame(mean=e.mean,sd=e.sd,n=e.n,se=e.se)
e2

# sum stats with summarySE
# richness
# calc means for each treat group, per site
rare.agg <-summarySE(data_alphadiv_rar, measurevar = "S.obs", groupvars = c("treatment","site"))
rare.agg
# calc overall means per treat, based on site means
rare.agg_gr <- summarySE(rare.agg, measurevar = "S.obs", groupvars = c("treatment"))
rare.agg_gr
# alt graph
rare_bar <- ggplot(rare.agg_gr, aes(x=treatment, y=S.obs)) +
  geom_bar(stat = "identity", col="black", lwd=0.35, fill=c("#319903","#f50505","#ea9a05","#1859f0"), width = 0.8) + ylim(c(0,250)) + 
  geom_errorbar(aes(ymin=S.obs-se, ymax=S.obs+se), width=0.2, size=0.4)
rare_bar

rare_bar + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))
# evenness
# calc means for each treat group, per site
even.agg <-summarySE(data_alphadiv_rar, measurevar = "data_evenness_rar", groupvars = c("treatment","site"))
even.agg
# calc overall means per treat, based on site means
even.agg_gr <- summarySE(even.agg, measurevar = "data_evenness_rar", groupvars = c("treatment"))
even.agg_gr
# alt graph
even_bar <- ggplot(even.agg_gr, aes(x=treatment, y=data_evenness_rar)) +
  geom_bar(stat = "identity", col="black", lwd=0.35, fill=c("#319903","#f50505","#ea9a05","#1859f0"), width = 0.8) + ylim(c(0,1))  +
  geom_errorbar(aes(ymin=data_evenness_rar-se, ymax=data_evenness_rar+se), width=0.2, size=0.4)
even_bar

even_bar + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))

###### b-div: Bray-Curtis NMDS
# create distance matrix. Bray-Curtis using the vegan package. Rarefied data
dim(data_otu_rar)
set.seed(100)
BC_rar <- as.matrix(vegdist(data_otu_rar, method = "bray")) 
BC_rar2 <- as.matrix(vegdist(OTU_rar2, method = "bray")) 
# a peek at the first five rows / columns
BC_rar[1:5, 1:5]
BC_rar2[1:5, 1:5]
# permanova (using vegan)
set.seed(100)
perma_BC_rar <- adonis2(BC_rar ~ treatment*site, data = data_grp2, permutations = 10000 )
perma_BC_rar
# pairwiseAdonis (using phyloseq object & pairwiseAdonis package)
set.seed(100)
pairwise_BC_rar <- pairwise.adonis(BC_rar, phyloseq::sample_data(data_phylo_rar)$treatment)
pairwise_BC_rar # also tried with *site included, but got error

# run NMDS ordination
NMDS_BC_rar <- ordinate(data_phylo_rar, "NMDS", "bray") # got warning
# get correct colors for treats
cols <- c("control" = "#319903","burn" = "#f50505","mechanical" = "#ea9a05","mech-burn" = "#1859f0")
# plot NMDS
plot_NMDS_BC_rar <- plot_ordination(data_phylo_rar, NMDS_BC_rar, color = "treatment") +
  geom_point(size = 3) + scale_color_manual(values=cols)
# adjust plot
plot_NMDS_BC_rar + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"))

###### b-div: UniFrac NMDS
# create distance matrix (phyloseq)
set.seed(100)
wUF_rar <- phyloseq::distance(data_phylo_rar, method = "wunifrac")
# Warning message:
# In matrix(tree$edge[order(tree$edge[, 1]), ][, 2], byrow = TRUE,  :
# data length [4959] is not a sub-multiple or multiple of the number of rows [2480]
# continuing anyway
# view
as.matrix(wUF_rar)[1:6, 1:6]
# permanova (vegan)
set.seed(100)
permanova_wUF_rar <- adonis2(wUF_rar ~ treatment*site, data = data_grp2, permutations = 10000 )
permanova_wUF_rar
# pairwiseAdonis (using phyloseq object & pairwiseAdonis package)
pairwise_wUF_rar <- pairwise.adonis(wUF_rar, phyloseq::sample_data(data_phylo_rar)$treatment)
pairwise_wUF_rar
# run ordination
NMDS_wUF_rar <- ordinate(data_phylo_rar, method = "NMDS", distance = wUF)
NMDS_wUF_rar
# plot 
plot_NMDS_wUF_rar <- plot_ordination(data_phylo_rar, NMDS_wUF_rar, color = "treatment") +
  geom_point(size = 3) + scale_color_manual(values=cols)
plot_NMDS_wUF_rar
# adjust plot
plot_NMDS_wUF_rar + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white"))

###### Abundance Pie charts & Bar plots
# used this Loue River tutorial for example code/ methodology 
# https://scienceparkstudygroup.github.io/microbiome-lesson/aio/index.html

# Using unfiltered ASVs normalized with rarefaction & data_taxo2 (with functional guild info)
data_otu_rar[1:5, 1:6] 
data_phylo_rar # 2483 taxa by 15 tax ranks (including functional columns)
sum(data_otu_rar) # 1180188, total num squ reads
dim(data_otu_rar) # [36, 2483]

# prep data
# add the rownames in a column in order to be able to use it later with dplyr
data_taxo2$OTU_id <- rownames(data_taxo2) 
# new temp grp table with grouping factors & OTU occurrences
# using original data_grp2 file
data_grp2_temp <- data_grp2
# to be able to use the column later with dplyr 
data_grp2_temp$sample_id <- rownames(data_grp2_temp) 
data_grp2_temp
# new temp otu table 
data_otu_rar_temp <- data_otu_rar
# to be able to use the column later with dplyr
data_otu_rar_temp$sample_id <- rownames(data_otu_rar_temp) 
dim(data_otu_rar_temp) # [36, 2484]
# join tables
data_otu_grp <- inner_join(data_grp2_temp, data_otu_rar_temp, by = "sample_id")
# see the new table (meta & asv frequ combined, and had "sample_id" which will be useful)
data_otu_grp[1:5, 1:7]
dim(data_otu_grp) # [36, 2487]
# aggregate the new data table, data_otu_grp, by treatment and
# add the aggregated data table to the data table that include taxonomic information, data_taxo2
# Calculate sum of counts per treatment
data_otu_rar_treat <- aggregate(data_otu_grp[, 5:2487], by = list(treatment=data_otu_grp$treatment), sum)
data_otu_rar_treat[, 1:7]
dim(data_otu_rar_treat) # 4 groups (treats) & 2484 OTUs 
# Transpose the data table, calculate & add total_count for all samples, and add a column with OTU_id
data_otu_rar_treat_temp <- as.data.frame(t(as.matrix(data_otu_rar_treat[,2:dim(data_otu_rar_treat)[2]])))
dim(data_otu_rar_treat_temp) # [2483, 4]
colnames(data_otu_rar_treat_temp) <- data_otu_rar_treat[,1]
data_otu_rar_treat_temp$total_counts <- rowSums(data_otu_rar_treat_temp)
data_otu_rar_treat_temp$OTU_id <- rownames(data_otu_rar_treat_temp)
dim(data_otu_rar_treat_temp) # [2483, 6] (2 new columns)
# confirm, has columns: 4 treats & total_counts & OTU_id
head(data_otu_rar_treat_temp)
# Transpose the full unmodified OTU data table (data_otu_rar) & add col for OTU_id
data_otu_rar_t <- as.data.frame(t(as.matrix(data_otu_rar)))
data_otu_rar_t$OTU_id <- rownames(data_otu_rar_t)
# Combine the two OTU tables, the aggregated OTU table (data_otu_rar_treat_temp) with data_otu_rar_t
data_otu_rar_temp <- inner_join(data_otu_rar_treat_temp, data_otu_rar_t, by = "OTU_id")
# viz
data_otu_rar_temp[1:5, 1:5] # looks really nice!
str(data_otu_rar_temp) # 2483 obs 42 vars (4 treats, total_count, OTU_id, all 36 sites)
dim(data_otu_rar_temp) # [2483, 42]
# Combine the taxonomic table, data_taxo_filt_rar, with the new OTU table, data_otu_rar_temp
data_otu_taxo2_rar <- inner_join(data_taxo2, data_otu_rar_temp, by = "OTU_id")
# viz
dim(data_otu_taxo2_rar) # 2483 OTUs & 57 cat vars (taxa groups, treats, other info, samples)
data_otu_taxo2_rar[1:5,] # has all the original tax & grp2 data with OTU data (all in 1)
str(data_otu_taxo2_rar)

### pie charts
# We will sum all the OTU that belonging to the same Phylum
# Filter abundant Phyla (those with >=1% of total proportion)
abundant_com_rar <- data_otu_taxo2_rar %>%
  select(phylum, total_counts) %>%
  group_by(phylum) %>%
  summarize(total_counts_per_phylum = sum(total_counts)) %>%
  ungroup() %>%
  mutate(total_counts_percentage = total_counts_per_phylum / sum(total_counts_per_phylum) * 100) %>%
  filter(total_counts_percentage >= 1) %>%
  select(phylum, total_counts_percentage)
abundant_com_rar # has 4 phyla
# Compile rare Phyla
# Data table global communities simplified
rare_com_rar <- data_frame(phylum=c("Rare_Phyla"),total_counts_percentage=100-sum(abundant_com_rar$total_counts_percentage))
rare_com_rar
global_com_rar <- bind_rows(abundant_com_rar,rare_com_rar)
global_com_rar
# Piechart
pie_rar_phylum<- ggplot(global_com_rar, aes(x="", y=total_counts_percentage, fill=phylum)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank(), legend.position="bottom", legend.text = element_text(size = 8))
pie_rar_phylum

# We will sum all the OTU that belonging to the same Genus
# Filter abundant Genera (those with >=1% of total proportion)
abundant_com_rarG <- data_otu_taxo2_rar %>%
  select(genus, total_counts) %>%
  group_by(genus) %>%
  summarize(total_counts_per_genus = sum(total_counts)) %>%
  ungroup() %>%
  mutate(total_counts_percentage = total_counts_per_genus / sum(total_counts_per_genus) * 100) %>%
  filter(total_counts_percentage >= 1) %>%
  select(genus, total_counts_percentage)
abundant_com_rarG # has 20 Genera
# Compile rare Genera
# Data table global communities simplified
rare_com_rarG <- data_frame(genus=c("Rare_Genera"),total_counts_percentage=100-sum(abundant_com_rarG$total_counts_percentage))
rare_com_rarG
global_com_rarG <- bind_rows(abundant_com_rarG,rare_com_rarG)
# Piechart
pie_rar_genus<- ggplot(global_com_rarG, aes(x="", y=total_counts_percentage, fill=genus)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank(), legend.position="bottom", legend.text = element_text(size = 8))
pie_rar_genus

# same for trophic_mode
# Filter abundant taxa in trophic_mode (those with >=1% of total proportion)
abundant_com_rarTM <- data_otu_taxo2_rar %>%
  select(trophic_mode, total_counts) %>%
  group_by(trophic_mode) %>%
  summarize(total_counts_per_trophic_mode = sum(total_counts)) %>%
  ungroup() %>%
  mutate(total_counts_percentage = total_counts_per_trophic_mode / sum(total_counts_per_trophic_mode) * 100) %>%
  filter(total_counts_percentage >= 1) %>%
  select(trophic_mode, total_counts_percentage)
abundant_com_rarTM # has 7 TMs
# Compile rare TMs
rare_com_rarTM <- data_frame(trophic_mode=c("Rare_Taxa"),total_counts_percentage=100-sum(abundant_com_rarTM$total_counts_percentage))
global_com_rarTM <- bind_rows(abundant_com_rarTM,rare_com_rarTM)
global_com_rarTM
# Piechart
pie_rar_TM<- ggplot(abundant_com_rarTM, aes(x="", y=total_counts_percentage, fill=trophic_mode)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank(), legend.position="right", legend.text = element_text(size = 8))
pie_rar_TM

# same for guild
# Filter abundant taxa in guild (GD) (those with >=1% of total proportion)
str(data_otu_taxo2_rar$trophic_mode) # has 97 guilds & 8 trophic modes
abundant_com_rarGD <- data_otu_taxo2_rar %>%
  select(guild, total_counts) %>%
  group_by(guild) %>%
  summarize(total_counts_per_guild = sum(total_counts)) %>%
  ungroup() %>%
  mutate(total_counts_percentage = total_counts_per_guild / sum(total_counts_per_guild) * 100) %>%
  filter(total_counts_percentage >= 1) %>%
  select(guild, total_counts_percentage)
abundant_com_rarGD # has 10 Guilds, only those with >= 1% relative abundance... 
# Compile rare TMs
rare_com_rarGD <- data_frame(guild=c("Rare_Taxa"),total_counts_percentage=100-sum(abundant_com_rarGD$total_counts_percentage))
global_com_rarGD <- bind_rows(abundant_com_rarGD,rare_com_rarGD)
global_com_rarGD
# Piechart
pie_rar_GD<- ggplot(abundant_com_rarGD, aes(x="", y=total_counts_percentage, fill=guild)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank(), legend.position="right", legend.text = element_text(size = 8)) +
  guides(color=guide_legend(nrow=6, ncol=2, byrow=TRUE))
pie_rar_GD

### bar charts
# phylum
# Sum per phylum by treat
view(data_otu_taxo2_rar) # (treatments are columns 17:20)
rar_treat_phylum_temp <- aggregate(data_otu_taxo2_rar[, 17:20], by=list(phylum=data_otu_taxo2_rar$phylum), sum)
view(rar_treat_phylum_temp) #14 phyla (1 unassigned)
# Tidying the data set
rar_treat_phylum <- rar_treat_phylum_temp %>%
  gather(key   = treatment, value = obs_values,-phylum)
view(rar_treat_phylum) # has x features
# filter 
rar_treat_phylum_filt <- filter(rar_treat_phylum, obs_values >=1000)
view(rar_treat_phylum_filt) # has only 5 phyla
# Barplot
bar_phylum <- rar_treat_phylum_filt %>%
  mutate(treatment = fct_relevel(treatment, "control", "mechanical", "burn", "mech-burn")) %>%
  ggplot(., aes(x = treatment, y = obs_values, fill = phylum)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8, color = "black", lwd=0.35) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title= 'Fungal community composition at Phylum level', x= '', y= 'Proportion of the total counts')

bar_phylum + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))

# genus
# Sum per genus 
str(data_otu_taxo2_rar) # (treatments are columns 17:20)
rar_treat_genus_temp <- aggregate(data_otu_taxo2_rar[, 17:20], by=list(genus=data_otu_taxo2_rar$genus), sum)
view(rar_treat_genus_temp) # has 435 genera in rar data
# Tidying the data set
rar_treat_genus <- rar_treat_genus_temp %>%
  gather(key   = treatment, value = obs_values,-genus)
dim(rar_treat_genus) # has 1740 features (all assigned genera?)
# filter 
rar_treat_genus_filt <- filter(rar_treat_genus, obs_values >=5000)
view(rar_treat_genus_filt) # has 157 features (genera) >1000 
# 110 >2000
# 79 >3000
# 48 >5000 (11.03 % of 435 rarefied genera)
# Barplot
bar_genus <- rar_treat_genus_filt %>%
  mutate(treatment = fct_relevel(treatment, "control", "mechanical", "burn", "mech-burn")) %>%
  ggplot(., aes(x = treatment, y = obs_values, fill = genus)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8, color = "black", lwd=0.35) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title= 'Fungal community composition at Genus level', x= '', y= 'Proportion of the total counts')
bar_genus + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))

# trophic_mode
# Sum per trophic_mode 
rar_treat_TM_temp <- aggregate(data_otu_taxo2_rar[, 17:20], by=list(trophic_mode=data_otu_taxo2_rar$trophic_mode), sum)
view(rar_treat_TM_temp)
# Tidying the data set
rar_treat_TM <- rar_treat_TM_temp %>%
  gather(key   = treatment, value = obs_values,-trophic_mode)
view(rar_treat_TM) # has 32 observations
# filter [skip]
# Barplot
bar_TM <- rar_treat_TM %>%
  mutate(treatment = fct_relevel(treatment, "control", "mechanical", "burn", "mech-burn")) %>%
  ggplot(., aes(x = treatment, y = obs_values, fill = trophic_mode)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8, color = "black", lwd=0.35) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title= 'Fungal community composition by Trophic Mode', x= '', y= 'Proportion of the total counts')
bar_TM + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))

# Guild
# Sum per guild 
rar_treat_guild_temp <- aggregate(data_otu_taxo2_rar[, 17:20], by=list(guild=data_otu_taxo2_rar$guild), sum)
str(rar_treat_guild_temp) # has 85 guilds
view(rar_treat_guild_temp)
# Tidying the data set
rar_treat_guild <- rar_treat_guild_temp %>%
  gather(key   = treatment, value = obs_values,-guild)
dim(rar_treat_guild) # has 340 features (all assigned genera?)
# filter 
rar_treat_guild_filt <- filter(rar_treat_guild, obs_values >=4000)
view(rar_treat_guild_filt) # has 65 features (guilds) >1000 
# 48 >2000
# 34 >4000  (40% of 85 guilds in rarefied data)
(34/85)*100
# Barplot
bar_guild <- rar_treat_guild_filt %>%
  mutate(treatment = fct_relevel(treatment, "control", "mechanical", "burn", "mech-burn")) %>%
  ggplot(., aes(x = treatment, y = obs_values, fill = guild)) +
  geom_bar(position = "fill", stat = "identity", width = 0.8, color = "black", lwd=0.35) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title= 'Fungal community composition at Guild level', x= '', y= 'Proportion of the total counts')
bar_guild + theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  linewidth = 0.5, linetype = "solid"),
  panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                  colour = "#d9d4d4"), 
  panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                  colour = "white"))



