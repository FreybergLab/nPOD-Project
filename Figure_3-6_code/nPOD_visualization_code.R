library(data.table)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(egg)
library(grid)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., nudge = 0, draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             # now nudge them apart
                             newdata$x <- ifelse(newdata$group %% 2 == 1,
                                                 newdata$x - nudge,
                                                 newdata$x + nudge)
                             
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              nudge = 0, draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, nudge = nudge, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#### START here ####
dt <- read.csv("../Data/quad_data.csv")
dt$unique.ID <- paste(dt$Object.Id, dt$ImageID, sep = "_")
dt <- dt[!duplicated(dt$unique.ID),]
dt$Gender <- dt$Sex
dim(dt)

# donor level
# transform raw data to 100 percentile value 
cell.type <- "soma" # could be ins, glu, soma, PP
cell.type.percent <- paste0("Percent.", cell.type)
cell.type.col.idx <- which(colnames(dt) == cell.type.percent)
donor.dt <- NULL
for(dd in unique(dt$Donor)){
  tmp.donor <- dt[dt$Donor == dd,]
  donor.summary <- data.frame(Donor = tmp.donor$Donor[1],
                              Age = tmp.donor$Age[1],
                              AgeCat = tmp.donor$AgeCat[1],
                              Gender = tmp.donor$Gender[1],
                              Diagnosis = tmp.donor$Diagnosis[1],
                              Disease.Duration = tmp.donor$Disease.Duration[1],
                              Region = rep(c("Head", "Body", "Tail"), each = 100))
  tmp.Head <- tmp.donor[tmp.donor$Region == "Head",]
  tmp.body <- tmp.donor[tmp.donor$Region == "Body",]
  tmp.Tail <- tmp.donor[tmp.donor$Region == "Tail",]
  donor.summary$quantile <- rep(c(1:100), 3)
  donor.summary$cell.proportion <- c(quantile(tmp.Head[,cell.type.col.idx], probs = seq(0.01,1,0.01)),
                                     quantile(tmp.body[,cell.type.col.idx], probs = seq(0.01,1,0.01)),
                                     quantile(tmp.Tail[,cell.type.col.idx], probs = seq(0.01,1,0.01)))
  donor.summary$cum.proportion <- sum(donor.summary$cell.proportion, na.rm = T)
  donor.dt <- rbind(donor.dt, donor.summary)
}
donor.dt <- donor.dt[order(donor.dt$Diagnosis, donor.dt$Disease.Duration, decreasing = TRUE),]
donor.dt$Donor <- factor(donor.dt$Donor, 
                         levels = unique(donor.dt$Donor))
donor.dt$Region <- factor(donor.dt$Region, levels = c("Head", "Body", "Tail"))
donor.dt$Diagnosis <- factor(donor.dt$Diagnosis, levels = c("T1D", "ND"))
donor.dt$rev_quantile <- 101 - donor.dt$quantile

## grouped heatmap
grouped.heatmap <- ggplot(donor.dt, aes(x = Donor, y = rev_quantile)) +
  geom_tile(aes(fill = cell.proportion), na.rm = T)  +
  facet_grid(Region ~ AgeCat*Gender,
             scales = "free", switch = "y") +
  theme_classic() +
  scale_fill_gradientn(colours = c("lightblue", "lightyellow", "red"),
                       values = c(0,0.5,1),
                       na.value = 'lightgrey') +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(y = "", fill = cell.type.percent) +
  theme(panel.spacing = unit(2, "lines"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.margin = unit(c(10,5.5,0,0), 'pt'), 
        strip.text = element_text(size = 12),
        legend.position = "none")

# change color of facet strip
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
grouped.heatmap.grob <- ggplot_gtable(ggplot_build(grouped.heatmap))
strips <- which(grepl('strip-', grouped.heatmap.grob$layout$name))[-5:-7] # first 3 are region groups
pal <- gg_color_hue(4)
for (i in seq_along(strips)) {
  k1 <- which(grepl('rect', grouped.heatmap.grob$grobs[[strips[i]]]$grobs[[1]]$childrenOrder)) # age group
  k2 <- which(grepl('rect', grouped.heatmap.grob$grobs[[strips[i]]]$grobs[[2]]$childrenOrder)) # gender group
  grouped.heatmap.grob$grobs[[strips[i]]]$grobs[[1]]$children[[k1]]$gp$fill <- pal[i]
  grouped.heatmap.grob$grobs[[strips[i]]]$grobs[[2]]$children[[k2]]$gp$fill <- pal[i]
}


side.bar.duration <- ggplot(donor.dt, aes(x = Donor, y = "1")) + 
  geom_tile(aes(fill = Disease.Duration)) +
  scale_fill_gradientn(colours = c("lightgreen", "purple4"),
                       na.value = 'white') +
  facet_grid(~ AgeCat * Gender ,
             scales = "free") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border=element_rect(colour="grey2",fill=NA, linewidth=1),
        plot.margin = unit(c(0,5.5,0,20), 'pt'))

side.bar.T1D <- ggplot(donor.dt, aes(x = Donor, y = "1")) + 
  geom_tile(aes(fill = Diagnosis)) +
  scale_fill_discrete(breaks = c("T1D","No.diagnosis"),
                      type = c("black", "white")) +
  facet_grid(~ AgeCat * Gender,
             scales = "free") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border=element_rect(colour="grey2",fill=NA, linewidth=1),
        plot.margin = unit(c(0,5.5,0,20), 'pt'))

donor.dt$gp <- factor(paste0(donor.dt$Gender, " (", donor.dt$AgeCat, ")"),
                      levels = c("Female (<20_years)", "Male (<20_years)", "Female (20-40_years)", "Male (20-40_years)"))

bottom.violin <- ggplot(donor.dt, aes(x = gp, y = cell.proportion, fill = Diagnosis, color = gp)) + 
  geom_split_violin(scale = "width", lwd = 0.8, nudge = 0.02) +
  facet_grid(Region ~ .,
             scales = "free") +
  scale_fill_discrete(breaks = c("T1D", "No.diagnosis"),
                      type = c("grey2", "white")) +
  theme_classic() +
  labs(y = "Cell Proportion") +
  scale_y_continuous(expand=c(0,0), position = "right") +
  theme(panel.spacing = unit(2, "lines"),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(size = 15),
        plot.margin = unit(c(50,5.5,0,0), 'pt'))


# output
graphics.off()   # closes every open graphics device, including the stale one
fig <- ggpubr::ggarrange(
  as.ggplot(grouped.heatmap.grob), bottom.violin,
  side.bar.T1D, NULL,
  side.bar.duration, NULL,
  nrow = 3, ncol = 2, heights = c(10, 0.3, 0.3), widths = c(10, 3),
  legend = "none"
)
ggsave(paste0(cell.type, "_grouped_v4.pdf"), fig, width = 15, height = 9)