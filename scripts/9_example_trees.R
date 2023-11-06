library(treeio)
library(ggtree)
library(ggplot2)
library(patchwork)
library(ape)
library(tibble)
library(tidyr)


tr2 <- read.beast.newick("simulated_trees/results/true_all/s1_h125.0_m4_6.tre")
tr3 <- read.beast.newick("simulated_trees/results/true_div/s1_h125.0_m4_6.tre")
tr4 <- read.beast.newick("simulated_trees/results/unknown_rates/s1_h125.0_m4_6.tre")

  
make_shiftplot <- function(tr, limits){
  if ("N" %in% colnames(tr@data)){
    tr@data$nshift <- sapply(tr@data$N, sum)
  }
  ## Number of shifts (N)
  p <- ggtree(tr, aes(color = nshift)) +
    scale_colour_gradient("", low = "gray", high = "red", limits = limits) +
    theme(legend.pos = c(0.125, 0.875))
    #ggtitle()
  return(p)
}

make_divplot <- function(tr, limits){
  if ("mean_netdiv" %in% colnames(tr@data)){
    tr@data$r <- tr@data$mean_netdiv  
  }
  ## net-diversification rate
  p <- ggtree(tr, aes(color = r)) +
    
    scale_colour_gradient("", low = "gray", high = "#619CFF", limits = limits) +
    theme(legend.pos = c(0.125, 0.875))
  #ggtitle()
  return(p)
}

################################
##
## tree that did not work well
##  trees: 100 Ma, model 1, replicate 13
################################
tr1 <- read.beast.newick("simulated_trees/s1_h100.0_m1_13.tre")
tr2 <- read.beast.newick("simulated_trees/results/unknown_rates/s1_h100.0_m1_13.tre")
min_max <- c(
  min(0, min(tr2@data$nshift)),
  max(1, max(tr2@data$nshift))
)

r_min_max <- c(
  min(min(tr1@data$r), min(tr2@data$mean_netdiv)),
  max(max(tr1@data$r), max(tr2@data$mean_netdiv))
)

## good trees 
## 6, 22 (100 Ma)
p1a <- make_shiftplot(tr1, min_max) + ggtitle("true shifts")
p1b <- make_shiftplot(tr2, min_max) + ggtitle("inferred shifts")
p1c <- make_divplot(tr1, r_min_max) + ggtitle("true net-diversification")
p1d <- make_divplot(tr2, r_min_max) + ggtitle("inferred net-diversification")


upper_panel <- p1a + p1b + p1c + p1d +
  plot_layout(ncol = 4, guides = "keep", tag_level = "new")

#ggsave("figures/fourpaneltree1.pdf", height = 125, width = 250, units = "mm")

################################
##
## tree that worked quite well
##
################################
tr1 <- read.beast.newick("simulated_trees/s1_h100.0_m4_22.tre")
tr2 <- read.beast.newick("simulated_trees/results/unknown_rates/s1_h100.0_m4_22.tre")
min_max <- c(
  min(0, min(tr2@data$nshift)),
  max(1, max(tr2@data$nshift))
)

## good trees 
## 6, 22 (100 Ma)
r_min_max <- c(
  min(min(tr1@data$r), min(tr2@data$mean_netdiv)),
  max(max(tr1@data$r), max(tr2@data$mean_netdiv))
)

p1a <- make_shiftplot(tr1, min_max) + ggtitle("true shifts")
p1b <- make_shiftplot(tr2, min_max) + ggtitle("inferred shifts")
p1c <- make_divplot(tr1, r_min_max) + ggtitle("true net-diversification")
p1d <- make_divplot(tr2, r_min_max) + ggtitle("inferred net-diversification")

lower_panel <- p1a + p1b + p1c + p1d +
  plot_layout(ncol = 4, guides = "keep", tag_level = "new")

upper_panel / lower_panel +
  plot_annotation(tag_levels = c("a","1"), tag_suffix = ")")

ggtree(tr2, aes(color = log(shift_bf)))
ggtree(tr2, aes(color = nshift))
tr2@data[tr2@data$shift_bf > 10,]
tr1@data[sapply(tr1@data$N, sum) > 0.5,]

ggsave("figures/tree_examples.pdf", height = 225, width = 270, units = "mm")
