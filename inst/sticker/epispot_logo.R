require(hexSticker)
require(dplyr)
require(ggpubr)
require(gridExtra)
require(GenomicRanges)
require(ggbio)

seed <- 123
set.seed(seed)

p <- 20
q <- 1000
df <- data.frame("snp" = 1:p, 
                 "hotspot_size" = q*rbeta(p, shape1 = 0.1, shape2 = 100))

col <- "grey30"
col_light <- "grey50"
logo_top <- ggdotchart(df, x = "snp", y = "hotspot_size",
                       color = ifelse(df$hotspot_size>1, "#E6AB02", "#66A61E"),
                       sorting = "none",
                       add = "segments",
                       add.params = list(color = col),
                       rotate = FALSE,
                       dot.size = 1,
                       ggtheme = theme_void()) + 
  theme_transparent() + 
  theme(legend.position = "none") 


seed <- 111
set.seed(seed)

N <- 10

gr <- GRanges(seqnames = rep("chr1", N),
              IRanges(
                start = sample(1:(20*N), size = N, replace = TRUE),
                width = sample(70:75, size = N, replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N,
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"),
                              size = N, replace = TRUE),
              pair = sample(letters, size = N,
                            replace = TRUE)) 

logo_bottom <- ggplot() + geom_alignment(gr) + theme_void() + scale_y_reverse()

dir.create("man/figures/", showWarnings = FALSE)

sticker(grid.arrange(logo_top, logo_bottom, heights=c(3.5,1)), 
        package="epispot", 
        p_size=4.5, 
        s_x=0.975, 
        s_y=1, 
        s_width=1.7, 
        s_height=1.3,
        p_x = 0.57, 
        p_y = 1.4, 
        u_color = "white", 
        u_size = 1,
        h_fill="grey65", 
        h_color="grey65",
        filename="man/figures/epispot_logo.png",
        dpi = 1200)

