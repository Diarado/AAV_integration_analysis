intSiteDistributionPlot <- function (d, chromosomeLengths, alpha = 1) 
{
  library(GenomicRanges)
  library(ggplot2)
  library(gtools)
  library(dplyr)
  d <- d[, c("start", "seqnames", "subject")]
  d <- suppressWarnings(bind_rows(d, bind_rows(lapply(names(chromosomeLengths)[!names(chromosomeLengths) %in% 
                                                                                 unique(as.character(d$seqnames))], function(x) {
                                                                                   data.frame(start = 1, seqnames = x)
                                                                                 }))))
  d$seqnames <- factor(d$seqnames, levels = mixedsort(names(chromosomeLengths)))
  
  d <- lapply(split(d, d$seqnames), function(x) {
    lines <- data.frame(x = rep(x$start, each = 2), 
                        y = rep(c(0, 1), nrow(x)), 
                        g = rep(1:nrow(x), each = 2), 
                        s = x$subject,
                        seqnames = x$seqnames[1])
    box <- data.frame(boxYmin = 0, boxYmax = 1, boxXmin = 1, 
                      boxXmax = chromosomeLengths[[as.character(x$seqnames[1])]], 
                      seqnames = x$seqnames[1])
    list(lines = lines, box = box)
  })
  sites <- do.call(rbind, lapply(d, "[[", 1))
  boxes <- do.call(rbind, lapply(d, "[[", 2))
  
  
  ggplot() + 
    theme_bw() + 
    scale_color_manual(name = 'Dog', values = colorRampPalette(brewer.pal(12, "Paired"))(6)) +
    geom_line(data = sites, 
              alpha = alpha, 
              size = 0.25,
              aes(x, y, group = g, color = s)) +
    geom_rect(data = boxes, 
              color = "black", 
              size = 0.25,
              alpha = 0,
              mapping = aes(xmin = boxXmin, 
                            xmax = boxXmax, 
                            ymin = boxYmin, 
                            ymax = boxYmax)) + 
    facet_grid(seqnames ~ ., switch = "y") + 
    #scale_x_continuous(expand = c(0, 0)) + 
    labs(x = "Genomic position", y = "") + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), 
          strip.text.y = element_text(size = 12, angle = 180), 
          strip.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

report$genomeMap <- intSiteDistributionPlot(filter(intSites, ! subject %in% c('HO2', 'M12')), chromosomeLengths, alpha = 1)
ggsave(report$genomeMap, file = 'tables_and_figures/chromosome_map.pdf', height = 10, width = 8, units = 'in')

