
require(dplyr)
require(ggplot2)
require(readr)

maincov <- read_tsv('main_coverage.tsv') %>%
    mutate(cov = omnipath / total * 100) %>%
    arrange(desc(total)) %>%
    mutate(cls = factor(cls, levels = cls, ordered = TRUE))

p <- ggplot(
    maincov, aes(x = cls, total)
) +
geom_col(fill = 'black') +
xlab('Classes') +
ylab('Number of proteins') +
theme_bw() +
theme(
    text = element_text(family = 'DINPro'),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
)

ggsave('class_sizes.pdf', device = cairo_pdf, width = 5, height = 4)



p <- ggplot(
    maincov, aes(x = cls, cov)
    ) +
geom_col(fill = 'black') +
xlab('Classes') +
ylab('Covered in OmniPath') +
theme_bw() +
theme(
    text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        )

ggsave('class_coverage.pdf', device = cairo_pdf, width = 5, height = 4)



cat_overlaps <- read_tsv('category_overlaps.tsv') %>%
    mutate(poverlap = overlap / pmin(size0, size1) * 100) %>%
    filter(cat0 != cat1)

p <- ggplot(cat_overlaps, aes(x = cat0, y = cat1, size = poverlap)) +
    geom_point(color = 'black') +
    scale_radius(
        guide = guide_legend(title = 'Overlap (%)'),
        range = c(0, 6)
    ) +
    xlab('Categories') +
    ylab('Categories') +
    theme_bw() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave('category_overlaps.pdf', device = cairo_pdf, width = 6, height = 5)


dens <- 0.00040319772927432604
connections <- read_tsv('connections.tsv')

connections <- bind_rows(
    connections,
    connections %>% rename(cat0 = cat1, cat1 = cat0)
) %>%
mutate(econn = (conn / (size0 * size1)) / dens)

p <- ggplot(connections, aes(x = cat0, y = cat1, size = conn)) +
    geom_point(color = 'black') +
    scale_radius(
        guide = guide_legend(title = 'Connections'),
        range = c(0, 6)
    ) +
    xlab('Categories') +
    ylab('Categories') +
    theme_bw() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave('numof_connections.pdf', device = cairo_pdf, width = 6, height = 5)


p <- ggplot(connections, aes(x = cat0, y = cat1, size = econn)) +
    geom_point(color = 'black') +
    scale_radius(
        guide = guide_legend(title = 'Enrichment'),
        range = c(0, 6)
    ) +
    xlab('Categories') +
    ylab('Categories') +
    theme_bw() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave('enrich_connections.pdf', device = cairo_pdf, width = 6, height = 5)

