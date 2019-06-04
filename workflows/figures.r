#!/usr/bin/env Rscript

# Denes Turei
# turei.denes@gmail.com


require(dplyr)
require(ggplot2)
require(readr)

d <- read_tsv('data/stats_by_resource_20190604.tsv') %>%
    mutate(
        coverage0 = size_cls0 / omnipath0 * 100,
        coverage1 = size_cls1 / omnipath1 * 100,
        group_coverage0 = size_cls0 / size_parent0 * 100,
        group_coverage1 = size_cls1 / size_parent1 * 100
    ) %>%
    arrange(desc(size_cls0)) %>%
    mutate(
        name_cls0 = factor(name_cls0, levels = unique(name_cls0), ordered = TRUE)
    )

pd <- d %>% filter(
    typ_cls0 == typ_cls1 &
    entity == 'protein' &
    typ_cls0 != 'misc' &
    typ_cls0 != 'small_main' &
    !is.na(name_cls0)
)
cd <- d %>% filter(
    typ_cls0 == typ_cls1 &
    entity == 'complex' &
    typ_cls0 != 'misc'
)

classes <- unique(pd$parent0)

for(parent in classes){
    
    pd_c <- pd %>% filter(parent0 == parent) %>%
        group_by(name_cls0) %>%
        summarize_all(first) %>%
        mutate(
            label0 = ifelse(is.na(label0), 'Total', label0)
        ) %>%
        arrange(desc(size_cls0)) %>%
        mutate(
            label0 = factor(label0, levels = unique(label0), ordered = TRUE)
        )
    
    p <- ggplot(pd_c, aes(x = label0, y = size_cls0)) +
        geom_col(fill = 'black') +
        xlab('Resources') +
        ylab('Number of proteins') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        ) +
        ggtitle(sprintf('Resources: %s', toupper(parent)))
    
    ggsave(sprintf('sizes_%s.pdf', parent), device = cairo_pdf, height = 3, width = 4)
    
}
