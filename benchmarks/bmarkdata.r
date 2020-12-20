#!/usr/bin/env Rscript

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://omnipathdb.org/
#

require(magrittr)
require(dplyr)
require(readr)
require(stringr)
require(ggplot2)
require(OmnipathR)

options(
    omnipath.print_urls = TRUE,
    width = 200,
    dplyr.width = 200
)
datadir <- 'data'
figdir <- 'figures'
intercell_variants <- c(
    'omnipath-all',
    'omnipath-q50',
    'omnipath-ligrec-only',
    'omnipath-ligrec-only-q50',
    'cellphonedb',
    'ramilowski'
)
intercell_resources <- c(
    'CellPhoneDB',
    'Ramilowski2015',
    'Baccin2019',
    'LRdb',
    'Kirouac2010',
    'ICELLNET',
    'iTALK',
    'EMBRACE',
    'HPMR',
    'Guide2Pharma'
)


.in_datadir <- function(fname){

    file.path(datadir, fname)

}


.pipe_message <- function(insert, msg){

    message(sprintf(msg, insert))
    flush.console()
    insert

}


retrieve_intercell <- function(
        causality,
        qthreshold = NULL,
        only = NULL,
        resource = NULL
    ){

    source <- `if`(is.null(resource), 'composite', 'resource_specific')

    import_omnipath_intercell(
        causality = causality,
        scope = 'generic',
        source = source,
        entity_type = 'protein',
        resources = resource
    ) %>%
    as_tibble() %>%
    {`if`(
        is.null(qthreshold),
        .,
        group_by(., parent) %>%
        filter(
            consensus_score >=
            quantile(consensus_score, qthreshold)
        ) %>%
        ungroup()
    )} %>%
    {`if`(
        is.null(only),
        .,
        filter(., parent == only)
    )}

}


intercell_param <- list(
    transmitter_all = list(
        causality = 'transmitter'
    ),
    transmitter_q50 = list(
        causality = 'transmitter',
        qthreshold = .5
    ),
    transmitter_ligrec_only = list(
        causality = 'transmitter',
        only = 'ligand'
    ),
    transmitter_ligrec_only_q50 = list(
        causality = 'transmitter',
        only = 'ligand',
        qthreshold = .5
    ),
    transmitter_cellphonedb = list(
        causality = 'transmitter',
        resource = 'CellPhoneDB'
    ),
    transmitter_ramilowski = list(
        causality = 'transmitter',
        resource = 'Ramilowski2015'
    ),
    receiver_all = list(
        causality = 'receiver'
    ),
    receiver_q50 = list(
        causality = 'receiver',
        qthreshold = .5
    ),
    receiver_ligrec_only = list(
        causality = 'receiver',
        only = 'receptor'
    ),
    receiver_ligrec_only_q50 = list(
        causality = 'receiver',
        only = 'receptor',
        qthreshold = .5
    ),
    receiver_cellphonedb = list(
        causality = 'receiver',
        resource = 'CellPhoneDB'
    ),
    receiver_ramilowski = list(
        causality = 'receiver',
        resource = 'Ramilowski2015'
    )
)


compile_intercell <- function(){

    intercell_all <- NULL

    for(dataset in names(intercell_param)){

        param <- intercell_param[[dataset]]
        database <- `if`(is.null(param$resource), 'omnipath-', '')
        label_parts <- dataset %>%
            str_split_fixed('_', n = 2)
        causality <- label_parts[1]
        variant <- label_parts[2] %>% str_replace_all('_', '-')
        variant <- sprintf('%s%s', database, variant)
        outfile <- sprintf('intercell_%s_%s.tsv', causality, variant) %>%
            .in_datadir()

        message(sprintf('Writing intercell data to %s', outfile))
        flush.console()

        intercell_all <- do.call(retrieve_intercell, param) %>%
            select(uniprot, genesymbol) %>%
            group_by(uniprot) %>%
            summarize_all(first) %>%
            ungroup() %>%
            write_tsv(outfile) %>%
            mutate(
                causality = causality,
                variant = variant
            ) %>%
            bind_rows(intercell_all)

    }

    intercell_all %>% write_tsv(.in_datadir('intercell_all.tsv'))

    invisible(NULL)

}


retrieve_networks <- function(){

    tf_target_all <- import_transcriptional_interactions(
            dorothea_levels = c('A', 'B', 'C', 'D'),
            entity_types = 'protein'
        ) %>%
        write_tsv(.in_datadir('tftarget_omnipath-all.tsv'))

    tf_target_ab <- import_transcriptional_interactions(
            entity_types = 'protein'
        ) %>%
        write_tsv(.in_datadir('tftarget_omnipath-hiconf.tsv'))

    ppi_omnipath <- import_post_translational_interactions(
            entity_types = 'protein'
        ) %>%
        write_tsv(.in_datadir('ppi_omnipath-all.tsv'))

    invisible(NULL)

}


compile_networks <- function(){

    for(infile in list.files(datadir)){

        if(!any(startsWith(infile, c('ppi_', 'tftarget_')))){

            next

        }

        network_variant <- infile %>%
            str_replace('\\..*$', '') %>%
            str_split_fixed('_', n = 2) %>%
            `[`(2)


        network_type <- `if`(
            startsWith(infile, 'ppi_'),
            'intracell',
            'tftarget'
        )

        network_path <- 'network__%s__%s.tsv' %>%
            sprintf(network_type, network_variant) %>%
            .in_datadir()

        message(sprintf('Writing network to %s', network_path))
        flush.console()

        this_network <- read_network(infile) %>%
            write_tsv(network_path)

        if(startsWith(infile, 'ppi_')){

            intercell_path_trunk <- 'network__intercell__%s' %>%
                sprintf(network_variant) %>%
                .in_datadir()

            for(ic_var in intercell_variants){

                intercell_path <- '%s__%s.tsv' %>%
                    sprintf(intercell_path_trunk, ic_var)

                transmitters <- intercell_components('transmitter', ic_var)
                receivers <- intercell_components('receiver', ic_var)

                this_network %>%
                filter(
                    source %in% transmitters &
                    target %in% receivers
                ) %>%
                write_tsv(intercell_path)

                message(sprintf(
                    'Writing intercell network to %s',
                    intercell_path
                ))
                flush.console()


            }

        }

    }

    invisible(NULL)

}


read_network <- function(infile){

    this_network <- infile %>%
        .in_datadir() %>%
        .pipe_message('Reading network data from %s') %>%
        read_tsv(col_types = cols()) %>%
        filter(
            as.logical(is_directed)
        ) %>%
        {`if`(
            'entity_type_source' %in% colnames(.),
            filter(
                .,
                entity_type_source == 'protein' &
                entity_type_target == 'protein'
            ),
            .
        )}

    this_network %>%
        filter(as.logical(is_inhibition)) %>%
        mutate(sign = -1) %>%
        bind_rows(
            this_network %>%
            filter(
                !as.logical(is_inhibition) | as.logical(is_stimulation)
            ) %>%
            mutate(sign = is_stimulation)
        ) %>%
        select(
            source,
            target,
            source_genesymbol,
            target_genesymbol,
            sign
        ) %>%
        group_by(source, target, sign) %>%
        summarize_all(first) %>%
        ungroup()

}


intercell_components <- function(causality, variant){

    'intercell_%s_%s.tsv' %>%
    sprintf(causality, variant) %>%
    .in_datadir() %>%
    .pipe_message('Reading intercell data from %s') %>%
    read_tsv(col_types = cols()) %>%
    pull(uniprot) %>%
    unique()

}


textbook_benchmark <- function(){

    dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
    .intercell_resources <- sort(intercell_resources)
    cov_data <- NULL
    totals <- NULL

    icn_textbook <- 'cytokine_receptor_book_translated.tsv' %>%
        .in_datadir() %>%
        read_tsv(col_types = cols()) %>%
        filter(ligand_uniprot != '' & receptor_uniprot != '') %>%
        group_by(ligand_uniprot, receptor_uniprot) %>%
        summarize_all(first) %>%
        ungroup()

    for(resource in .intercell_resources){

        param <- list(resources = resource)

        this_icn <- import_intercell_network(
                interactions_param = param,
                transmitter_param = param,
                receiver_param = param
            ) %>%
            group_by(source, target) %>%
            summarize_all(first) %>%
            ungroup()

        totals %<>% append(this_icn %>% nrow)

        cov_data %<>%
            append(
                this_icn %>%
                intercell_coverage(icn_textbook)
            )

    }

    icn_total <- import_intercell_network(
            interactions_param = list(
                datasets = c('omnipath', 'ligrecextra', 'pathwayextra')
            )
        ) %>%
        group_by(source, target) %>%
        summarize_all(first) %>%
        ungroup()

    totals %<>% append(icn_total %>% nrow)

    cov_data %<>%
        append(
            icn_total %>%
            intercell_coverage(icn_textbook)
        )

    coverages <- tibble(
        resource = c(.intercell_resources, 'OmniPath'),
        coverage = cov_data / nrow(icn_textbook) * 100,
        size = totals
    ) %>%
    arrange(desc(coverage)) %>%
    mutate(resource = factor(resource, levels = resource, ordered = TRUE))

    plot_textbook_benchmark(coverages)
    plot_textbook_benchmark_bar(coverges)

    invisible(coverages)

}


plot_textbook_benchmark <- function(coverages){

    cairo_pdf(
        'textbook_coverage.pdf',
        width = 3,
        height = 3,
        family = 'DINPro'
    )

        {
            ggplot(coverages, aes(x = coverage, y = resource)) +
            geom_point(size = 4, color = '#4268B3') +
            geom_text(aes(
                x = coverage - 10,
                y = resource,
                label = sprintf('%.1f%%', coverage)),
                size = 2
            ) +
            scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
            expand_limits(x = c(18, 102)) +
            xlab('Coverage [%]') +
            ylab('Resources') +
            ggtitle(paste0(
                'Coverage of intercellular communi-\ncation resources ',
                'on a textbook\nknowledge of 131 interactions'
            )) +
            theme_bw() +
            theme(
                plot.title = element_text(size = 8)
            )
        } %>% print()

    dev.off()

}


plot_textbook_benchmark_bar <- function(coverages){

    coverages %<>% mutate(
            size = ifelse(resource == 'HPMR', 650, size)
        ) %>%
        arrange(desc(size)) %>%
        mutate(
            cov_label = sprintf('%i%%', round(coverage)),
            resource = factor(
                resource,
                levels = unique(resource),
                ordered = TRUE
            )
        )

    cairo_pdf(
        'textbook_coverage_bar.pdf',
        width = .8,
        height = 3,
        family = 'DINPro'
    )

        {
            ggplot(coverages, aes(x = coverage, y = resource)) +
            geom_col(fill = '#4268B3') +
            scale_x_continuous(breaks = c(0, 50, 100)) +
            scale_y_discrete(
                labels = coverages$cov_label,
                position = 'right'
            ) +
            theme_minimal() +
            theme(
                plot.title = element_blank(),
                text = element_text(
                    family = 'HelveticaNeueLT Std Lt Cn',
                    face = 'bold'
                ),
                axis.title = element_blank()
            )
        } %>% print()

    dev.off()

}


intercell_coverage <- function(icn, icn_textbook){

    icn %>%
    select(source, target) %>%
    inner_join(
        icn_textbook,
        by = c(
            'source' = 'ligand_uniprot',
            'target' = 'receptor_uniprot'
        )
    ) %>%
    group_by(source, target) %>%
    summarize_all(first) %>%
    ungroup() %>%
    nrow()

}


bmarkdata_main <- function(){

    compile_intercell()
    retrieve_networks()
    compile_networks()
    textbook_benchmark()

}