
#' @importFrom dplyr full_join group_by mutate n relocate rename_with
#' summarise
#' @importFrom tidyr replace_na
#' @export
StatFilterCells <- function(raw, filtered, group.by = "Samples") {
  df.list <- list(raw = raw, filtered = filtered)
  for (i in names(df.list)) {
    df.list[[i]] <- df.list[[i]] %>%
      dplyr::group_by(across(all_of(group.by))) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        across(
          where(is.numeric) & !matches("n_cells"),
          median,
          .names = "median {.col}"
        )
      ) %>%
      dplyr::rename_with(~ paste(i, .x), where(is.numeric))
  }
  df <- Reduce(dplyr::full_join, df.list) %>%
    dplyr::mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) %>%
    dplyr::relocate(
      all_of(c("raw n_cells", "filtered n_cells")),
      .after = group.by
    ) %>%
    dplyr::mutate(
      pct = paste0(round(`filtered n_cells`/`raw n_cells` * 100, 2) , "%"),
      .after = "filtered n_cells"
    )
  df
}

#' @importFrom dplyr summarize
#' @importFrom ggplot2 geom_bar theme_light
#' @export
StatMarker <- function(markers, colors = NULL, outpfx = "DeGene.stat", ...) {
  stat <- markers %>%
    dplyr::group_by(Clusters) %>%
    dplyr::summarize(Number = n())
  writeTable(stat, paste0(outpfx, ".xls"))

  p <- ggplot(stat, aes(x = Clusters, y = Number, fill = Clusters)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(x = "", y = "Number of DE features") +
    bar_theme_default() +
    theme(legend.position = "none")
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  h <- 5
  w <- nrow(stat) * 0.4 + 1.2
  ggsave(p, file = paste0(outpfx, ".pdf"), height = h, width = w)

  invisible(stat)
}
