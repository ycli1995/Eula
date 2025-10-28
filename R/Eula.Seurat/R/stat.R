
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


