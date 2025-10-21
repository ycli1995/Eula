
#' @export statByGroups
statByGroups <- function(x, ...) {
  UseMethod("statByGroups", x)
}

#' @export
#' @method statByGroups data.frame
statByGroups.data.frame <- function(
    x,
    group.by,
    sum = NULL,
    sum.pfx = "Total ",
    mean = NULL,
    mean.pfx = "Mean ",
    median = NULL,
    median.pfx = "Median ",
    n.name = "Number",
    add.total = "Total",
    add.pct = "Proportion",
    ...
) {
  n.name <- n.name[1] %||% "Number"
  stat <- x %>%
    group_by(across(all_of(group.by))) %>%
    summarise(
      !!n.name := n(),
      across(
        all_of(sum),
        ~ sum(., na.rm = TRUE),
        .names = paste0(sum.pfx, "{.col}")
      ),
      across(
        all_of(mean),
        ~ mean(., na.rm = TRUE),
        .names = paste0(mean.pfx, "{.col}")
      ),
      across(
        all_of(median),
        ~ median(., na.rm = TRUE),
        .names = paste0(median.pfx, "{.col}")
      )
    )
  if (validCharacters(add.pct)) {
    stat <- stat %>%
      ungroup() %>%
      mutate(
        across(all_of(n.name), ~ . / sum(.), .names = add.pct),
        .after = all_of(n.name)
      )
  }
  stat <- as.data.frame(stat, check.names = FALSE)
  if (!validCharacters(add.total)) {
    return(stat)
  }
  tot <- x %>%
    summarise(
      across(all_of(group.by), ~ add.total),
      !!n.name := n(),
      across(
        all_of(sum),
        ~ sum(., na.rm = TRUE),
        .names = paste0(sum.pfx, "{.col}")
      ),
      across(
        all_of(mean),
        ~ mean(., na.rm = TRUE),
        .names = paste0(mean.pfx, "{.col}")
      ),
      across(
        all_of(median),
        ~ median(., na.rm = TRUE),
        .names = paste0(median.pfx, "{.col}")
      )
    )
  if (validCharacters(add.pct)) {
    tot <- tot %>%
      ungroup() %>%
      mutate(
        across(all_of(n.name), ~ . / sum(.), .names = add.pct),
        .after = all_of(n.name)
      )
  }
  tot <- as.data.frame(tot, check.names = FALSE)
  stat <- rbind(tot, stat)
  stat
}

#' @export statCrossTab
statCrossTab <- function(x, ...) {
  UseMethod("statCrossTab", x)
}

#' @export
#' @method statCrossTab data.frame
statCrossTab.data.frame <- function(
    x,
    group.by,
    stat.what,
    add.pct = c("none", "paste", "column"),
    ...
) {
  add.pct <- match.arg(add.pct)
  group.by <- group.by[1]
  stat.what <- stat.what[1]
  x[, group.by] <- as.factor(x[, group.by])
  x[, stat.what] <- as.factor(x[, stat.what])
  stat <- x %>%
    group_by(across(all_of(c(group.by, stat.what))), .drop = FALSE) %>%
    summarise(y = n())
  tot <- x %>%
    group_by(across(all_of(group.by)), .drop = FALSE) %>%
    summarise(across({{stat.what}}, ~ factor("Total")), y = n())
  stat <- full_join(x = tot, y = stat) %>%
    dcast(formula = as.formula(paste0(stat.what, " ~ ", group.by)), fill = 0)
  if (add.pct == "paste") {
    stat <- stat %>%
      mutate_if(
        is.numeric,
        list(~paste0(., " (", round(./.[1] * 100,2), "%)"))
      )
  }
  if (add.pct == "column") {
    ll <- levels(x[, group.by])
    stat <- stat %>%
      mutate(
        across({{ll}}, ~ round(./.[1] * 100, 2), .names = paste0("{.col}_pct"))
      )
  }
  return(stat)
}

#' @export featureSetPct
featureSetPct <- function(object, ...) {
  UseMethod("featureSetPct", object)
}

#' @export
#' @method featureSetPct default
featureSetPct.default <- function(object, features, ...) {
  features <- intersect(features, rownames(object))
  if (length(features) == 0) {
    fastWarning("No feature found.")
    return(rep.int(0, ncol(object)))
  }
  csum0 <- Matrix::colSums(object)
  csum <- Matrix::colSums(object[features, , drop = FALSE])
  csum / csum0
}
