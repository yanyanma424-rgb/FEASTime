# visualisation.R

utils::globalVariables(c("Stage", "Contribution", "Source"))

#' Plot cumulative source contributions as stacked bar chart
#'
#' @param feastime_result Output from run_feastime().
#' @param stage_order Optional character vector controlling x-axis order.
#' @param title Plot title.
#' @param palette RColorBrewer palette name (default "Set3").
#' @return A ggplot object.
#' @export
plot_contributions <- function(feastime_result,
                                stage_order = NULL,
                                title = "Cumulative source contributions",
                                palette = "Set3") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required. Install with: install.packages(ggplot2)")
  cumul <- feastime_result$cumulative_contributions
  if (is.null(stage_order)) stage_order <- names(cumul)
  df <- do.call(rbind, lapply(stage_order, function(st) {
    v <- cumul[[st]]
    data.frame(Stage = st, Source = names(v),
               Contribution = as.numeric(v), stringsAsFactors = FALSE)
  }))
  df$Stage <- factor(df$Stage, levels = stage_order)
  ggplot2::ggplot(df, ggplot2::aes(x = Stage, y = Contribution, fill = Source)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
    ggplot2::scale_fill_brewer(type = "qual", palette = palette) +
    ggplot2::labs(title = title, x = "Stage",
                  y = "Cumulative contribution", fill = "Source") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"))
}

#' Plot source contribution trends as line chart
#'
#' @inheritParams plot_contributions
#' @return A ggplot object.
#' @export
plot_trends <- function(feastime_result,
                         stage_order = NULL,
                         title = "Source contribution trends",
                         palette = "Set3") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.")
  cumul <- feastime_result$cumulative_contributions
  if (is.null(stage_order)) stage_order <- names(cumul)
  df <- do.call(rbind, lapply(stage_order, function(st) {
    v <- cumul[[st]]
    data.frame(Stage = st, Source = names(v),
               Contribution = as.numeric(v), stringsAsFactors = FALSE)
  }))
  df$Stage <- factor(df$Stage, levels = stage_order)
  ggplot2::ggplot(df, ggplot2::aes(x = Stage, y = Contribution,
                                   color = Source, group = Source)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_brewer(type = "qual", palette = palette) +
    ggplot2::labs(title = title, x = "Stage",
                  y = "Cumulative contribution", color = "Source") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"))
}
