# From a tibble data, with 2 variables varx adn vary, = an extra label, 
# this function calculates "type" statistic and produces a scatterplot 
# with this stat

scatterplot_with_stats <- function(data, varx, vary, label = FALSE, type, title){
  correlation <- rcorr(data[[varx]], data[[vary]], type = type)
  correlation_txt <- paste0(type, ": R = ", round(correlation$r["x","y"], 2),
                            ", pval = ", round(correlation$P["x","y"], 4))
  ggplot(data) +
    aes_string(x = varx, y = vary, label = label) +
    geom_point() +
    geom_smooth(method=lm) +
    theme_bw() +
    geom_text_repel() +
    labs(title=title, subtitle = correlation_txt)
}