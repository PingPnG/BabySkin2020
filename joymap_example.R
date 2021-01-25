#https://cran.r-project.org/web/packages/ggridges/vignettes/gallery.html
#https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
#https://stackoverflow.com/questions/45384281/ggjoy-facet-with-ggtree
library(ggplot2movies)

ggplot(movies[movies$year>1912,], aes(x = length, y = year, group = year)) +
  geom_density_ridges(scale = 10, size = 0.25, rel_min_height = 0.03) +
  theme_ridges() +
  scale_x_continuous(limits = c(1, 200), expand = c(0, 0)) +
  scale_y_reverse(
    breaks = c(2000, 1980, 1960, 1940, 1920, 1900),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off")