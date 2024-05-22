library(ggplot2)



make_expr <- function(group, sample, n_fov, prob=0.01, gene="gene") {
  data.frame(cbind(gene, group, sample, expression=log2(rbinom(n_fov,10000,prob))))
}


data <- bind_rows(
  make_expr("Treated","s1", 3, 0.013),
  make_expr("Treated","s2", 3, 0.013),
  make_expr("Treated","s3", 2, 0.013),
  make_expr("Control","s4", 1),
  make_expr("Control","s5", 2),
  make_expr("Control","s6", 3)
)
data$expression <- as.numeric(data$expression)

ggplot(data, aes(x=sample, fill=group, y=expression)) +
  geom_boxplot() +
  theme_void() +
  scale_fill_manual(values=c("#91bec9",'#e3a274'))



