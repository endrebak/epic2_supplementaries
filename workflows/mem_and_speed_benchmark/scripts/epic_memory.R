library(ggplot2)

df = read.table(snakemake@input[[1]], sep=" ", header=1)
print(head(df))
## df$Size = round(df$Size, 2)

g = ggplot(df, aes(x=Intervals, y=MaxRSSGB)) + geom_line(aes(colour=Software)) + geom_point(aes(colour=Software)) + xlab("Library size") + ylab("Memory usage.") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle("B")
g$layout$l[g$layout$name == "title"] <- 1
ggsave(snakemake@output[[1]], height=3, width=3)
