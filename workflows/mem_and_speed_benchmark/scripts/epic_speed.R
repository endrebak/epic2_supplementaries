library(ggplot2)
library(gridExtra)

df = read.table(snakemake@input[[1]], sep=" ", header=1)
df$Software = factor(df$Software, levels=c("SICER", "MACS2", "epic2"))

f = ggplot(df, aes(x=Intervals, y=Minutes)) + geom_line(aes(colour=Software)) + geom_point(aes(colour=Software)) + xlab("Library Size") + ylab("Running time in minutes.") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + ggtitle("A") # + theme() # + labs(tag = "A")
# theme(legend.position="none")
f$layout$l[f$layout$name == "title"] <- 0

g = ggplot(df, aes(x=Intervals, y=MaxRSSGB)) + geom_line(aes(colour=Software)) + geom_point(aes(colour=Software)) + xlab("Library size") + ylab("Memory usage in GB.") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ggtitle("B")
g$layout$l[g$layout$name == "title"] <- 1

ggsave(snakemake@output[[1]], arrangeGrob(f, g, nrow=1, ncol=2), height=3, width=6)
