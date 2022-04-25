library(ggplot2)
library(ggpubr)
setwd("~/Patera/RADseq/SNAQ/r100-update/")
#ggplot command for R. For some reason it would not show in Julia, but it works in R.

net1<-read.csv("~/Patera/RADseq/SNAQ/r100-update/fittedCF_bs10_net1.csv")
net0<-read.csv("~/Patera/RADseq/SNAQ/r100-update/fittedCF_bs10_net0.csv")
net2<-read.csv("~/Patera/RADseq/SNAQ/r100-update/fittedCF_bs10_net2.csv")
net0_plot<-ggplot(net0, aes(x=obsCF,y=expCF)) + theme_classic() +
  geom_segment(x=0,y=0,xend=1,yend=1, color="#008080", size=0.3) + # diagonal line
  geom_point(alpha=0.5, color="#008080", position=position_jitter(width=0.005, height=0.005)) +
  ylab("quartet CF expected from network") + ggtitle("Network with 0 reticulations") + xlab("quartet CF observed in gene trees") + coord_equal(ratio=1);
net1_plot<-ggplot(net1, aes(x=obsCF,y=expCF)) + theme_classic() +
  geom_segment(x=0,y=0,xend=1,yend=1, color="#008080", size=0.3) + # diagonal line
  geom_point(alpha=0.5, color="#008080", position=position_jitter(width=0.005, height=0.005)) +
  ylab("quartet CF expected from network") + ggtitle("Network with 1 reticulation") + xlab("quartet CF observed in gene trees") + coord_equal(ratio=1);
net2_plot<-ggplot(net2, aes(x=obsCF,y=expCF)) + theme_classic() +
  geom_segment(x=0,y=0,xend=1,yend=1, color="#008080", size=0.3) + # diagonal line
  geom_point(alpha=0.5, color="#008080", position=position_jitter(width=0.005, height=0.005)) +
  ylab("quartet CF expected from network") + ggtitle("Network with 2 reticulation") + xlab("quartet CF observed in gene trees") + coord_equal(ratio=1);

net0_plot
net1_plot
net2_plot
pdf("observed-vs-expected_CF-plots.pdf")
plot_for_saving<-ggarrange(net0_plot, net1_plot + rremove("ylab"), net2_plot + rremove("ylab"), labels = c("A", "B", "C"), ncol = 3, nrow = 1, align = "v")
dev.off()
ggexport(plot_for_saving, filename="observed-vs-expected_CF-plots.pdf")
