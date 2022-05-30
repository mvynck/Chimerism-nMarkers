#######################################
#
# Load libraries
#
#######################################
library(ggplot2)
library(scales)
library(cowplot)
library(reshape2)
library(qvalue)

#######################################
#
# Define general parameters
#
#######################################

#number of simulations per setting
nsim <- 10000

# maximum number of markers in assay
max.markers <- 100


#######################################
#
# Setting 1:
# MAFs sampled from U[0.2, 0.4]
# Autosomes (i.e. no X/Y)
#
#######################################
source("simMAF0.2_0.4.R")

#######################################
#
# Setting 2:
# MAFs all 0.5 ("ideal scenario")
# Diploid (i.e. no X/Y)
#
#######################################
source("simMAF0.5.R")

#######################################
#
# Setting 3:
# MAFs sampled from U[0.3, 0.5]
# Diploid (i.e. no X/Y)
#
#######################################
source("simMAF0.3_0.5.R")

#######################################
#
# Setting 4:
# MAFs from Devyser assay
# Diploid (i.e. no X/Y)
#
#######################################
source("simDevyser.R")

#######################################
#
# Setting 5:
# MAF 0.5
# X chromosome (donor-recipient)
#  sibling: brother-brother, sister-sister, brother-sister
#  or unrelated: female-female, female-male, male-male
# Y chromosome, male-male
#   i.e. male-male (unrelated)
#   
# combined fully informative and potentially informative
#   markers considered
#
#######################################
source("simXY.R")

#######################################
#
# Setting 6:
# MAF 0.5, double transplant
#   
# fully informative markers
#
#######################################
source("simDoubleTransplant.R")

#######################################
#
# Plot:
#   Relation numer of markers - 
#       informativity rate
#
#######################################


# data frame for ggplot of ideal scenario
#  marker type / relatedness impact
df1 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related.ideal,
                  Related = "Sibling",
                  Type = "Informative")
df2 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated.ideal,
                  Related = "Unrelated",
                  Type = "Informative")
df3 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc.ideal,
                  Related = "Parent-child",
                  Type = "Informative")

df4 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related.ideal.homo,
                  Related = "Sibling",
                  Type = "Homozygous informative")
df5 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated.ideal.homo,
                  Related = "Unrelated",
                  Type = "Homozygous informative")
df6 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc.ideal.homo,
                  Related = "Parent-child",
                  Type = "Homozygous informative")

df7 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related.ideal.allinf,
                  Related = "Sibling",
                  Type = "Informative and potentially informative")
df8 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated.ideal.allinf,
                  Related = "Unrelated",
                  Type = "Informative and potentially informative")
df9 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc.ideal.allinf,
                  Related = "Parent-child",
                  Type = "Informative and potentially informative")

df.final.ideal <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

# data frame for ggplot of informative markers only
#  MAF and relatedness impact
df1 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related.ideal,
                  Related = "Sibling",
                  MAF = "all 0.5")
df2 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated.ideal,
                  Related = "Unrelated",
                  MAF = "all 0.5")
df3 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc.ideal,
                  Related = "Parent-child",
                  MAF = "all 0.5")

df4 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related,
                  Related = "Sibling",
                  MAF = "range [0.2 to 0.4]")
df5 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated,
                  Related = "Unrelated",
                  MAF = "range [0.2 to 0.4]")
df6 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc,
                  Related = "Parent-child",
                  MAF = "range [0.2 to 0.4]")

df7 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.related.real1,
                  Related = "Sibling",
                  MAF = "range [0.3 to 0.5]")
df8 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.unrelated.real1,
                  Related = "Unrelated",
                  MAF = "range [0.3 to 0.5]")
df9 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.pc.real1,
                  Related = "Parent-child",
                  MAF = "range [0.3 to 0.5]")

df10 <- data.frame(Markers=3:24,
                  Percentage = means.vec.related.devyser,
                  Related = "Sibling",
                  MAF = "Devyser 24-plex assay")
df11 <- data.frame(Markers=3:24,
                  Percentage = means.vec.unrelated.devyser,
                  Related = "Unrelated",
                  MAF = "Devyser 24-plex assay")
df12 <- data.frame(Markers=3:24,
                  Percentage = means.vec.pc.devyser,
                  Related = "Parent-child",
                  MAF = "Devyser 24-plex assay")

df13 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.mud.mud,
                  Related = "Double (2x MUD)",
                  MAF = "all 0.5")
df14 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.sib.sib,
                  Related = "Double (2x sibling)",
                  MAF = "all 0.5")
df15 <- data.frame(Markers=3:max.markers,
                  Percentage = means.vec.mud.sib,
                  Related = "Double (sibling + MUD)",
                  MAF = "all 0.5")

df.comp <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15)

grey_palette <- c("#FCBF49", "#003049", "#D62828")
grey_palette.ext <- c("#000000", "#000000", "#000000", "#FCBF49", "#003049", "#D62828")
g1 <- ggplot(df.final.ideal, aes(x=Markers, y = Percentage, col = Related, type = Type))+
        geom_line(aes(linetype = Type), size =0.5)+
        theme_linedraw()+
        theme(panel.grid.minor = element_line(color = "#DDDDDD",
                                                size = 0.1,
                                                linetype = 2),
              panel.grid.major = element_line(color = "#BBBBBB"),
              legend.position = "top",
              legend.direction = "vertical",
              legend.title = element_blank(),
              panel.border = element_blank(),
              plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))+
        scale_y_continuous(breaks = seq(0,  1.00, .20), minor_breaks = seq(0, 1, 0.05), labels = label_percent(), limits = c(-0.001, 1), expand=c(0, 0))+
        scale_x_continuous(breaks= seq(0, 100, 10),
                     minor_breaks = seq(0, 100, 1),
                     expand = c(0, 0),
                     limits=c(0, 100))+
        xlab("Number of markers analysed")+ylab("Transplants with at least 3 informative markers")+
        scale_colour_manual(values = grey_palette)+
        scale_linetype_manual(values=c("solid", "longdash", "dotted"))
g2 <- ggplot(df.comp, aes(x=Markers, y = Percentage, col = Related, type = MAF))+
  geom_line(aes(linetype = MAF), size =0.5)+
  theme_linedraw()+
  theme(panel.grid.minor = element_line(color = "#DDDDDD",
                                          size = 0.1,
                                          linetype = 2),
        panel.grid.major = element_line(color="#BBBBBB"),
        legend.position = "top",
        legend.direction = "vertical",
        legend.title = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))+
  scale_y_continuous(breaks = seq(0,  1.00, .01), 
                     minor_breaks = seq(0, 1, 0.005), 
                     labels = label_percent(), 
                     limits = c(0.85, 1), 
                     expand=c(0, 0))+
  scale_x_continuous(breaks= seq(0, 100, 10),
                     minor_breaks = seq(0, 100, 1),
                     expand = c(0, 0),
                     limits=c(10, 100))+
  xlab("Number of markers analysed")+ylab("Transplants with at least 3 informative markers")+
  scale_colour_manual(values = grey_palette.ext)+
  scale_linetype_manual(values=c("dotted", "longdash", "solid", "dotdash"))

pdf("../Figures/sim-nmarkers2.pdf", height = 11)
plot_grid(plotlist=list(g2, g1), labels = "AUTO", nrow = 2)
dev.off()




#######################################
#
# Compare theoretical and empirical
#   results
# Graphical: density plots
# Statistical: Poisson regression
#
#######################################

# compare theoretical vs experimental
MarkerFrequency <- read.csv("Data/NumberInformativeMarkers.csv", sep = ";")
MarkerFrequency$MUD[MarkerFrequency$MUD=="nee"] <- "Related"
MarkerFrequency$MUD[MarkerFrequency$MUD=="Y"] <- "Unrelated"
MarkerFrequency$MUD[MarkerFrequency$MUD=="y"] <- "Unrelated"
MarkerFrequency$MUD[MarkerFrequency$MUD=="ouder-kind"] <- "Parent-child"
MarkerFrequency <- MarkerFrequency[!c(MarkerFrequency$MUD == "haplo?"),]
MarkerFrequency <- MarkerFrequency[complete.cases(MarkerFrequency),]
dataForDensityInf <- list("Sibling Empirical"=MarkerFrequency$Informative.markers[MarkerFrequency$MUD=="Related"],
                          "Unrelated Empirical"=MarkerFrequency$Informative.markers[MarkerFrequency$MUD=="Unrelated"],
                       "Sibling Theoretical"=rowSums(res.related.devyser),
                       "Unrelated Theoretical"=rowSums(res.unrelated.devyser))
dataForDensityHomo <- list("Rel. Emp. Homo."=MarkerFrequency$Homozygous.informative[MarkerFrequency$MUD=="Related"],
                           "Unrel. Emp. Homo."=MarkerFrequency$Homozygous.informative[MarkerFrequency$MUD=="Unrelated"],
                       "Related Th. Homo."=rowSums(res.related.devyser.homo),
                       "Unrelated Th. Homo"=rowSums(res.unrelated.devyser.homo))
dataForDensityAll <- list("Rel. Emp. All"=MarkerFrequency$Informative.markers[MarkerFrequency$MUD=="Related"]+MarkerFrequency$Potentially.informative.markers[MarkerFrequency$MUD=="Related"],
                          "Unrel. Emp. All"=MarkerFrequency$Informative.markers[MarkerFrequency$MUD=="Unrelated"]+MarkerFrequency$Potentially.informative.markers[MarkerFrequency$MUD=="Unrelated"],
                       "Related Th. All"=rowSums(res.related.devyser.allinf),
                       "Unrelated Th. All"=rowSums(res.unrelated.devyser.allinf))
dataForDensityInf <- melt(dataForDensityInf)
dataForDensityHomo <- melt(dataForDensityHomo)
dataForDensityAll <- melt(dataForDensityAll)

grey_palette <- c("#003049", "#D62828", "#F77F00", "#FCBF49")
gInf <- ggplot(dataForDensityInf,
               aes(x = value, y=L1, fill=L1))+
  geom_violin(adjust=3,
              trim = FALSE,
              color = NA,
              alpha = 0.6)+
  geom_boxplot(outlier.shape = NA,
               alpha = 0.5)+
  theme_minimal()+
  xlab("")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank())+
  ylab("")+
  xlim(0,24)+
  scale_fill_manual(values = grey_palette)
gHomo <- ggplot(dataForDensityHomo,
                aes(x = value, y=L1, fill=L1))+
  geom_violin(adjust=3,
              trim = FALSE,
              color = NA,
              alpha = 0.6)+
  geom_boxplot(outlier.shape = NA,
               alpha = 0.5)+
  theme_minimal()+
  xlab("")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank())+
  ylab("")+
  xlim(0,24)+
  scale_fill_manual(values = grey_palette)
gAll <- ggplot(dataForDensityAll,
               aes(x = value, y=L1, fill=L1))+
  geom_violin(adjust=3,
              trim = FALSE,
              color = NA,
              alpha = 0.6)+
  geom_boxplot(outlier.shape = NA,
               alpha = 0.5)+
  theme_minimal()+
  xlab("Number of informative markers")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_blank())+
  ylab("")+
  xlim(0,24)+
  scale_fill_manual(values = grey_palette)
legend <- get_legend(
  gInf + 
    theme(legend.position = "top")
)
plots <- plot_grid(
  legend,
  gInf + theme(legend.position="none"),
  gHomo + theme(legend.position="none"),
  gAll + theme(legend.position="none"),
  labels = c("", "A", "B", "C"),
  nrow = 4, rel_heights = c(.2, 1, 1, 1)
)

pdf("Figures/empthe.pdf", height = 6)
plots
dev.off()

# compare theoretical to empirical using poisson regression
#   on the observed frequencies of different types of markers

dataForDensityInfPoisson <- dataForDensityInf
dataForDensityInfPoisson$Related <- substr(dataForDensityInfPoisson$L1, 1, 3) == "Rel"
dataForDensityInfPoisson$Theoretical <- unlist(strsplit(dataForDensityInfPoisson$L1, " "))[seq(2, nrow(dataForDensityInfPoisson)*2, 2)] == "Theoretical"
InfReg <- glm(value ~ Related + Theoretical, data = dataForDensityInfPoisson, family = "poisson")
summary(InfReg)

dataForDensityHomoPoisson <- dataForDensityHomo
dataForDensityHomoPoisson$Related <- substr(dataForDensityHomoPoisson$L1, 1, 3) == "Rel"
dataForDensityHomoPoisson$Theoretical <- unlist(strsplit(dataForDensityHomoPoisson$L1, " "))[seq(2, nrow(dataForDensityHomoPoisson)*3, 3)] == "Th."
HomoReg <- glm(value ~ Related + Theoretical, data = dataForDensityHomoPoisson, family = "poisson")
summary(HomoReg)

dataForDensityAllPoisson <- dataForDensityAll
dataForDensityAllPoisson$Related <- substr(dataForDensityAllPoisson$L1, 1, 3) == "Rel"
dataForDensityAllPoisson$Theoretical <- unlist(strsplit(dataForDensityAllPoisson$L1, " "))[seq(2, nrow(dataForDensityAllPoisson)*3, 3)] == "Th."
AllReg <- glm(value ~ Related + Theoretical, data = dataForDensityAllPoisson, family = "poisson")
summary(AllReg)

# save all results
save.image("Data/results.RData")
