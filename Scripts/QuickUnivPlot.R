library(ggplot)
library(stringr)
pltrds = readRDS("../FULLResults/EDA/bio_dist_impLOG.rds")
data = pltrds$data

data$Biomarker = str_replace_all(data$Biomarker, "\\.", " ")
data$Biomarker = str_replace_all(data$Biomarker, "Glycated haemoglobin HbA1c", "HbA1c")
data$Biomarker = str_replace_all(data$Biomarker, "Gamma glutamyltransferase", "GGT")
data$Biomarker = str_replace_all(data$Biomarker, "Alanine aminotransferase", "Alanine aminotrans.")
data$Biomarker = str_replace_all(data$Biomarker, "Alkaline phosphatase", "Alkaline phosph.")
data$Biomarker = str_replace_all(data$Biomarker, "Aspartate aminotransferase", "Aspartate aminotrans.")
data$CVD_status = plyr::revalue(data$CVD_status, c("0"="Controls", "1"="Cases") )
fig = data %>%
  ggplot(aes(x = log10(Amount), color = CVD_status)) +
  geom_density(alpha = 0)+
  facet_wrap(Biomarker~.,
             scales = "free")+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  ggtitle("Biomarker distribution for imputed data(log scaled)")+
  labs(color = "CVD status")
fig




# quick univariate plots from data
library(readr)
library(ggrepel)
univtable = read_csv("../FULLResults/UnivariateAnalysis.csv")
colnames(univtable)[4:5] = c("lowerCI", "upperCI")
univtable= univtable %>%
  mutate_at(c("lowerCI", "upperCI", "OR"), ~round(., 4))
write_csv(univtable, "../FULLResults/UnivariateAnalysis.csv")
univtable %>% filter(Data == "KNN") %>%
ggplot(aes(y = OR, x = -log10(p.value),
           label = ifelse(-log10(p.value)>50|
                            OR<0.7|
                            (OR>1.05 & -log10(p.value)>40),
                          Biomarkers, "")))+
  geom_point(fill = "grey", color = "black", shape = 21)+
  geom_label_repel(min.segment.length = 0)+
  ylim(0,NA)+
  theme_bw()+xlab(expression("-log10(p-value)"))+
  geom_vline(xintercept = -log10(0.05/28), linetype = "dashed")+
  annotate("text", x = -log10(0.05/28)+3, y = 2.8,
            label = paste0("Bonferroni\nadjusted\nsignificance\nlevel: p =",
                           round(0.05/28,4)),
           hjust = "inward", vjust = "inward")+
  ylab("Odds Ratio(OR)")+theme(text = element_text(size = 16))

