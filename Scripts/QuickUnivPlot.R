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
           label = ifelse(-log10(p.value)>50|OR<0.7,
                          Biomarkers, "")))+
  geom_point()+geom_label_repel()+ylim(0,NA)+
  theme_bw()
