library(ggplot2)
library(tidyr)
setwd("E:/PycharmProjects/DRS_viral_RNA_modification_scripts/QC")
# 创建数据框
data <- read.csv("data/mis_heatmap_all.csv")

# 将数据框转换为长格式
data_long <- gather(data, key = Reference_base, value = percentage, -Basecall.base,-group)
data_long$percentage <- round(data_long$percentage, 2)
data_long$Basecall.base <- factor(data_long$Basecall.base, levels = c("A", "C", "G", "T", "delete"), ordered = TRUE)
# 使用 ggplot2 创建热图
plot <- ggplot(data_long, aes(x = Reference_base, y = Basecall.base, fill = percentage)) +
  geom_tile() +
  geom_text(aes(label = paste(percentage, "%")), color = "black", size = 2) +  # 添加文本标签
  scale_fill_gradient(low = "white", high = "#F0616E", breaks = c(0, 2, 5), limits = c(0, 5)) +
  labs(x = "Reference base", y = "Basecalled base")+
  theme_bw()+
  facet_wrap(. ~ group, scales = "free_y",ncol=3)+
  theme(
  strip.background = element_blank(),
  legend.position = 'none',
  text = element_text(size = 12),
  panel.grid.minor =element_blank()
  )


ggsave("plots/all.pdf", plot, width =9, height = 3)