#Code adapted from https://www.sciencedirect.com/science/article/pii/S2666166721007899

#!/usr/bin/Rscript

# init

library(readr)

library(ggplot2)

library(RColorBrewer)

# plot panel B

df <- read_table("Downloads/Brain_Volume_Measurement_MAGMA/image_tx_subtypes.gsa.out", skip = 4)

df$`-log10P` <- 0 - log10(df$P)

df$VARIABLE <- factor(df$VARIABLE, levels = c("Tx1", "Tx10", "Tx11", "Tx12", "Tx13", "Tx14", "Tx15", "Tx16", "Tx17", "Tx18", "Tx19", "Tx2", "Tx20", "Tx3", "Tx4", "Tx5", "Tx6", "Tx7", "Tx8", "Tx9"))

df$TYPE <- "brain_volume_measurement"

ggplot(df, aes(x = VARIABLE, y = TYPE, size = BETA, fill = `-log10P`)) +
  
  scale_fill_gradient(low = "white",
                      
                      limits = c(0, 5),
                      
                      high = "darkgreen") +
  
  scale_size_area(limits = c(0, 0.5)) +
  
  geom_point(shape = 21) +
  
  xlab("Tx Subtype") +
  
  ylab("") +
  
  theme_bw() +
  
  theme(legend.position = "bottom",
        
        legend.box = "horizontal",
        
        axis.text = element_text(size = 14),
        
        axis.title = element_text(size = 16)) +
  
  scale_x_discrete(labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20))


