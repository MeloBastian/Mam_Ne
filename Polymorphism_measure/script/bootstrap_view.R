library(ggplot2)

a=read.table("/home/mbastian/data/Enard_postbusco/144spsummary_genomique_withn50",header=T) #table with genome coverage

b=read.delim("~/data/Enard_postbusco/bootstrap_ps/6002genes/bootstrap_q5%_psandpnps_6002genes_masked")
b$qual <- ifelse(
  b$sp %in% c("Muscardinus_avellanarius", "Beatragus_hunteri", "Sigmodon_hispidus", "Diceros_bicornis", "Sousa_chinensis", "Alouatta_palliata"), "red",
  ifelse(
    b$sp %in% c("Acomys_cahirinus", "Przewalskium_albirostris", "Mastomys_coucha", "Litocranius_walleri", "Cheirogaleus_medius", "Cephalophus_harveyi"), "blue",
    "black"
  )
)


brem<-b[b$qual != "blue",] #"blue" species are the one with the polymorphism systematically removed from the study
c=merge(a,brem, by="sp" )


ggplot(c) +
  geom_errorbar(aes(ymin = q025_pNpS, ymax = q975_pNpS, x = pS, color = coverage)) +
  geom_errorbarh(aes(xmin = q025_pS, xmax = q975_pS, y = pNpS, color = coverage)) +
  geom_point(data = subset(c, qual == "red"), aes(pS, pNpS), shape = 15, size = 2) +
  geom_point(data = subset(c, qual == "black"), aes(pS, pNpS), shape = 1, size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10", limits = c(1e-06, 8e-03)) +
  scale_color_viridis_c(
    trans = "log", 
    breaks = c(min(c$coverage), max(c$coverage)),  
    labels = round(c(min(c$coverage), max(c$coverage)), 2), 
    option = "H"
  ) +
  labs(title = "", 
       x = expression(log[10](Pi[S])),  
       y = expression(log[10](Pi[N]/Pi[S])),  
       color = "Sequencing Depth") +
  theme_minimal()
