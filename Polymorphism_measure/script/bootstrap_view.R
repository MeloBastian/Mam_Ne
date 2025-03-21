library(ggplot2)

b=read.delim("~/data/Enard_postbusco/bootstrap_ps/7446genes/bootstrap_quantile_psandpnps")
b=read.delim("~/data/Enard_postbusco/bootstrap_ps/bootsrtap_quantile_psandpnps")
b=read.delim("/home/mbastian/Bureau/data_florianL/bootstrap_quantile_psandpnps_6002genes_masked")
b$qual=ifelse(b$sp %in% c("Canis_lupus", "Cheirogaleus_medius", "Acomys_cahirinus","Przewalskium_albirostris","Cephalophus_harveyi" ), "red", "black")


ggplot(b)+
  geom_point(aes(pS, pN, color=qual ), size=.25)+
  #scale_y_continuous(trans = "log10") +
  #scale_x_continuous(trans = "log10") +
  geom_errorbar(aes(ymin=q1_pN, ymax=q3_pN, x=pS, color=qual))+
  geom_errorbarh(aes(xmin=q1_pS, xmax=q3_pS, y=pN, color=qual))










b1 <- b[order(b$pS, decreasing = T), ]
b1 <- b1[!is.na(b1$pS), ]
b1$sp <- factor(b1$sp, levels = b1$sp)



b2 <- b[order(b$pNpS, decreasing = T), ]
b2 <- b2[!is.na(b2$pNpS), ]
b2$sp <- factor(b2$sp, levels = b2$sp)

b3 <- b[order(b$pN, decreasing = T), ]
b3 <- b3[!is.na(b3$pN), ]
b3$sp <- factor(b3$sp, levels = b3$sp)


ggplot(b1) +
  geom_point(aes(x = sp, y = pS)) +
  geom_errorbar(aes(x = sp, ymin = q1_pS, ymax = q3_pS)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggplot(b2) +
  geom_point(aes(x = sp, y = pNpS)) +
  geom_errorbar(aes(x = sp, ymin = q1_pNpS, ymax = q3_pN.pS)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggplot(b3) +
  geom_point(aes(x = sp, y = pN)) +
  geom_errorbar(aes(x = sp, ymin = q1_pN, ymax = q3_pN)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



#########################################################################
#decoupling pS and pN
c=read.delim("~/data/Enard_postbusco/decoupling_pSpN/summary_quantile123")


c1 <- c[order(c$med_pS, decreasing = T), ]
c1 <- c1[!is.na(c1$med_pS), ]
c1$sp <- factor(c1$sp, levels = c1$sp)

c2 <- c[order(c$med_pNpS, decreasing = T), ]
c2 <- c2[!is.na(c2$med_pNpS), ]
c2$sp <- factor(c2$sp, levels = c2$sp)

c3 <- c[order(c$med_pN, decreasing = T), ]
c3 <- c3[!is.na(c3$med_pN), ]
c3$sp <- factor(c3$sp, levels = c3$sp)

ggplot(c1) +
  geom_point(aes(x = sp, y = med_pS)) +
  geom_errorbar(aes(x = sp, ymin = q1_pS, ymax = q3_pS)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggplot(c2) +
  geom_point(aes(x = sp, y = med_pNpS)) +
  geom_errorbar(aes(x = sp, ymin = q1_pNpS, ymax = q3_pN.pS)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggplot(c3) +
  geom_point(aes(x = sp, y = med_pN)) +
  geom_errorbar(aes(x = sp, ymin = q1_pN, ymax = q3_pN)) +
  scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

plot(log(b$pNpS), log(c$med_pNpS))
