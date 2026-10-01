# Generates figures/threshold_demo.png: the neuron's response to stimulus
# currents below, near, and above threshold.
source("hodgkin_huxley.R")
library(ggplot2)

currents <- c(2.0, 2.2, 2.3, 5, 15)
sims <- do.call(rbind, lapply(currents, function(I0) {
  s <- simulate_hh(I0)
  s$label <- sprintf("%.1f uA/cm^2  (%d spike%s)", I0, count_spikes(s),
                     ifelse(count_spikes(s) == 1, "", "s"))
  s
}))
sims$label <- factor(sims$label, levels = unique(sims$label))

p <- ggplot(sims, aes(time, V)) +
  annotate("rect", xmin = 10, xmax = 40, ymin = -Inf, ymax = Inf,
           fill = "#F2C572", alpha = 0.25) +
  geom_line(colour = "#0B6E6E", linewidth = 0.7) +
  facet_wrap(~label, ncol = 1) +
  labs(
    title = "Hodgkin-Huxley neuron: all-or-nothing firing",
    subtitle = "Shaded band = stimulus current step (30 ms).\nThe firing threshold is about 2.24 uA/cm^2: 2.2 gives no spike, 2.3 gives one.",
    x = "Time (ms)", y = "Membrane potential (mV)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"))

ggsave("figures/threshold_demo.png", p, width = 7, height = 9, dpi = 150, bg = "white")
cat("Saved figures/threshold_demo.png\n")
