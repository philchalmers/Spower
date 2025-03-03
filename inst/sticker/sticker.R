library(pwr)
library(hexSticker)

N <- 250

ds <- seq(0, 2, by=.1)
ps <- sapply(ds, \(d) pwr.t.test(n=10, d=d)$power)
se <- sqrt(ps * (1 - ps) / N )
plot(ds, ps, las=1, type='b', pch=16)

library(ggplot2)
dd <- data.frame(ds, ps, se.l=ps-3*se, se.u=ps+3*se)
ggplot(dd, aes(ds, ps)) +
	geom_ribbon(aes(ymin=se.l, ymax=se.u), fill='grey70') +
	geom_point() + geom_line() +
	xlab("") + ylab("") + #ggtitle("Power") +
	#theme(plot.title = element_text(hjust = 0.5)) +
	theme_void() + theme_transparent() -> gg
gg

# main_col <- "#FDE725"
main_col <- "white"

sticker(gg, package='Spower', filename='S.png',
		s_x = 1,
		s_y = 1.1,
		s_width = 1.4,
		s_height = 1,
		p_x = 1,
		p_y = 0.5,
		p_color = main_col,
		p_size = 80,
		p_family = "sans",
		h_size = 1,
		h_fill = "#212529",
		h_color = main_col,
		dpi = 1200,
		spotlight = TRUE, l_y=1, l_alpha = .3)
