# ==============================================================================

# Title: "Diglis Fish Pass Plot"
# Author: "Graham Sellers"
# Date: "19 April 2024"

# ==============================================================================

# read in the data:
dat = read.csv('data/Diglis_Counts_2023.tsv', header = T, sep = '\t', stringsAsFactors = F)

# get position of every other entry:
pos = seq(1, nrow(dat), 2)
# make the list of relevant dates:
dates = sub('/2023', '', dat$Date[pos])

# get entries with DNA reads:
dna = dat[!is.na(dat$Sabrina_eDNA_reads),]  # not NAs but actual numbers!
xpos = as.numeric(rownames(dna)) - 0.5  # nudge the plotting position

# ==============================================================================

# PLOT THE PLOT:

png('output/shad_fish_pass.png', width = 1920, height = 1080, units = 'px', pointsize = 30)

par(mar = c(5,3,2,4))

barp = barplot(dat$Shad,
        col = '#440154',
        border = NA,
        las = 1,
        space = 0,
        xpd = F,
        axes = F,
        ylim = c(0, 45))

# vertical lines and rects for reads:
abline(v = xpos, lty = 3, lwd = 5)
rect(xpos - 0.3, rep(0, length(xpos)), xpos + 0.3, dna$Sabrina_eDNA_reads/100, col = '#70cf57', xpd = T)

# 3 axes added:
axis(side = 1, at = seq(0, length(barp)), labels = NA, tcl = -0.3, pos = 0, lwd = 3)
axis(side = 2, at = seq(0, 40, 5), labels = seq(0, 40, 5), tcl = -0.5, pos = 0, las = 1, lwd = 3)
axis(side = 4, at = seq(0, 40, 5), labels = seq(0, 40, 5)*100, tcl = -0.5, pos = length(barp), las = 1, lwd = 3)

# dates for x axis:
text(barp[pos], rep(-3, length(pos)), labels = dates, srt = 90, cex = 1, xpd = T)
# dna read counts per sample:
text(xpos, rep(46, length(pos)), labels = dna$Sabrina_eDNA_reads, cex = 1, xpd = T)

# axis labels:
mtext('Diglis fish pass count', side = 2, at = 20, line = 0.7, cex = 1.5)
mtext('Shad eDNA read count', side = 4, at = 20, line = 2, cex = 1.5)
mtext('Date', side = 1, line = 3.5, cex = 1.5)

dev.off()


