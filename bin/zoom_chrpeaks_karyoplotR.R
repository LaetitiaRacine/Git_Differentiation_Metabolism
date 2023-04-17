
library(karyoploteR)
library(dplyr)

directory = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/"

path_bam_2DG = paste0(directory, "data/2DG_possorted_bam.bam")

pp = getDefaultPlotParams(plot.type=1)
pp$ideogramheight = 20
pp$data1height = 500
pp$topmargin = 20

chr = "chr2"
kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp)
kpAddBaseNumbers(kp, tick.dist=10000000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)
rm(kp)

kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp, zoom=toGRanges("chr2:32600000-33100000"))
kpAddBaseNumbers(kp, tick.dist=100000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)


chr = "chr21"
kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp)
kpAddBaseNumbers(kp, tick.dist=10000000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)
rm(kp)

kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp, zoom=toGRanges("chr21:5000000-9000000"))
kpAddBaseNumbers(kp, tick.dist=100000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)


chr = "chrY"
kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp)
kpAddBaseNumbers(kp, tick.dist=10000000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)
rm(kp)

kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp, zoom=toGRanges("chrY:56200000-57200000"))
kpAddBaseNumbers(kp, tick.dist=100000, add.units=TRUE, cex=0.5)
kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)
