
library(karyoploteR)
library(dplyr)

directory = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/"
directory_output = paste0(directory, "exp/karyoploteR_plots/")

color_code = c(
  "CTRL" = "#009E73",
  "DON" = "#0072B2",
  "2DG" = "#E69F00",
  "AOA" = "#D30707"
)
chrom_gr = data.frame(seqnames = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                                   "chr7", "chr8", "chr9", "chr10", "chr11", 
                                   "chr12", "chr13", "chr14", "chr15", "chr16", 
                                   "chr17", "chr18", "chr19", "chr20", "chr21", 
                                   "chr22", "chrX", "chrY"),
                      start = rep(0, 24),
                      end = c(248956422, 242193529, 198295559, 190214555,
                              181538259,170805979,159345973,145138636,138394717,
                              133797422, 135086622,133275309,114364328,107043718,
                              101991189,90338345,83257441,80373285,58617616,
                              64444167,46709983,50818468,156040895,57227415))

# https://rdrr.io/bioc/karyoploteR/man/kpPlotBAMDensity.html 
path_bam_2DG = paste0(directory, "data/2DG_possorted_bam.bam")
path_bam_AOA = paste0(directory, "data/AOA_possorted_bam.bam")
path_bam_CTRL = paste0(directory, "data/CTRL_possorted_bam.bam")
path_bam_DON = paste0(directory, "data/DON_possorted_bam.bam")

pp = getDefaultPlotParams(plot.type=1)
pp$ideogramheight = 20
pp$data1height = 500
pp$topmargin = 20

for(i in 21:nrow(chrom_gr)){
  
  print(i)
  chr = chrom_gr$seqnames[i]

  # find max value of peak to put the same axe on all coverages
  kp <- plotKaryotype(genome = "hg38", chromosomes = chr)
  temp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 2e7)
  ymax_axis = temp$latest.plot$computed.values$max.density
  temp <- kpPlotBAMDensity(kp, data=path_bam_AOA, window.size = 2e7)
  if (ymax_axis < temp$latest.plot$computed.values$max.density) {ymax_axis = temp$latest.plot$computed.values$max.density}
  temp <- kpPlotBAMDensity(kp, data=path_bam_DON, window.size = 2e7)
  if (ymax_axis < temp$latest.plot$computed.values$max.density) {ymax_axis = temp$latest.plot$computed.values$max.density}
  temp <- kpPlotBAMDensity(kp, data=path_bam_CTRL, window.size = 2e7)
  if (ymax_axis < temp$latest.plot$computed.values$max.density) {ymax_axis = temp$latest.plot$computed.values$max.density}
  rm(kp)

  # draw coverage
  png(filename = paste0(directory_output, "coverage_kpbam_", chr, ".png"),
      width = 1000, height = 1400)
  
  kp <- plotKaryotype(genome = "hg38", chromosomes = chr, plot.type=1, plot.params = pp)
  kpAddBaseNumbers(kp, tick.dist=10000000, add.units=TRUE, cex=0.5)
  
  kp <- kpPlotBAMDensity(kp, data=path_bam_2DG, window.size = 5e5, col = color_code["2DG"], r0=0,  r1=0.23)
  kp <- kpPlotBAMDensity(kp, data=path_bam_AOA, window.size = 5e5, col = color_code["AOA"], r0=0.27, r1=0.48)
  kp <- kpPlotBAMDensity(kp, data=path_bam_DON, window.size = 5e5, col = color_code["DON"], r0=0.52,  r1=0.73)
  kp <- kpPlotBAMDensity(kp, data=path_bam_CTRL, window.size = 5e5, col = color_code["CTRL"], r0=0.77,  r1=1)
  
  kpAddLabels(kp, "2DG", label.margin = 0.01, side="right", col = color_code["2DG"], r0=0,  r1=0.15, cex=1)
  kpAddLabels(kp, "AOA", label.margin = 0.01, side="right", col = color_code["AOA"], r0=0.27, r1=0.4, cex=1)
  kpAddLabels(kp, "DON", label.margin = 0.01, side="right", col = color_code["DON"], r0=0.52,  r1=0.6, cex=1)
  kpAddLabels(kp, "CTRL", label.margin = 0.01, side="right", col = color_code["CTRL"], r0=0.77,  r1=0.9, cex=1)
  
  kpAxis(kp, ymax=ymax_axis, r0=0,  r1=0.23, cex=0.8)
  kpAxis(kp, ymax=ymax_axis, r0=0.27, r1=0.48, cex=0.8)
  kpAxis(kp, ymax=ymax_axis, r0=0.52,  r1=0.73, cex=0.8)
  kpAxis(kp, ymax=ymax_axis, r0=0.77,  r1=1, cex=0.8)
   
  dev.off()
  
}


