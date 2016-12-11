#####
# You can use this script to delineate super-enhancers from a .bed, *only* after it has been
# pre-processed by bamliquidator. Conveniently this can be done by riesling.py, which also
# allows for alternative stratifications of super-enhancers.
#
#
# Run as: Rscript get-SuperEnhancers.R input.bed output_directory
#
#
# **Input: A .bed from bamliquidator
#    (A .bed with 7 columns, where the 7th column contains the normalized score for ranking.)
#
#  ## This input might look like:
#      V1      V2      V3 V4    V5 V6         V7
#   1 chr1 3514643 3515351  1 30438  . 0.17688850
#   2 chr1 4426753 4427110  2  4119  . 0.04747231
#
#
# **Output:
#    Super-enhacners (defined by the point where the tangent reaches 1) and
#    Stretch enhancers (>3 kb)
#    Super-Stretch enhancers (> tangent cutoff AND >3kb)
#
#
# Copyright (c) 2014-2016 Nick Semenkovich <semenko@alum.mit.edu>.
#   https://nick.semenkovich.com/
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
#   https://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  http://opensource.org/licenses/MIT
#
# Source: https://github.com/GordonLab/riesling-pipeline
#
# Includes approaches inspired by:
#  https://stackoverflow.com/questions/29642867/drawing-a-tangent-to-the-plot-and-finding-the-x-intercept-using-r

# Force sane graphics output
options(bitmapType='cairo')

#### Cleanup
# dev.off(dev.list()["RStudioGD"]) # Clear the graphs
rm(list = ls(all = TRUE)) # Clear all ojects from workspace
cat("\014") # Reset the console

# argv[1] is the input bed file we're operating on. 
args <- commandArgs(TRUE)
input_filename = args[1]
output_dir = args[2]


print(paste('Working on:', input_filename))
print(paste('Output dir:', output_dir))

bed <- read.table(input_filename, sep="\t", header=FALSE)

## Looks like:
#    V1      V2      V3 V4    V5 V6         V7
# 1 chr1 3514643 3515351  1 30438  . 0.17688850
# 2 chr1 4426753 4427110  2  4119  . 0.04747231

if (!is.na(output_dir)) {
  print(paste('Current directory is:', getwd()))
  print(paste('Setting output directory to:', output_dir))
  setwd(output_dir)
}


## Sort and scale axes
y = sort(bed[,c(7)]*(bed[,c(3)]-bed[,c(2)]))  # normalized_counts * width
x = c(1:length(y))
ynorm = y*length(x)/max(y)
# plot(x, ynorm)
# plot(x, ynorm, log="xy")

spl <- smooth.spline(x, ynorm)
# pred <- predict(spl)
# lines(pred, col=2)

ynorm.prime <- diff(ynorm)/diff(x)
# plot(ynorm.prime)
pred.prime <- predict(spl, deriv=1)
# lines(pred.prime$y, col=2)

## Find where the tangent line first crosses 1
se_cutoff <- min(which(pred.prime$y > 1))
print(paste('Inflection at entry:', se_cutoff))
print(paste('Corresponding cutoff score:', y[[se_cutoff]]))

# Use the spline models to plot tangent to that point
pred0 <- predict(spl, x=se_cutoff, deriv=0)
pred1 <- predict(spl, x=se_cutoff, deriv=1)

# And compute intercepts for graphing
yint <- pred0$y - (pred1$y*se_cutoff)
xint <- -yint/pred1$y


################ Save subpopulation beds
y_se_cutoff = y[[se_cutoff]]
se_population = bed[bed$V7*(bed$V3-bed$V2) >= y_se_cutoff,]
te_population = bed[bed$V7*(bed$V3-bed$V2) < y_se_cutoff,]

write.table(se_population, file='0-se-population.R.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(te_population, file='0-te-population.R.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

# >3000 stretch
stretch_cutoff = 3000
stretch_population = bed[bed$V3-bed$V2 >= stretch_cutoff,]
stretch_se_population = se_population[se_population$V3-se_population$V2 >= stretch_cutoff,]
write.table(stretch_population, file='0-stretch-population.R.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(stretch_se_population, file='0-stretch-se-population.R.bed', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)


## Save a diagnostic plot
png(paste('se-cutoff.R.png', sep=''), width=800, height=800)
plot(x, ynorm, cex=0.5)
abline(h=0, col=8) # baseline zero
lines(spl, col=2) # spline
abline(h=ynorm[[se_cutoff]], col=8)
abline(v=se_cutoff, col=8)

points(pred0, col=2, pch=19) # point to predict tangent 
lines(x, yint + pred1$y*x, col=3) # tangent (1st deriv. of spline at se_cutoff)
points(xint, 0, col=3, pch=19) # x intercept
dev.off()

###
# Export pie charts & histograms
###

## First, make some data frames to use ggplot2
se_sizes = se_population$V3-se_population$V2
se_signals = se_population$V7*se_sizes
te_sizes = te_population$V3-te_population$V2
te_signals = te_population$V7*te_sizes

# size_df = data.frame(c(se_sizes, te_sizes), factor(rep(c('se','te'), c(length(se_sizes), length(te_sizes)))))
# colnames(size_df) <- c('size', 'type')
# signal_df = data.frame(c(se_signals, te_signals), factor(rep(c('se','te'), c(length(se_signals), length(te_signals)))))
# colnames(signal_df) <- c('signal', 'type')


# Pie chart of % SE vs % TE [count]
png(paste('se-vs-te-count-pie.R.png', sep=''), width=800, height=800)
pie(c(nrow(se_population), nrow(te_population)),
    labels=c(paste('Super Enhancers\n', nrow(se_population)),
             paste('Traditional Enhancers\n', nrow(te_population))),
    main="Number of Super- vs Traditional Enhancers")
dev.off()

# Stretch vs Non-Stretch
png(paste('stretch-vs-nonstretch-count-pie.R.png', sep=''), width=800, height=800)
pie(c(nrow(stretch_population), nrow(bed)-nrow(stretch_population)),
    labels=c(paste('Stretch Enhancers\n', nrow(stretch_population)),
             paste('Non-Stretch Enhancers\n', nrow(bed)-nrow(stretch_population))),
    main="Number of Stretch vs Non-Stretch Enhancers")
dev.off()

# And SE stretch / traditional stretch, etc.
png(paste('se-te-stretch-vs-nonstretch-count-pie.R.png', sep=''), width=800, height=800)
pie(c(nrow(se_population)-nrow(stretch_se_population), nrow(stretch_se_population),
      nrow(te_population)-(nrow(stretch_population)-nrow(stretch_se_population)),
      nrow(stretch_population)-nrow(stretch_se_population)),
    labels=c(paste('Super Enhancers\n', nrow(se_population)-nrow(stretch_se_population)),
             paste('Super-Stretch Enhancers\n', nrow(stretch_se_population)),
             paste('Traditional Enhancers\n', nrow(te_population)-(nrow(stretch_population)-nrow(stretch_se_population))),
             paste('Traditional-Stretch Enhancers\n', nrow(stretch_population)-nrow(stretch_se_population))),
    main="Number of Super- vs Traditional Enhancers w/ Stretch >3kb")
dev.off()


# Pie chart of % SE Signal vs % TE Signal
total_signal = sum(as.numeric(bed$V5))
se_signal = sum(as.numeric(se_population$V5))
# te_signal = sum(as.numeric(te_population$V5))
se_fraction = round(se_signal/total_signal*100, digits = 1)
te_fraction = 100-se_fraction

png(paste('se-vs-te-signal-pie.R.png', sep=''), width=800, height=800)
pie(c(se_fraction, te_fraction),
    labels=c(paste('Super Enhancers\n', se_fraction, '%', sep = ''),
             paste('Traditional Enhancers\n', te_fraction, '%', sep = '')),
    main="Signal in Super- vs Traditional Enhancers")
dev.off()


###### Histograms
# Histogram of SE sizes
png(paste('se-size-histogram.R.png', sep=''), width=800, height=800)
hist(se_signals, breaks=10,
     xlab = 'Super-Enhancer Sizes', ylab="Counts", main="Super-Enhancer Size Distribution")
dev.off()

png(paste('te-size-histogram.R.png', sep=''), width=800, height=800)
hist(te_signals, breaks=10,
     xlab = 'Super-Enhancer Sizes', ylab="Counts", main="Traditional-Enhancer Size Distribution")
dev.off()

# Inverse hockeystick, honestly
# qplot(signal, data=signal_df[signal_df$type == 'se',], geom="histogram",
#      xlab="Normalized ATAC Signal", ylab="Count")

######## Text Diagnostics
# Print some raw statistics to a text file

fh <- file("0-enhancer-stats.txt", "w")
writeLines(paste("Statistics for:", input_filename), con=fh)
writeLines(paste("SE Signal %:", round(se_fraction,2)), con=fh)
writeLines(paste("TE Signal %:", round(te_fraction,2)), con=fh)
writeLines(paste("SE Count:", nrow(se_population)), con=fh)
writeLines(paste("TE Count:", nrow(te_population)), con=fh)
writeLines(paste("SE Count %:", round(nrow(se_population)/nrow(bed)*100, 2)), con=fh)
writeLines(paste("TE Count %:", round(nrow(te_population)/nrow(bed)*100, 2)), con=fh)
writeLines(paste("Mean SE Size:", round(mean(se_sizes), 2)), con=fh)
writeLines(paste("Mean TE Size:", round(mean(te_sizes), 2)), con=fh)
writeLines(paste("Median SE Size:", median(se_sizes)), con=fh)
writeLines(paste("Median TE Size:", median(te_sizes)), con=fh)
close(fh)

