#genome plots script

#set working directory
home.dir = "/path/to/current/directory";
setwd(home.dir);

#load source file
lib.file = "/Volumes/GRAID/seqTools/R.Libraries/genomePlot/Plot.region.lib.R";
source(lib.file);

#set output format and directory
output.dir = home.dir

#options:"pdf", "svg", "screen"
output = "pdf"; 

#set species
#org = "Mmusculus.UCSC.mm9";
#org = "Mmusculus.UCSC.mm10";
#org = "Hsapiens.UCSC.hg19";
org = "Hsapiens.UCSC.hg38";

#Define regions to plot
regions.dir = home.dir;
regions.file = "regions.txt";
regions = read.table(paste(regions.dir, regions.file, sep = ""), header = TRUE, sep = "\t", as.is = TRUE);
#limit file for those to plot
regions.plot = regions[regions$plot, ];

#define tracks to plot
track.dir = home.dir
track.file = "tracks.txt";
track.file = paste(track.dir, track.file, sep = "");

#annotate genomic sequence
plot.pattern = FALSE
pattern = "GAATTC"  #this could be CG or any other restriction site

#plot options
width = 6
arrow.length = 0.05
xlab.size = 1
xaxis.size = 1
xlab.line = 1
xaxis.win.height = 0.05
outer.margin = c(0.4, 0.8, 0.25, 0.25)
track.margin = c(0.02, 0, 0.02, 0)
legend.horiz = TRUE


#loop through all the regions and plot
for (region in 1:dim(regions.plot)[1]) {

	#set output file
	outFile = paste(output.dir, regions.plot$name[region], ".", output, sep = "");

	#plot call
	plotRegion(track.file, chr = as.character(regions.plot$chr[region]), start = as.numeric(regions.plot$start[region]), end = as.numeric(regions.plot$end[region]), 
	output = output, output.file = outFile, org = org, plot.fwd = regions.plot$plot.fwd[region], width = width, plot.pattern = plot.pattern, pattern = pattern,  arrow.length = arrow.length, 
	xlab.size = xlab.size, xaxis.size = xaxis.size, xlab.line = xlab.line, xaxis.win.height = xaxis.win.height, outer.margin = outer.margin, track.margin = track.margin, legend.horiz = legend.horiz);

}
