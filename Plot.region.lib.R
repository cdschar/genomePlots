#genomePlot v2.1
#10/16/18

library("Biostrings");
library("BSgenome");
library("rtracklayer");

plotRegion = function(track.file, chr, start, end, output = c("screen", "svg", "pdf"), output.file = c("plot.pdf"), plot.fwd = c(TRUE), plot.pattern = c(FALSE), pattern = c(NA), pattern.height = c(0.5),
pattern.col = c(rgb(0, 0.5, 0)), org = c("Mmusculus.UCSC.mm9"), fixed.width = c(TRUE), width = c(6), bp.length = c(0.001), xlab.inc = c(5000), xlab.inc.on = c(FALSE),
xaxis.win.height = c(0.20), xaxis.size = c(2), xlab.size = c(0.2), xlab.line = c(4), ylab.inc = c(0.5), ylab.inc.on = c(FALSE), arrow.lwd = c(1), arrow.height = c(0.8), arrow.length = c(0.01), arrowhead.size = c(0.1),
outer.margin = c(c(0.75, 0.75, 0.5, 0.75)), track.margin = c(c(0.05, 0, 0.05, 0)), legend.loc = c("topleft"), legend.horiz = c(FALSE), legend.size = c(1), mgp = c(3, 1, 0), ptsPch = c(19), ptsCex = c(0.5), axis.lwd = c(0.5)) {

	#Debug
	#track.file = track.file; org = "Mmusculus.UCSC.mm9"; chr = as.character(regions.plot$chr[region]); start = as.numeric(regions.plot$start[region]); end = as.numeric(regions.plot$end[region]); output = "screen"; output.file = output.file; gene.win.height = 1; fixed.width = TRUE; width = 8; plot.fwd = TRUE; xaxis.size = 2; xlab.size = 0.8; xlab.inc = 30000; xaxis.win.height = 0.20; arrow.lwd = 3; xaxis.size = 2; ylab.inc = 0.5; ylab.inc.on = FALSE; arrowhead.size = 0.05; arrow.length = 0.01; arrow.height = 0.6; outer.margin = c(0.75, 0.75, 0.5, 0.75); track.margin = c(0.05, 0, 0.05, 0); xlab.inc.on = FALSE; legend.loc = c("topleft"); legend.horiz = c(FALSE); xlab.line = 4; legend.size = 1.5;

	library(paste("BSgenome", org, sep = "."), character.only = TRUE);	#load appropriate library

	org.common = gsub("\\.[[:alnum:]]+", "", org);	#Parse organism common name
	seq = subseq(eval(parse(text = org.common))[[which(seqnames(eval(parse(text = org.common))) == chr)]], start = start, end = end);
	seq.length = nchar(seq);

	#Find Pattern
	if (plot.pattern) pat = matchPattern(pattern, seq);

	#Read in track data
	tracks <- read.table(track.file, sep = "\t", header = TRUE, as.is = TRUE);
	tracks.inc = tracks[tracks$show, ];
	tracks.inc = tracks.inc[order(tracks.inc$order, decreasing = FALSE), ];

	track.heights = NULL;

	for (track in unique(tracks.inc$track)) track.heights = c(track.heights, max(tracks.inc$height[tracks.inc$track == track]) + track.margin[1] + track.margin[3]);

	track.heights[1] = track.heights[1] + outer.margin[3];
	track.heights = c(track.heights, xaxis.win.height + outer.margin[1]);

	if (fixed.width) pic.width = width else pic.width = ((end - start) * bp.length) + outer.margin[2] + outer.margin[4];

	if (output == "svg") {
		svg(file = output.file, width = pic.width, height = sum(track.heights));
		par(family = "Arial");
	} else if (output == "pdf") {
		cairo_pdf(file = output.file, width = pic.width, height = sum(track.heights));
		par(family = "Arial");
	} else if (output == "screen") {
		X11(width = pic.width, height = sum(track.heights));
		par(family = "Arial");
	}

	plot.format <- layout(matrix(c(1:length(track.heights)), length(track.heights), 1, byrow=TRUE), heights = track.heights, widths = pic.width)
	#layout.show(plot.format);

	#plot tracks
	for (track in unique(tracks.inc$track)) {

		this.track = tracks.inc[tracks.inc$track == track, ];
		this.track = this.track[order(this.track$order, decreasing = TRUE), ];

		if (min(this.track$order) == min(tracks.inc$order)) {
			par(mai = c(track.margin[1], outer.margin[2], outer.margin[3] + track.margin[3], outer.margin[4]), bty = "n");
		} else {
			par(mai = c(track.margin[1], outer.margin[2], track.margin[3], outer.margin[4]), bty = "n");
		}

		if (is.na(max(this.track$max))) {
			ymax.track = readTrack(this.track$type[1], this.track$file[1], chr, start, end);
			ymax = max(bedY(ranges(ymax.track)), na.rm = TRUE) + 1;	#Determine y max
		} else {
			ymax = max(this.track$max);
		}

		if (plot.fwd) {
			plot(NA, ylim = c(0, ymax), xlim = c(start, end), xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.lab = 1.5);
		} else {
			plot(NA, ylim = c(0, ymax), xlim = rev(c(start, end)), xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.lab = 1.5);
		}

		legend = NULL;
		legend.col = NULL;
		legend.lwd = NULL;
		legend.lty = NULL;

		for (i in 1:dim(this.track)[1]) {
	
			print(paste("Reading in track: ", this.track$name[i], " - ", Sys.time())); 
			track.data = readTrack(this.track$type[i], this.track$file[i], chr, start, end);	

			if (this.track$legend[i]) {
				legend = c(legend, this.track$name[i]);
				legend.lwd = c(legend.lwd, this.track$line.lwd[i]);
				legend.lty = c(legend.lty, this.track$line.lty[i]);
			
				if (!is.na(eval(parse(text = this.track$color[i])))) {
					legend.col = c(legend.col, eval(parse(text = this.track$color[i])));
				} else if (!is.na(eval(parse(text = this.track$border[i])))) {
					legend.col = c(legend.col, eval(parse(text = this.track$border[i])));
				} else {
					legend.col = c(legend.col, rgb(0, 0, 0));
				}
			}

			if (this.track$type[i] == "bed") {
				if (this.track$contains.color[i] & !is.null(track.data$itemRgb)) {
					track.cols = track.data$itemRgb;
					track.borders = track.data$itemRgb;
				} else {
					track.cols = rep(eval(parse(text = this.track$color[i])), length(ranges(track.data)));
					track.borders = rep(eval(parse(text = this.track$border[i])), length(ranges(track.data)));
				}

				if (length(ranges(track.data)) > 0) {	#if there are bed tracks				
					y.mean = bedY(ranges(track.data));	#Determine y location
					for (j in 1:length(ranges(track.data))) {		#for each bed track

						rect(start(track.data)[j], y.mean[j] - (this.track$bed.height[i]) / 2, end(track.data)[j], y.mean[j] + (this.track$bed.height[i]) / 2, col = track.cols[j], border = track.borders[j]);

						#Draw block regions (e.g. exons)
						if (!is.null(track.data$blocks) & !is.na(this.track$block.height)) {
							for (block in 1:length(start(track.data$blocks[[j]]))) {
								if(end(track.data$blocks[[j]])[block] > start & start(track.data$blocks[[j]])[block] < end) {
									rect(start(track.data$blocks[[j]])[block], y.mean[j] - (this.track$block.height[i] / 2), end(track.data$blocks[[j]])[block], y.mean[j] + (this.track$block.height[i] / 2), col = track.cols[j], border = track.borders[j]);
								}
							}
						}
						#Label bed elements
						if (this.track$bed.name && as.character(strand(track.data))[j] == "*") {
							text(start(track.data)[j], y.mean[j], track.data$name[j], pos = 2, cex = this.track$label.size[i], font = 0.1); 
						}
						#Draw arrows & label tracks
						if (!is.null(as.character(strand(track.data))[j]) && as.character(strand(track.data))[j] == "+") {
							if (!is.na(this.track$bed.arrow) && this.track$bed.arrow) {
								lines(c(start(track.data)[j], start(track.data)[j]), c(y.mean[j], y.mean[j] + arrow.height), lwd = arrow.lwd, col = track.cols[j]);
								arrows(start(track.data)[j], y.mean[j] + arrow.height, start(track.data)[j] + ((end - start) * arrow.length), length = arrowhead.size, lwd = arrow.lwd, col = track.cols[j])
							}
							if (!is.na(this.track$bed.name) & this.track$bed.name) {
								if (start(track.data)[j] > start) {
									if (plot.fwd) pos = 2 else pos = 4;
									text(start(track.data)[j], y.mean[j], track.data$name[j], pos = pos, cex = this.track$label.size[i], font = 3); 
								} else if (end(track.data)[j] < end) {
									if (plot.fwd) pos = 4 else pos = 2;
									text(end(track.data)[j], y.mean[j], track.data$name[j], pos = pos, cex = this.track$label.size[i], font = 3); 
								}
							}
						} else if (!is.null(as.character(strand(track.data))) && as.character(strand(track.data))[j] == "-") {
							if (!is.na(this.track$bed.arrow) && this.track$bed.arrow) {
								lines(c(end(track.data)[j], end(track.data)[j]), c(y.mean[j], y.mean[j] + arrow.height), lwd = arrow.lwd, col = track.cols[j]);
								arrows(end(track.data)[j], y.mean[j] + arrow.height, end(track.data)[j] - ((end - start) * arrow.length), length = arrowhead.size, lwd = arrow.lwd, col = track.cols[j]);
							}
							if (!is.na(this.track$bed.name) & this.track$bed.name) {
								if (end(track.data)[j] > start) {
									if (plot.fwd) pos = 4 else pos = 2;
									text(end(track.data)[j], y.mean[j], track.data$name[j], pos = pos, cex = this.track$label.size[i], font = 3); 
								} else if (start(track.data)[j] < end) {
									if (plot.fwd) pos = 2 else pos = 4;
									text(start(track.data)[j], y.mean[j], track.data$name[j], pos = pos, cex = this.track$label.size[i], font = 3); 
								}
							}
						}
					}
				}
			} else if (this.track$type[i] == "bigWig" | this.track$type[i] == "bedGraph") {

				#Plot Tag data
				tag.col = eval(parse(text = this.track$color[i]));
				border.col = eval(parse(text = this.track$border[i]));
				track.lwd = this.track$line.lwd[i];
				track.lty = this.track$line.lty[i];

				if (this.track$plot.type[i] == "polygon" || is.na(this.track$plot.type[i])) {
					#Fill in zero data on track
					zeroTrack = zeroTrackData(start(track.data), end(track.data), score(track.data));

					if (!is.na(this.track$smooth[i]) & this.track$smooth[i]) {
						#sLine = smooth.spline(x = zeroTrack$locs, y = zeroTrack$scores, spar = 0.1);
						sLine = smooth.spline(x = zeroTrack$locs, y = zeroTrack$scores, spar = 0.2);
						sLine$y[sLine$y < 0] = 0;
						polygon(c(min(sLine$x), sLine$x, max(sLine$x)), c(0, sLine$y, 0), col = tag.col, border = border.col, lwd = track.lwd, lty = track.lty);

					} else {
						polygon(c(min(zeroTrack$locs), zeroTrack$locs, max(zeroTrack$locs)), c(0, zeroTrack$scores, 0), col = tag.col, border = border.col, lwd = track.lwd, lty = track.lty);
					}
					#plot polygon
					
				} else if (this.track$plot.type[i] == "line") {

					zeroTrack = zeroTrackData(start(track.data), end(track.data), score(track.data));				
						
					if (!is.na(this.track$smooth[i]) & this.track$smooth[i]) {
						sLine = smooth.spline(x = zeroTrack$locs, y = zeroTrack$scores, spar = 0.1);
						sLine$y[sLine$y < 0] = 0;
						lines(sLine$x, sLine$y, col = border.col, lwd = track.lwd, lty = track.lty);
					} else {
						lines(zeroTrack$locs, zeroTrack$scores, col = border.col, lwd = track.lwd, lty = track.lty);
					}
				
					#plot 0 line
					lines(c(start, end), c(0, 0));
				} else if (this.track$plot.type[i] == "points") {

					if (!is.na(this.track$ptsPch[i])) ptsPch = this.track$ptsPch[i]
					if (!is.na(this.track$ptsCex[i])) ptsCex = this.track$ptsCex[i]
	
					points(start(track.data), score(track.data), col = tag.col, lwd = track.lwd, pch = ptsPch, cex = ptsCex);
					
				} else if (this.track$plot.type[i] == "bar") {

					#plot 0 line
					lines(c(start, end), c(0, 0));
		
					for (j in 1:length(start(track.data))) rect(start(track.data)[j], 0, end(track.data)[j], score(track.data)[j], col = tag.col, border = border.col, lwd = track.lwd, lty = track.lty);
				}
			}
		}

		if (this.track$type[1] == "bigWig" | this.track$type[1] == "bedGraph") {
			if (ylab.inc.on) {
				axis(side = 2, at = seq(0, max(this.track$max), by = ylab.inc * max(this.track$max)), labels = seq(0, max(this.track$max), by = ylab.inc * max(this.track$max)), cex.axis = max(this.track$label.size, na.rm = TRUE), font = 1, mgp = mgp, lwd = axis.lwd);
			} else {
				axis(side = 2, at = c(0, max(this.track$max)), labels = c(0, max(this.track$max)), cex.axis = max(this.track$label.size, na.rm = TRUE), font = 1, mgp = mgp, lwd = axis.lwd);
			}

			if(any(this.track$legend)) {
				legend.unique = unique(cbind(legend, legend.col, legend.lwd, legend.lty)); #make legend unique
				legend(legend.loc, legend.unique[, "legend"], lwd = as.numeric(legend.unique[, "legend.lwd"]), lty = as.numeric(legend.unique[, "legend.lty"]), text.col = legend.unique[, "legend.col"], col = legend.unique[, "legend.col"], bty = "n", cex = legend.size, horiz = legend.horiz, text.font = 1, x.intersp = 0.25);
			}
		}

		if (!is.na(this.track$label.orientation[1])) lab.las = this.track$label.orientation[1] else lab.las = 0;
		if (!is.na(this.track$label.line[1])) lab.line = this.track$label.line[1] else lab.line = 0;
		if (!is.na(this.track$label.adj[1])) lab.adj = this.track$label.adj[1] else lab.adj = 0;
		if (!is.na(this.track$label.track[1]) && this.track$label.track) { 
			mtext(this.track$label.name[1], side = 2, line = lab.line, cex = this.track$label.size[1], las = lab.las, adj = lab.adj, font = 1);
		}
	}

	par(mai = c(outer.margin[1], outer.margin[2], 0, outer.margin[4]), bty = "n");
	if (plot.fwd) {
		plot(c(start, end), xlim = c(start, end), ylim = c(0, xaxis.win.height), col = rgb(1, 1, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "");
	} else {
		plot(c(start, end), xlim = rev(c(start, end)), ylim = c(0, xaxis.win.height), col = rgb(1, 1, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "");
	}		

	#draw bottom x-axis and annotate with genome position information
	if (xlab.inc.on) {
		axis(side = 1, at = c(start, (start:end)[(start:end - start) %% xlab.inc == 0], end), mgp = mgp, lwd = axis.lwd, labels = FALSE);
	} else {
		axis(side = 1, at = c(start, end), mgp = mgp, lwd = axis.lwd, labels = FALSE);
	}
	mtext(paste0(chr, ":", start, "-", end), side = 1, cex.axis = xlab.size, line = xlab.line, font = 1);

	#replaced above to only draw line and change label format
	#if (xlab.inc.on) {
	#	axis(side = 1, at = c(start, (start:end)[(start:end - start) %% xlab.inc == 0], end), labels = format(c(start, (start:end)[(start:end - start) %% xlab.inc == 0], end), scientific = FALSE), cex.axis = xaxis.size, mgp = mgp, font = 1);
	#} else {
	#	axis(side = 1, at = c(start, end), labels = format(c(start, end), scientific = FALSE), cex.axis = xaxis.size, mgp = mgp, font = 1);
	#}
	#mtext(gsub("chr", "Chromosome ", chr), side = 1, cex.axis = xlab.size, line = xlab.line, font = 1);

	#plot Pattern
	if (plot.pattern) suppressWarnings(rug(start(pat) + start, ticksize = pattern.height, col = pattern.col));


	if (output != "screen") {
		dev.off();
	}
}


#y.mean = bedY(ranges(track.data));	#Determine y location
	
#Function to calculate the bed region based on overlaps
bedY = function(ranges) {
	#Debug: ranges = ranges(ymax.track)

	order = order(start(ranges)) #order bed tracks by start
	y.mean = rep(NA, length(ranges));

	if (length(ranges) > 0) {
		for (i in 1:length(ranges)) {	#for each bed track
			overlaps = findOverlaps(ranges[order][1:i], ranges[order][i]);
			#y.mean[i] = lowestUnfilledInt(y.mean[overlaps@queryHits]);
			y.mean[i] = lowestUnfilledInt(y.mean[queryHits(overlaps)]);
		}
		return(y.mean[match(1:length(order), order)]);	#Return non-ordered y.mean values
	} else {
		return(1);
	}
}

#Function to return the lowest unfilled integer
lowestUnfilledInt = function(ints) {
	int = NA;
	i = 1;

	while (is.na(int)) {
		if (!(i %in% ints)) int = i else i = i + 1;
	}
	return(int);
}

readTrack = function(type, file, chr, start, end, trackLine = c(FALSE)) {

	#Debug
	#type = "bigWig"; file = "/home/bbarwick/Documents/tools/Genome/mm9/phastCons/PhastCons30.placental.bw"; chr = "chr10"; start = 44175000; end = 44180000

	if (type == "bed") {
		track = import.bed(file, which = GenomicRanges::GRanges(chr, IRanges(start, end)));
	} else if (type == "bed15") {
		track = import.bed15(file, which = GenomicRanges::GRanges(chr, IRanges(start, end)));
	} else if (type == "bigWig") {
		track = import.bw(file, which = GenomicRanges::GRanges(chr, IRanges(start, end)));
	} else if (type == "bedGraph") {
		track = import.bedGraph(file, which = GenomicRanges::GRanges(chr, IRanges(start, end)));

	}
	return(track);
}

zeroTrackData = function(starts, ends, scores) {

	#Debug
	#starts = start(track.data[chr]); ends = end(track.data[chr]); scores = score(track.data[chr]);
	zeroLocs = NULL;
	zeroScores = NULL;

	if (length(starts) > 1) {
		for (i in 1:(length(starts) - 1)) {

			zeroLocs = c(zeroLocs, starts[i], ends[i]);
			zeroScores = c(zeroScores, rep(scores[i], 2));

			if (starts[i + 1] != ends[i] + 1) {
				zeroLocs = c(zeroLocs, ends[i] + 1, starts[i + 1] - 1);
				zeroScores = c(zeroScores, rep(0, 2))
			}
		}

		zeroLocs = c(zeroLocs, starts[i + 1], ends[i + 1]);
		zeroScores = c(zeroScores, rep(scores[i + 1], 2));
	}

	track = list(locs = zeroLocs, scores = zeroScores);

	return(track);
}





