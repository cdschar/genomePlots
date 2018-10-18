#specify the input refFlat table from UCSC browser download pages
file = "refFlat.txt.gz"

#specify the output bed file name by changing the genome
out_file = "hg38.refFlat.bed"

if(grepl("gz", file)){
	gz = gzfile(file, 'rt')
	dat = read.csv(gz,header=F,sep="\t")
} else {
	dat = read.csv(file,header=F,sep="\t")
}

col1 = c()
col2 = c()
for(i in 1:dim(dat)[1]){
	x1 = unlist(strsplit(as.vector(dat[i,11]),","))
	x2 = unlist(strsplit(as.vector(dat[i,10]),","))
	len1 = paste(as.numeric(x1)-as.numeric(x2), collapse=",")
	len2 = paste(as.numeric(x1)-dat[i,5], collapse=",")
	col1 = c(col1, len1)
	col2 = c(col2, len2)
}

bed_out = data.frame( dat[, c(3,5,6,1)],0,dat[,c(4,7,8)],0,dat[,9], col1, col2)
write.table(bed_out, out_file, quote=F, sep="\t", col.names=F, row.names=F)
