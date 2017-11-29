library(Rsamtools)
bamcoverage <- function (bamfile) {
  # read in the bam file
  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
  # filter reads without match position
  ind <- ! is.na(bam$pos)
  ## remove non-matches, they are not relevant to us
  bam <- lapply(bam, function(x) x[ind])
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ## names of the bam data frame:
  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
  ## construc: genomic ranges object containing all reads
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
  return (mean(coverage(ranges)))      
}
bams = list.files()
bams = bams[grep('.filter.sort.bam',bams)]
bams = bams[-grep('bai',bams)]

test = bamcoverage(bams[1])
coverage = matrix(0, nrow=length(test), ncol=length(bams))
colnames(coverage) = gsub(".filter.sort.bam","",bams)
rownames(coverage) = names(test)
coverage[names(test),1] = test
for(i in 2:length(bams)){
  print(i)
  test = bamcoverage(bams[i])
  coverage[names(test),i] = test
}

write.csv(coverage,file="coverage773_pediatric.csv")

