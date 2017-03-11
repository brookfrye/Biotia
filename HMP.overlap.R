meta.hmp <- read.delim("~/Dropbox/Ambulance_Study/overlap_meta_clark/HMP_annotations/overlap_v1_metaphlan_HMP_annotation_v1_R.txt")
colnames(meta.hmp)[1] <- "sample.ID"
meta.olf <- read.csv("~/Desktop/biotia/spp_outlap_files_meta_clark_v1/metaphlan.overlap.files.2")
colnames(meta.olf)[1] <- "sample.ID"
hmp.olf <- merge(meta.hmp, meta.olf, by = "sample.ID", all.x = FALSE)

try <- meta.hmp[,2:19]
Bin <- select(try, ends_with("_Binary"))
Sum <- select(try, ends_with("_Sum"))
rownames(Bin) <- hmp.olf$sample.ID
rownames(Sum) <- hmp.olf$sample.ID

#helinger transformation
Sum.scaled[is.na(Sum.scaled)] <- 0 
Sum.Veg <- vegdist(Sum.scaled) 
Sum.Veg[is.na(Sum.Veg)] <- 0

mds.sum <- dudi.pco(Sum.Veg, scannf=F)
scatter(mds.sum)

ppp + geom_point(data=data.frame(mds.sum$li, sample_data(data.frame(meta.full))), 
                 aes(x=A1, y=A2, col=meta.full$sample.Surface), size = 2, alpha=.6) +
  labs(title="MDS: Surfaces") + scale_fill_hue(c=45, l=80)
