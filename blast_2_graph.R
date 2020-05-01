############# This is blast format 6
# 1. 	 qseqid 	 query (e.g., unknown gene) sequence id
# 2. 	 sseqid 	 subject (e.g., reference genome) sequence id
# 3. 	 pident 	 percentage of identical matches
# 4. 	 length 	 alignment length (sequence overlap)
# 5. 	 mismatch 	 number of mismatches
# 6. 	 gapopen 	 number of gap openings
# 7. 	 qstart 	 start of alignment in query
# 8. 	 qend 	 end of alignment in query
# 9. 	 sstart 	 start of alignment in subject
# 10. 	 send 	 end of alignment in subject
# 11. 	 evalue 	 expect value
# 12. 	 bitscore 	 bit score


##################### This code takes the blast output and generates a graph

library(stringr)
library(sqldf)

dat <- read.csv("homology/mouse_all2lall_blast_prot.csv",sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(dat) <- c(
  "from_gene","to_gene","pident","length","mismatch","gapopen",
  "qstart","qend","sstart","send","evalue","bitscore")

dat <- dat[,c("from_gene","to_gene","pident","length","evalue","bitscore")]

head(dat)


#### Remove self-comparisons
dat <- dat[dat$from_gene!=dat$to_gene,]
# Might need to force blast to return more values!


#### Make symmetric - no idea what it is currently otherwise
dat <- rbind(
  data.frame(
    from_gene=dat$from_gene,
    to_gene=dat$to_gene,
    pident=dat$pident,
    length=dat$length,
    evalue=dat$evalue,
    bitscore=dat$bitscore
  ),
  data.frame(
    from_gene=dat$to_gene,
    to_gene=dat$from_gene,
    pident=dat$pident,
    length=dat$length,
    evalue=dat$evalue,
    bitscore=dat$bitscore
  )
)


####  Keep the best comparison for each to-from
dat2 <- sqldf("select distinct from_gene, to_gene, max(pident) as pident from dat group by from_gene, to_gene")


#### For each "from", only retain the top 5 neigbours


#### TODO - maybe - cut-off for E-value. But keep at minimum one neighbour




hist(dat2$pident, breaks = 100)
dat3 <- dat2[dat2$pident>30,]
hist(dat3$pident, breaks = 100)

#Save the graph
dim(dat)
dim(dat2)
write.table(dat3, "homology/homology_graph.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
