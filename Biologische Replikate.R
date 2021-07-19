Genome1 <- read.csv("Coverage 03-05 NEW SeqA5 E.coli low S.csv")
SynVic1 <- read.csv("Coverage 06-05 NEW SeqA5 synVic low S.csv")
SynVic2 <- read.csv("Coverage 06-05 NEW SeqA8 synVic low S.csv")
Genome2 <- read.csv("Coverage 06-05 NEW SeqA8 Ecoli low S.csv")
names(Genome1) <- c("Position","Genome1")
names(Genome2) <- c("Position","Genome2")
names(SynVic1) <- c("Position","SynVic1")
names(SynVic2) <- c("Position","SynVic2")

Plasmid <- read.fasta("E.coli Genome.fasta", as.string = TRUE, seqtype = c("DNA", "AA"), seqonly = TRUE,)
Positions <- str_locate_all(Plasmid, "GATC")
PosFrame <- as.data.frame.list(Positions)

PosFrame <- select(PosFrame,start)
PosFrame <- left_join(PosFrame,Genome1,by = c("start" = "Position"))
PosFrame <- left_join(PosFrame,Genome2,by = c("start" = "Position"))


summary(lm(Genome1~Genome2, data = PosFrame))
plot(Genome1~Genome2, data = PosFrame)
abline(lm(Genome1~Genome2, data = PosFrame))

