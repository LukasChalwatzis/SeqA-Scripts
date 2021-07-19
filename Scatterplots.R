install.packages("stringr")
install.packages("seqinr")
install.packages("dplyr")
install.packages("slider")
library(slider)
library(dplyr)
library(seqinr)
library(stringr)

#Laden der simulierten Coverage aus dem Arbeitsverzeichnis
load("SimulCoverage.Rda")
names(SimulCoverage) <- c("Position","SimCov")


#Das was vorher Fuzznuc gemacht hat
Plasmid <- read.fasta("E.coli Genome.fasta", as.string = TRUE, seqtype = c("DNA", "AA"), seqonly = TRUE,)
Positions <- str_locate_all(Plasmid, "GATC")
PosFrame <- as.data.frame.list(Positions)

#Zordnung der realen Coverage Daten zu den GATC Positionen
Coverage <- read.csv("CoverGenome.csv", sep = ",", header = TRUE)
PosFrame <- select(PosFrame,start)
PosFrame <- left_join(PosFrame,Coverage,by = c("start" = "Position"))
PosFrame <- left_join(PosFrame,SimulCoverage,by = c("start" = "Position"))

#Normierung der Coverage Daten auf ihr Maximum
PosFrame <- mutate(PosFrame,NormCov = PosFrame[,2]/max(PosFrame[,2]))
PosFrame <- mutate(PosFrame,NormSim = PosFrame[,3]/max(PosFrame[,3]))

#Bestimmtheitsmaß, Scatterplot und Ausgleichsgerade
summary(lm(NormSim~NormCov, data = PosFrame))
plot(NormSim~NormCov, data = PosFrame)
abline(lm(NormSim~NormCov, data = PosFrame))

