install.packages("stringr")
install.packages("seqinr")
install.packages("dplyr")
install.packages("slider")
library(slider)
library(dplyr)
library(seqinr)
library(stringr)

#Das was vorher Fuzznuc gemacht hat
Plasmid <- read.fasta("E.coli Genome.fasta", as.string = TRUE, seqtype = c("DNA", "AA"), seqonly = TRUE,)
Positions <- str_locate_all(Plasmid, "GATC")
PosFrame <- as.data.frame.list(Positions)

#Distanzen der GATC Stellen zu den Nachbarn
PosFrame <- mutate(PosFrame, distance = start - lag(start))
PosFrame <- mutate(PosFrame, distance2 = (start -lead(start))*-1)
PosFrame <- mutate(PosFrame, distance3 = start - lag(start,n=2))
PosFrame <- mutate(PosFrame, distance4 = (start -lead(start,n=2))*-1)

#Anzahl der Nachbarn
PosFrame <- mutate(PosFrame,Counter =1)
Nachbarn <- slide_index_sum(PosFrame[,7],PosFrame[,1], before = 185, after = 185)
PosFrame <- mutate(PosFrame,Nach60=Nachbarn)

#Import der Coverage Daten als Data Frame
Coverage <- read.csv("CoverGenome.csv", sep = ",", header = TRUE)

#Zuordnung der Realen Coverages der GATC Positionen
PosFrame <- select(PosFrame,start,distance,distance2,distance3,distance4,Nach60)
PosFrame <- left_join(PosFrame,Coverage,by = c("start" = "Position"))

#Korrektur der Coverage auf die Anzahl der Nachbarn (nur zu Testzwecken)
PosFrame <- mutate(PosFrame, KorrNach = PosFrame[,7]/PosFrame[,6])

#Loop der Für jedes Spacing von 4-34 alle passenden Gatc Stellen auswertet und den Median der Coverage ausgibt
List = list()
for(i in 4:34)
{
  tab <- filter(PosFrame, distance == i|distance2 == i|distance3 == i|distance4 == i)
  LM <- median(tab[,7])
  List[[length(List)+1]] = LM
}
List <- as.numeric(List)

#Plot der Mediancoverage jedes einzelnen Spacings
plot(4:34,List,"l",xlab = "Spacing(bp)", ylab = "Affinity", main="Lukaswerte (Median)")

#Die neuen Werte werden für die Verwendung in der Simulation gespeichert
BrendlerNeu <- data.frame(4:34)
BrendlerNeu <- mutate(BrendlerNeu, Brendler = List)
save(BrendlerNeu,file="BrendlerNeu.Rda")

