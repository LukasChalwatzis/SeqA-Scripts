install.packages("stringr")
install.packages("seqinr")
install.packages("dplyr")
install.packages("slider")
library(slider)
library(dplyr)
library(seqinr)
library(stringr)
# Alle Packages die für das Skript gebraucht werden


# synVic2-4GATCs.fa
# E.coli Genome.fasta
# MinisynVicII.fasta
#Alternative Sequenzen die hier reinkopiert werden können
Plasmid <- read.fasta("E.coli Genome.fasta", as.string = TRUE, seqtype = c("DNA", "AA"), seqonly = TRUE,)
Positions <- str_locate_all(Plasmid, "GATC")
PosFrame <- as.data.frame.list(Positions)
# macht das was vorher fuzznuc gemacht hat
load("BrendlerNeu.Rda")

Brendlerx <- BrendlerNeu
#Brendlerx <- read.csv("Mappe1.csv", sep = ";", header = TRUE)
#Alternativ können hier die originalen Brendlerwerte hier (Mappe1) ausgewählt werden

names(Brendlerx) <- c("distance", "affinity")
# importiert Nokis Brendlerwerte als Dataframe aus Excel

PosFrame <- mutate(PosFrame, distance = start - lag(start))
PosFrame <- mutate(PosFrame, distance2 = (start -lead(start))*-1)
PosFrame <- mutate(PosFrame, distance3 = start - lag(start,n=2))
PosFrame <- mutate(PosFrame, distance4 = (start -lead(start,n=2))*-1)
# rechnet die Abstände vor und nach jeder GATC Position aus


Test <- left_join(PosFrame,Brendlerx)
colnames(Test)[7] <- "aff1"
Test <- left_join(Test,Brendlerx, by = c("distance2" = "distance"))
colnames(Test)[8] <- "aff2"
Test <- left_join(Test,Brendlerx, by = c("distance3" = "distance"))
colnames(Test)[9] <- "aff3"
Test <- left_join(Test,Brendlerx, by = c("distance4" = "distance"))
colnames(Test)[10] <- "aff4"
# automatische Zuordnung der Brendlerwerte zu den GATC Stellen


allwert <- select(Test,start,aff1,aff2,aff3,aff4)
allwert[is.na(allwert)] = 0
allwert <- mutate(allwert,affmax = pmax(aff1,aff2,aff3,aff4))

#automatische Maximalwert Bestimmung der Brendlerwertes

Wert60 <- slide_index_sum(allwert[,6],allwert[,1], before = 60, after = 60)
Han60 <- (slide_index_sum(allwert[,6],allwert[,1], before = 60, after = 60)-allwert[,4])*0.5
Han185 <- (slide_index_sum(allwert[,6],allwert[,1], before = 185, after = 185) - Wert60)*0.25

#Berechnung der Han werte

Sumall <- allwert[,6]+Han60+Han185
allwert <- mutate(allwert,Endresult = Sumall)


#Berücksichtigung der Lokalen SeqA Konzentration

ExpWert <- slide_index_sum(allwert[,7],allwert[,1], before = 350, after = 350)
allwert <- mutate(allwert, Exp = ExpWert)


#Ausrechnen der Endsumme

allwert <- select(allwert,start,Exp)

#Dataframe für alle Basen der Sequenz
start <- c(1:nchar(Plasmid))
Fullplot <- data.frame(start)

#Zuordnung der GATC Stellen in die gesamte Sequenz
Fullplot <- left_join(Fullplot,allwert)
Fullplot[is.na(Fullplot)] = 0

#Simulation der 75bp großen Illumina Fragmente
SimChip <- slide_sum(Fullplot[,2], before = 75, after = 75)

#Plot der Simulierten Coverage
plot(1:nchar(Plasmid), SimChip,"l",xlab = "Position(bp)", ylab = "Coverage(AU)", main="E.coli Genom(Sim)")
SimulCoverage <- data.frame(1:nchar(Plasmid), SimChip)
save(SimulCoverage,file="SimulCoverage.Rda")









