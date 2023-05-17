library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

###################################
# proteinGroups.txt file curation #
###################################
proteins <- read.table("proteinGroups.txt", header = TRUE, sep = "\t")

proteins <- proteins[proteins$Potential.contaminant != "+", ]
proteins <- proteins[proteins$Reverse != "+", ]
proteins <- proteins[proteins$Only.identified.by.site != "+", ]

cnames <- c("^Majority.protein.IDs$", "^Protein.names$", "^Gene.names$", "^Sequence.length$", "^MS.MS.count.")

proteins.1 <- proteins[, unlist(lapply(cnames, function(x) grep(x, colnames(proteins))))]

#In case your MAXQUANT Names and the experiment IDs do not match, change column names in proteins.1 dataset----
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_6d_r1"] <- "MS.MS.count.CTRL_6d_r1"
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_6d_r2"] <- "MS.MS.count.CTRL_6d_r2"
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_6d_r3"] <- "MS.MS.count.CTRL_6d_r3"
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_r1"] <- "MS.MS.count.CTRL_6d_r4"
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_r2"] <- "MS.MS.count.CTRL_6d_r5"
names(proteins.1)[names(proteins.1) == "MS.MS.count.Ctrl_r3"] <- "MS.MS.count.CTRL_6d_r6"

names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_6d_r1"] <- "MS.MS.count.IP_6d_r1"
names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_6d_r2"] <- "MS.MS.count.IP_6d_r2"
names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_6d_r3"] <- "MS.MS.count.IP_6d_r3"
names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_r1"] <- "MS.MS.count.IP_6d_r4"
names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_r2"] <- "MS.MS.count.IP_6d_r5"
names(proteins.1)[names(proteins.1) == "MS.MS.count.TET2_r3"] <- "MS.MS.count.IP_6d_r6"

names(proteins.1)
#----

#######################################
## Creating files for SAINT analysis ##
#######################################

## NOTE: these files must be MANUALLY created for every comparison you want to perform...

## MOFWT vs CTRL

# bait file example

#    IP.name Bait T_C
#    MOFWT1  MOF   T
#    MOFWT2  MOF   T
#    MOFWT3  MOF   T
#    CTRL1  IgG   C
#    CTRL2  IgG   C
#    CTRL3  IgG   C

#AUTOMATITZAR LA CREACIÓ DEL BAIT.FILE !!!!!!

bait.file <- data.frame("IP name" = paste(rep(c("IP_6d_r", "CTRL_6d_r"), each = 6), 1:6, sep = ""), "Bait" = c(rep(c("TET2", "IgG"), each = 6)), "T_C" = rep(c("T", "C"), each = 6))

# prey file example

#       Protein.Name   Seq..Length    Gene.Name
#       Q9Y6J9          622           TAF6L
#       Q9Y6K1          912           DNMT3A
#       Q9Y6M1          599           IGF2BP2
#       Q9Y6V7          483           DDX49
#       Q9Y6X4          670           FAM169A
#       Q9Y6X9;Q86VD1  1032           MORC2;MORC1

# NOTE: to create the 'prey file' we need the object called 'proteins.1' 
prey.file <- data.frame("Protein Name" = proteins.1$Majority.protein.IDs, "Seq. Length" = proteins.1$Sequence.length, "Gene Name" = proteins.1$Gene.names)

# interaction file example
# NOTE: THIS FILE IS GENERATED AUTOMATICALLY ONCE THE 'bait file' and 'prey file' ARE MANUALLY CREATED:

#   MOFWT1	MOF	A0A075B6Q5	0
#   MOFWT1	MOF	A0A0A0MS14	3
#   MOFWT1	MOF	A0A0B4J1X5	0
#   MOFWT1	MOF	A5YKK6	0
#   MOFWT1	MOF	A6NHQ2	5
#   MOFWT1	MOF	A6NHR9	14
#   MOFWT1	MOF	A6NIH7	2
#   MOFWT1	MOF	P62834;P61224;A6NIZ1	1
#   MOFWT1	MOF	A6NJ78;P0C7V9	1
#   MOFWT1	MOF	A8MW92	0
#   MOFWT1	MOF	P62308;A8MWD9	1

if.ipname <- as.character(rep(bait.file$IP.name, each = nrow(prey.file)))
if.baitname <- as.character(rep(bait.file$Bait, each = nrow(prey.file)))
if.preyname <- as.character(rep(prey.file$Protein.Name, nrow(bait.file)))
interaction.file <- cbind(if.ipname, if.baitname, if.preyname)

# the last column of the 'interaction file' needs the object created with the 'PSM_Calculation.R' script which contains all the PSMs for all the samples (proteins.1)

if.spc <- c()
for (i in 1:nrow(bait.file)) {
  p <- proteins[, grep(paste("MS.MS.count.", as.character(bait.file$IP.name[i]), sep = ""), colnames(proteins.1))]
  if.spc <- c(if.spc, p) 
}

interaction.file <- cbind.data.frame(interaction.file, if.spc)

write.table(bait.file, "bait.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(prey.file, "prey.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(interaction.file, "interaction.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


