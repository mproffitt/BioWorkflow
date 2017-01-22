
working.dir <- "/home/mproffitt/Sequences/run"
setwd(working.dir)

dest.dir <- paste0(working.dir,"/bam_merged")
# dir.create(dest.dir)

input.dir <- paste0(working.dir,"/bam_lanewise")
bamfiles <- list.files(input.dir, pattern=".bam$")

commands.dir <- paste0(working.dir,"/log")
command.file <- file.path(commands.dir, "merge_bam_commands.txt")

output.dir <- paste0(working.dir,"/alignment_stats")

pnames <- gsub("_L00().bam", "", bamfiles)
treatment <- vector("character", length(pnames))


treatment[grep("-Input_", pnames)] <- "GL30_Input"
treatment[grep("Hd2lox_Hd1_", pnames)] <- "GL30_Hd2lox_Hd1"
treatment[grep("Hd2lox-Hd2_", pnames)] <- "GL30_Hd2lox_Hd2"
treatment[grep("Hd2delta_Hd1_", pnames)] <- "GL30_Hd2delta_Hd1"
treatment[grep("Hd2delta_Hd2_", pnames)] <- "GL30_Hd2delta_Hd2"
treatment[grep("_RT7_S8_", pnames)] <- "GL30_CoREST"

groups <- factor(treatment, levels=unique(treatment))

sample.desc <- data.frame(filename=bamfiles, key=groups)
fac <- unique(as.character(sample.desc$key))


for(i in 1:length(fac)){
    ID = fac[i]
    cat("Merging ", ID, "\n", sep="")
    sel <- which(as.character(sample.desc[, "key"]) == ID)
    in.bam = paste(file.path(input.dir, as.character(sample.desc[sel, "filename"])), collapse=" ")
    outpath.sampleName = paste0(file.path(dest.dir, ID), ".bam")

    mergeBam.cmd <- paste("samtools merge - ", in.bam, " | samtools sort -o ", outpath.sampleName, " - ", sep=" ")
    print(mergeBam.cmd)
    indexBam.cmd <- paste("samtools index ", outpath.sampleName, sep="")

    system(mergeBam.cmd)
    system(indexBam.cmd)

    write(ID, file=command.file, sep="\n", append=TRUE)
    write(mergeBam.cmd, file=command.file, sep="\n", append=TRUE)
    write(indexBam.cmd, file=command.file, sep="\n", append=TRUE)
    write("\n", file=command.file, sep="\n", append=TRUE)
}

