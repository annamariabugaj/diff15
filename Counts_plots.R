library("DESeq2")
library("ggplot2")

#df1 <-read.csv('subtable_significant_genes_WT_VSD.csv', header = TRUE, sep = ",")
head(df1)

meta1 <- read.csv('meta_allis.csv', header = TRUE, sep = ',')
head(meta1)

countData <- df1
head(countData)
metaData <- meta1
head(metaData)

design<- (~ genotype+age+region)

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = design, tidy = TRUE)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)


############################# PLOT COUNTS ################################

# geneCounts <- plotCounts(dds, gene="Egr3", intgroup=c("age", "genotype", "region"), returnData=TRUE)
# Plot the data using ggplot2
# colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
# ggplot(geneCounts, aes(x=interaction(genotype, age), y=count, colour=genotype, group=genotype)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Tdg")
# ggplot(geneCounts, aes(x=interaction(genotype, age, region), y=count, colour=genotype, group=genotype)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Tdg")
# ggplot(geneCounts, aes(x=interaction(age, region, genotype), y=count, colour=genotype)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Tdg")

Counts <- plotCounts(dds, gene="Tmem159", intgroup=c("genotype", "age", "region"), returnData=TRUE)
# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
p <- ggplot(Counts, aes(x=age, y=count, colour=genotype, shape = genotype, group=genotype)) + geom_point()  + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Tmem159 pink module")
p <- p  + ylab("expression counts") 
p
p + facet_grid(.~region)

Counts <- plotCounts(dds, gene="Pde11a", intgroup=c("genotype", "age", "region"), returnData=TRUE)
# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
p <- ggplot(Counts, aes(x=age, y=count, colour=genotype, shape = genotype, group=genotype)) + geom_point()  + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Pde11a pink module")
p <- p  + ylab("expression counts") 
p
p + facet_grid(.~region)

Counts <- plotCounts(dds, gene="Corob1", intgroup=c("genotype", "age", "region"), returnData=TRUE)
# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
p <- ggplot(Counts, aes(x=age, y=count, colour=genotype, shape = genotype, group=genotype)) + geom_point()  + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Corob1 black module")
p <- p  + ylab("expression counts") 
p
p + facet_grid(.~region)

Counts <- plotCounts(dds, gene="Coro1b", intgroup=c("genotype", "age", "region"), returnData=TRUE)
# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
p <- ggplot(Counts, aes(x=age, y=count, colour=genotype, shape = genotype, group=genotype)) + geom_point()  + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Coro1b magenta module")
p <- p  + ylab("expression counts") 
p
p + facet_grid(.~region)

Counts <- plotCounts(dds, gene="Sipa1l1", intgroup=c("genotype", "age", "region"), returnData=TRUE)
# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e", "#463d25")
p <- ggplot(Counts, aes(x=age, y=count, colour=genotype, shape = genotype, group=genotype)) + geom_point()  + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Sipa1l1 magenta module")
p <- p  + ylab("expression counts") 
p
p + facet_grid(.~region)
