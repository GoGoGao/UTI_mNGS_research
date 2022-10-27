library(vegan)

##读入文件
#OTU 丰度表
otu <- read.delim('kraken2.abs.g.sort.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#样本分组文件
group <- read.delim('man.group.list', sep = '\t', stringsAsFactors = FALSE)

##PERMANOVA 分析（所有分组间比较，即整体差异）
#根据 group$site 这一列分组进行 PERMANOVA 分析检验组间差异，基于 999 次置换，详情 ?adonis

#（1）直接输入 OTU 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
adonis_result <- adonis(otu~site, group, distance = 'bray', permutations = 999)

otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
