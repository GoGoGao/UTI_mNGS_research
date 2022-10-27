mkdir pcoa.output
perl /data/Scripts/Jiangwei/PCoA/lib/Beta_diversity_index.pl --index BCD --rank --matrix pcoa.output/BCD.mat $1
perl -ne '{chomp;if ($_=~/#/){my @or=split /\t/;shift @or;print "\t",join("\t",@or),"\n";}else{print $_,"\n";}}' pcoa.output/BCD.mat > pcoa.output/all.even.mat.xls
perl /data/Scripts/Jiangwei/PCoA/lib/PCoAclust.pl pcoa.output/all.even.mat.xls $2 pcoa.output/
