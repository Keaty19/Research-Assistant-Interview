#!/bin/bash
#SBATCH -p htc
#SBATCH --mem=40G
#SBATCH --ntasks=2
#SBATCH --tasks-per-node=4
#SBATCH -t 24:00:00
#SBATCH -o DpTlsP.%J
#SBATCH -e DpTlsP.%J
#SBATCH --job-name=DpTls2
#SBATCH --account=scw1448

module load deeptools

computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName wt1coverage.bigWig \
--outFileName wt1Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep
computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName wt2coverage.bigWig \
--outFileName wt2Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep
computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName wt3coverage.bigWig \
--outFileName wt3Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep
computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName ko4coverage.bigWig \
--outFileName ko4Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep
computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName ko5coverage.bigWig \
--outFileName ko5Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep
computeMatrix reference-point --regionsFileName genesPval.pad1000.bed --scoreFileName ko6coverage.bigWig \
--outFileName ko6Matrix.tab.gz --referencePoint center --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions keep

plotHeatmap --matrixFile wt1Matrix.tab.gz --outFileName wt1pValHeatMap.png --sortRegions no
plotHeatmap --matrixFile wt2Matrix.tab.gz --outFileName wt2pValHeatMap.png --sortRegions no
plotHeatmap --matrixFile wt3Matrix.tab.gz --outFileName wt3pValHeatMap.png --sortRegions no
plotHeatmap --matrixFile ko4Matrix.tab.gz --outFileName ko4pValHeatMap.png --sortRegions no
plotHeatmap --matrixFile ko5Matrix.tab.gz --outFileName ko5pValHeatMap.png --sortRegions no
plotHeatmap --matrixFile ko6Matrix.tab.gz --outFileName ko6pValHeatMap.png --sortRegions no