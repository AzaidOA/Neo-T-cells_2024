### ------------------------------ Libraries ------------------------------ ###
library(data.table)
library(tidyr)
library(stringr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(colorRamps)
library(scales)

### ------------------------------ Functions ------------------------------ ###
plotTiff <- function(file.name,output.dir,ggplot.obj,blank.element=NULL){
  tmp.file.name <- paste0(output.dir,'/',file.name,'.C.tiff')
  tiff(file=tmp.file.name, height=1000, width=1000)
  print(ggplot.obj)
  dev.off()

  if(!is.null(blank.element)){
    ggplot.obj <- ggplot.obj +
                  blank.element
    tmp.file.name <- paste0(output.dir,'/',file.name,'.B.tiff')
    tiff(file=tmp.file.name, height=1000, width=1000)
    print(ggplot.obj)
    dev.off()
  }
}
### ------------------------------- Arguments ------------------------------ ###
srt.obj.tpos.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/cancer/paper_developments/Tfh-assesment/seurat/outs/07-21-2023/220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN-CellType-TetPos/pcs-20_var-25/seurat_objects/SeuratObjectForPrj220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN-CellType-TetPos_WithArgs_NoPCs_20_WithTCRTags_ChangedOn_2023-07-25.RDS'

srt.obj.tmerge.file <- '/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/cancer/paper_developments/Tfh-assesment/seurat/outs/05-17-2023/220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN/pcs-20_var-25/seurat_objects/SeuratObjectForPrj220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN_WithArgs_NoPCs_20_WithTCRTags_ChangedOn_2023-05-24.RDS'

cells.info <- c('/mnt/bioadhoc-temp/Groups/vd-vijay/moarias/sequencing_data/09-21-2022/aggr_vdj/220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN/220901_SPS_Tet_and_TfH_Kinetics_Ryan_220921_Batches-1_Tissue-LN_1.0_09-26-2023')

cells.clons.info.files <- paste0(cells.info,'/filtered_contig_annotations_aggr.csv')

clonotypes.of.int.file <- 'data/TCRs_of_interest_03-2024.csv'

reports.dir <-  'results'

genes.of.int <- c('Cd3e','Cd4','Icos','Foxp3','Junb','Klrg1','Tcf7','Cxcr4','Ccr4','Ccr6','Gzmb','Lef1','Tigit','Gata3','Cd69','Tnf','Cxcr3','Ccr7','Ifng','Ctla4','Il2','Il2rb','Il10','Il10ra','Pdcd1','Cxcr5','Il21','Bcl6')

tags.of.int <- c(
    'cell.type.tag','clonotype.tag','TRA.nt.chains.tag','TRB.nt.chains.tag',
    'TRA.aa.chains.tag','TRB.aa.chains.tag')

feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')

res.of.int <- 'RNA_snn_res.0.8'

blank.complement.3 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_blank(), axis.line=element_blank())

### --------------------- Loading data --------------------- ###
# ---> Loading data.
cells.clons.info <- fread(cells.clons.info.files,stringsAsFactors=FALSE)
clonotypes.of.int <- fread(clonotypes.of.int.file, na.strings="\"NA\"")
srt.obj.tpos <- readRDS(srt.obj.tpos.file)
srt.obj.tmerge <- readRDS(srt.obj.tmerge.file)

### --------------------- Extracting TCRS --------------------- ###
# Keeping good quality cells.
cells.to.keep <- cells.clons.info$is_cell==T & cells.clons.info$high_confidence==T & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive==T
cells.clons.info <- cells.clons.info[cells.to.keep, ..feats.to.keep]

vdj.gene.info <- cells.clons.info[
  ,
  .(
    v.gene=.SD[,
      .(count=.N),
      by=v_gene
    ][count==max(count), paste0(unique(v_gene), collapse='|')],
    j.gene=.SD[,
      .(count=.N),
      by=j_gene
    ][count==max(count), paste0(unique(j_gene), collapse='|')]
  ),
  by=.(
    raw_clonotype_id,
    chain
  )
]

# For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]


# Spread values according to clonotype ID (i.e., to disregard chain information)
vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
vdj.gene.info <- spread(data=vdj.gene.info[, .(raw_clonotype_id, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
# Separate values according to gene type for each chain.
vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)

### --------------------- Downstream analysis: Merged data set --------------------- ###
# ---> Merging TCR data and GEX data
meta.data <- as.data.table(srt.obj.tmerge@meta.data[,tags.of.int],keep.rownames='barcode')

gex.info <- srt.obj.tmerge@assays$RNA@data[genes.of.int,]

gex.info <- as.data.table(t(as.data.frame(gex.info)),keep.rownames='barcode')

meta.data <- merge(x=meta.data, y=gex.info, by='barcode')
meta.data <- meta.data[!is.na(clonotype.tag)]

clonotypes <- meta.data[,.(clone.size.tag=.N),by=tags.of.int]

for(tmp.gene in genes.of.int){
  tmp.meta.data <- meta.data[,c(tags.of.int,tmp.gene),with=FALSE]
  colnames(tmp.meta.data)[colnames(tmp.meta.data) == tmp.gene] <- 'tmp.tag'
  tmp.clonotypes <- tmp.meta.data[,.(tmp.exp=mean(tmp.tag)),by=tags.of.int]
  clonotypes[tmp.clonotypes,on=tags.of.int,tmp.exp:=tmp.exp]
  colnames(clonotypes)[colnames(clonotypes) == 'tmp.exp'] <- tmp.gene
}

clonotypes <- clonotypes[cell.type.tag == 'TetPos']
colnames(clonotypes)[1] <- 'Tet.tag'

clonotypes <- merge(
  x=clonotypes,
  y=vdj.gene.info,
  by.x='clonotype.tag', by.y='raw_clonotype_id'
)

tmp.output.file <- paste0(reports.dir,
  '/tables/Clonotypes_Tet-Positive.csv')
fwrite(x=clonotypes, file=tmp.output.file,na='NA')
tmp.output.file <- paste0(reports.dir,
  '/tables/Clonotypes_Tet-Positive_Foxp3-Positive.csv')
fwrite(x=clonotypes[Foxp3 > 0], file=tmp.output.file,na='NA')
tmp.output.file <- paste0(reports.dir,
  '/tables/Clonotypes_Tet-Positive_ClonalSize-4.csv')
fwrite(x=clonotypes[clone.size.tag >= 4], file=tmp.output.file,na='NA')

# V and j configuration
meta.data <- as.data.table(srt.obj.tmerge@meta.data[,tags.of.int],keep.rownames='barcode')
meta.data <- meta.data[!is.na(clonotype.tag)]
clonotypes <- meta.data[,.(clone.size.tag=.N),by=tags.of.int]
clonotypes <- merge(
  x=clonotypes,
  y=vdj.gene.info,
  by.x='clonotype.tag', by.y='raw_clonotype_id'
)

# TRA
tra.genes <- clonotypes[!is.na(tra.v) & !is.na(tra.j), .(count=.N), by=.(cell.type.tag,tra.v,tra.j)]

tmp.output.file <- paste0(reports.dir,
  '/figures/VDJUsagePieChart_Chain-TRA_Tet-Negative.pdf')
pdf(tmp.output.file)
tmp.p <- pie(tra.genes[order(count)][cell.type.tag == 'TetNeg', count],labels=NA)
print(tmp.p)
dev.off()

tmp.output.file <- paste0(reports.dir,
  '/figures/VDJUsagePieChart_Chain-TRA_Tet-Positive.pdf')
pdf(tmp.output.file)
p1 <- pie(tra.genes[order(count)][cell.type.tag == 'TetPos', count],labels=NA)
print(p1)
dev.off()

# TRB
trb.genes <- clonotypes[!is.na(trb.v) & !is.na(trb.j), .(count=.N), by=.(cell.type.tag,trb.v,trb.j)]

tmp.output.file <- paste0(reports.dir,
  '/figures/VDJUsagePieChart_Chain-TRB_Tet-Negative.pdf')
pdf(tmp.output.file)
tmp.p <- pie(trb.genes[order(count)][cell.type.tag == 'TetNeg', count],labels=NA)
print(tmp.p)
dev.off()

tmp.output.file <- paste0(reports.dir,
  '/figures/VDJUsagePieChart_Chain-TRB_Tet-Positive.pdf')
pdf(tmp.output.file)
p1 <- pie(trb.genes[order(count)][cell.type.tag == 'TetPos', count],labels=NA)
print(p1)
dev.off()

### --------------------- Downstream analysis: Tet-positive data set --------------------- ###

# ---> Dot plot
tmp.output.file <- paste0(reports.dir, 
  '/figures/DotPlot_Resolution-',
  str_remove(res.of.int,'RNA_snn_res.'),'.pdf')
pdf(tmp.output.file, width = 13, height = 4)
d1 <- DotPlot(object=srt.obj.tpos, features=genes.of.int, group.by=res.of.int, dot.scale = 12, col.min = -1.5, col.max = 1.5) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(d1)
dev.off()

# ---> UMAP with colored clonotypes.
meta.data <- as.data.table(cbind(
  srt.obj.tpos@meta.data,
  srt.obj.tpos@reductions$umap@cell.embeddings
))

meta.data[clonotypes.of.int, on=.(clonotype.tag), col.tag:=color.code]
meta.data[is.na(col.tag),col.tag:='#C0C0C0']
meta.data[,col.tag:=factor(col.tag,levels=clonotypes.of.int$color.code)]
meta.data <- meta.data[order(col.tag, decreasing=T)]

# Selected clonotypes bigger, as triangles, and less transparent than background clonotypes. Background clonotypes as hollow dots.
meta.data[col.tag == '#C0C0C0', size := "1"]
meta.data[col.tag != '#C0C0C0', size := "2"]
meta.data[, size := as.factor(size)]
tmp.ggplot <- ggplot(data=meta.data, aes(x=UMAP_1, y=UMAP_2, col=col.tag, size=size, alpha=size, shape=size)) +
  geom_point(stroke = 2) +
  scale_size_manual(values = c("1" = 4, "2" = 6), guide = FALSE) +
  scale_alpha_manual(values = c("1" = 0.7, "2" = 0.9), guide = FALSE) +
  scale_shape_manual(values = c("1" = 1, "2" = 19), guide = FALSE) +
  scale_color_manual(values=clonotypes.of.int$color.code,
    labels=clonotypes.of.int$name) +
  scale_y_continuous(expand=expansion(add=0.2)) + 
  scale_x_continuous(expand=expansion(add=0.2)) +
  labs(x='UMAP 1', y='UMAP 2', col='Clonotype:') +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.key.size = unit(2, "lines"),
  legend.key.width = unit(2, "lines"),
  legend.key.height = unit(2, "lines"),
  text = element_text(size = 16),
  panel.background = element_blank(),
  axis.line = element_line(color = "black"),
  legend.background = element_blank(),
  legend.key = element_blank())

tmp.file.name <- 'UMAP-SelectedClonotypes-Bigger_Background-Hollow_All'
plotTiff(file.name=tmp.file.name, 
  output.dir=paste0(reports.dir,'/figures'), 
  ggplot.obj=tmp.ggplot, 
  blank.element=blank.complement.3)

# ---> UMAP with clonotype affinities.
meta.data[clonotypes.of.int, on=.(clonotype.tag),
avg.affinity:=as.numeric(avg.affinity)]

# Color scale: green to red.
# Selected clonotypes bigger, and less transparent than background clonotypes. Background clonotypes as hollow dots.
tmp.ggplot <- ggplot(data=meta.data, aes(x=UMAP_1, y=UMAP_2, col=avg.affinity, size=size, alpha=size, shape=size)) +
  geom_point(stroke=2) +
  scale_size_manual(values = c("1" = 4, "2" = 6), guide = FALSE) +
  scale_alpha_manual(values = c("1" = 0.7, "2" = 0.9), guide = FALSE) +
  scale_shape_manual(values = c("1" = 1, "2" = 19), guide = FALSE) +
  scale_color_gradientn(colours=green2red(1024), na.value='#C0C0C0') +
  scale_y_continuous(expand=expansion(add=0.2)) + 
  scale_x_continuous(expand=expansion(add=0.2)) +
  labs(x='UMAP 1', y='UMAP 2', col='Average affinity:') +
  theme(legend.key.size = unit(2, "lines"),
  legend.key.width = unit(2, "lines"),
  legend.key.height = unit(2, "lines"),
  text = element_text(size = 16),
  panel.background = element_blank(),
  axis.line = element_line(color = "black"))

tmp.file.name <- 'UMAP_Heatmap_SelectedClonotypes-Bigger_Background-Hollow_All_GR'
plotTiff(file.name=tmp.file.name, output.dir=paste0(reports.dir,'/figures'), ggplot.obj=tmp.ggplot, blank.element=blank.complement.3)

