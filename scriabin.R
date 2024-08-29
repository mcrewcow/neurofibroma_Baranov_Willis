
#load Seurat object
deeann_control1 <- NormalizeData(deeann_control)
deeann_control_ALRA <- SeuratWrappers::RunALRA(deeann_control1)
deeann_control_ALRA_ip <- FindAllInteractionPrograms(deeann_control_ALRA, iterate.threshold = 300, group.by = "cell.type", assay = "alra", sim_threshold = 0.4)
deeann_control_ALRA_ip_sig <- InteractionProgramSignificance(deeann_control_ALRA_ip, n.replicate = 500)

ip_pvals <- deeann_control_ALRA_ip_sig %>% as_tibble() %>%
  dplyr::select(name,ends_with("pval")) %>% unique() %>%
  pivot_longer(!name, names_to = "cell.type", values_to = "pval") %>%
  group_by(name) %>% dplyr::mutate(min_p = min(pval)) %>%
  dplyr::select(name,min_p) %>% unique() %>% 
  dplyr::filter(min_p < 0.05) %>% pull(name)

library(dplyr)
deeann_control_ALRA_ip_sig %>% dplyr::filter(name %in% ip_pvals)

deeann_control_ALRA <- ScoreInteractionPrograms(deeann_control_ALRA, deeann_control_ALRA_ip_sig)

panc_id_ip_lig <- as.matrix(deeann_control_ALRA[["IPligands"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = deeann_control_ALRA$cell.type) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))
library(ComplexHeatmap)
library(cowplot)

Heatmap(panc_id_ip_lig, show_column_names = F, name = "Ligands")

panc_id_ip_rec <- as.matrix(deeann_control_ALRA[["IPreceptors"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = deeann_control_ALRA$cell.type) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))
Heatmap(panc_id_ip_rec, show_column_names = F, name = "Receptors")

act_stellate_ip <- panc_id_ip_lig["Schwann cells",]
poi <- gsub("ligands_","",names(which(act_stellate_ip==max(act_stellate_ip))))

#Seurat's FeaturePlot has a nice option to blend expression of two features together on the same plot
IPFeaturePlot(deeann_control_ALRA, ip = poi)
DimPlot(deeann_control_ALRA, group.by = "cell.type", label = T, repel = T) + NoLegend()

moi <- reshape2::melt(deeann_control_ALRA_ip_sig %>% dplyr::filter(name==poi) %>%
                        select("lr_pair",contains("connectivity"))) %>% arrange(-value)
moi$lr_pair <- factor(moi$lr_pair, levels = unique(moi$lr_pair))
ggplot(moi, aes(x = lr_pair, y = value, color = variable)) + 
  geom_point() + theme_cowplot() + ggpubr::rotate_x_text() + labs(x = NULL, y = "Intramodular\nconnectivity")

table(deeann_control_ALRA$cell.type)
deeann_control_ALRA <- RenameIdents(deeann_control_ALRA,"Pericytes" = 'rest','Endothelial cells 2'= 'rest', 'Endothelial cells 1'= 'rest', 'Fibroblasts'= 'rest',
                                    'Keratinocytes'= 'rest', 'Myeloid cells'= 'rest', 'Melanocytes'= 'rest' )

deeann_control_ALRA$cell.type1 <- deeann_control_ALRA@active.ident

ip_by_celltype <- IPCellTypeSummary(deeann_control_ALRA, group.by = "cell.type1")
ip_by_celltype %>% group_by(sender) %>% top_n(n = 1, wt = additive.score)

soi <- "Schwann cells" #sender of interest
roi <- "rest"#receiver of interest
scriabin::load_nichenet_database()
variant_genes <- IDVariantGenes(deeann_control_ALRA, group.by = "cell.type1", assay = 'alra')
gene_signature <- GenerateCellSignature(deeann_control_ALRA, variant_genes = variant_genes)

active_ligands <- RankActiveLigands(deeann_control_ALRA, signature_matrix = gene_signature, assay = 'alra')


TopLigandsByIdent(deeann_control_ALRA, active_ligands = active_ligands, 
                  sender = soi, receiver = roi, group.by = "cell.type1", assay = 'alra')

receiver_cells <- colnames(deeann_control_ALRA)[deeann_control_ALRA$cell.type1==roi]
# calculates the predicted target genes within a set of receiver cells (the ductal cells) 
PlotLigandTargetAlluvium(deeann_control_ALRA, signature_matrix = gene_signature,
                         active_ligands = active_ligands, receiver_cells = receiver_cells,
                         ligands_of_interest = c("TGFA","THPO","ADIPOQ",'SPP1','CXCL8','NTNG1','SLIT2','PLAT'))

roi <- "Schwann cells" #sender of interest
soi <- "rest"#receiver of interest
scriabin::load_nichenet_database()
variant_genes <- IDVariantGenes(deeann_control_ALRA, group.by = "cell.type1", assay = 'alra')
gene_signature <- GenerateCellSignature(deeann_control_ALRA, variant_genes = variant_genes)

active_ligands <- RankActiveLigands(deeann_control_ALRA, signature_matrix = gene_signature, assay = 'alra')


TopLigandsByIdent(deeann_control_ALRA, active_ligands = active_ligands, 
                  sender = soi, receiver = roi, group.by = "cell.type1", assay = 'alra')

receiver_cells <- colnames(deeann_control_ALRA)[deeann_control_ALRA$cell.type1==roi]
# calculates the predicted target genes within a set of receiver cells (the ductal cells) 
PlotLigandTargetAlluvium(deeann_control_ALRA, signature_matrix = gene_signature,
                         active_ligands = active_ligands, receiver_cells = receiver_cells,
                         ligands_of_interest = c("LAMA2","TGFB1","APP",'CALR','IGF1','NLGN3','FGF2','DKK2'))



gc()
deeann_treated1 <- NormalizeData(deeann_treated)
deeann_treated_ALRA <- SeuratWrappers::RunALRA(deeann_treated1)
deeann_treated_ALRA_ip <- FindAllInteractionPrograms(deeann_treated_ALRA, iterate.threshold = 300, group.by = "cell.type", assay = "alra", sim_threshold = 0.4)
deeann_treated_ALRA_ip_sig <- InteractionProgramSignificance(deeann_treated_ALRA_ip, n.replicate = 500)

ip_pvals <- deeann_treated_ALRA_ip_sig %>% as_tibble() %>%
  dplyr::select(name,ends_with("pval")) %>% unique() %>%
  pivot_longer(!name, names_to = "cell.type", values_to = "pval") %>%
  group_by(name) %>% dplyr::mutate(min_p = min(pval)) %>%
  dplyr::select(name,min_p) %>% unique() %>% 
  dplyr::filter(min_p < 0.05) %>% pull(name)

library(dplyr)
deeann_treated_ALRA_ip_sig %>% dplyr::filter(name %in% ip_pvals)

deeann_treated_ALRA <- ScoreInteractionPrograms(deeann_treated_ALRA, deeann_treated_ALRA_ip_sig)

panc_id_ip_lig <- as.matrix(deeann_treated_ALRA[["IPligands"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = deeann_treated_ALRA$cell.type) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))
library(ComplexHeatmap)
library(cowplot)



Heatmap(panc_id_ip_lig, show_column_names = F, name = "Ligands", row_order = c('Fibroblasts',
          'Schwann cells', 'Melanocytes','Keratinocytes', 'Pericytes', 'Myeloid cells',
          'Endothelial cells 2', 'Endothelial cells 1'))

panc_id_ip_rec <- as.matrix(deeann_treated_ALRA[["IPreceptors"]]@data %>% t() %>%
                              as.data.frame() %>% add_column(celltype = deeann_treated_ALRA$cell.type) %>%
                              group_by(celltype) %>%
                              summarise_if(is.numeric, mean) %>% column_to_rownames("celltype"))
Heatmap(panc_id_ip_rec, show_column_names = F, name = "Receptors", row_order = c(
                                                                                 'Schwann cells','Fibroblasts', 'Pericytes', 'Melanocytes','Keratinocytes',  'Myeloid cells',
                                                                                 'Endothelial cells 2', 'Endothelial cells 1'))

act_stellate_ip <- panc_id_ip_lig["Schwann cells",]
poi <- gsub("ligands_","",names(which(act_stellate_ip==max(act_stellate_ip))))

#Seurat's FeaturePlot has a nice option to blend expression of two features together on the same plot
IPFeaturePlot(deeann_treated_ALRA, ip = poi)
DimPlot(deeann_treated_ALRA, group.by = "cell.type", label = T, repel = T) + NoLegend()

moi <- reshape2::melt(deeann_treated_ALRA_ip_sig %>% dplyr::filter(name==poi) %>%
                        select("lr_pair",contains("connectivity"))) %>% arrange(-value)
moi$lr_pair <- factor(moi$lr_pair, levels = unique(moi$lr_pair))
ggplot(moi, aes(x = lr_pair, y = value, color = variable)) + 
  geom_point() + theme_cowplot() + ggpubr::rotate_x_text() + labs(x = NULL, y = "Intramodular\nconnectivity")

table(deeann_treated_ALRA$cell.type)
deeann_treated_ALRA <- RenameIdents(deeann_treated_ALRA,"Pericytes" = 'rest','Endothelial cells 2'= 'rest', 'Endothelial cells 1'= 'rest', 'Fibroblasts'= 'rest',
                                    'Keratinocytes'= 'rest', 'Myeloid cells'= 'rest', 'Melanocytes'= 'rest' )

deeann_treated_ALRA$cell.type1 <- deeann_treated_ALRA@active.ident

ip_by_celltype <- IPCellTypeSummary(deeann_treated_ALRA, group.by = "cell.type1")
ip_by_celltype %>% group_by(sender) %>% top_n(n = 1, wt = additive.score)

soi <- "Schwann cells" #sender of interest
roi <- "rest"#receiver of interest
scriabin::load_nichenet_database()
variant_genes <- IDVariantGenes(deeann_treated_ALRA, group.by = "cell.type1", assay = 'alra')
gene_signature <- GenerateCellSignature(deeann_treated_ALRA, variant_genes = variant_genes)

active_ligands <- RankActiveLigands(deeann_treated_ALRA, signature_matrix = gene_signature, assay = 'alra')


TopLigandsByIdent(deeann_treated_ALRA, active_ligands = active_ligands, 
                  sender = soi, receiver = roi, group.by = "cell.type1", assay = 'alra')

receiver_cells <- colnames(deeann_treated_ALRA)[deeann_treated_ALRA$cell.type1==roi]
# calculates the predicted target genes within a set of receiver cells (the ductal cells) 
PlotLigandTargetAlluvium(deeann_treated_ALRA, signature_matrix = gene_signature,
                         active_ligands = active_ligands, receiver_cells = receiver_cells,
                         ligands_of_interest = c("TGFA","ADIPOQ",'CXCL8','AGT','THBS1','BMP7','IL22','CCL17'))

roi <- "Schwann cells" #sender of interest
soi <- "rest"#receiver of interest
scriabin::load_nichenet_database()
variant_genes <- IDVariantGenes(deeann_treated_ALRA, group.by = "cell.type1", assay = 'alra')
gene_signature <- GenerateCellSignature(deeann_treated_ALRA, variant_genes = variant_genes)

active_ligands <- RankActiveLigands(deeann_treated_ALRA, signature_matrix = gene_signature, assay = 'alra')


TopLigandsByIdent(deeann_treated_ALRA, active_ligands = active_ligands, 
                  sender = soi, receiver = roi, group.by = "cell.type1", assay = 'alra')

receiver_cells <- colnames(deeann_treated_ALRA)[deeann_treated_ALRA$cell.type1==roi]
# calculates the predicted target genes within a set of receiver cells (the ductal cells) 
PlotLigandTargetAlluvium(deeann_treated_ALRA, signature_matrix = gene_signature,
                         active_ligands = active_ligands, receiver_cells = receiver_cells,
                         ligands_of_interest = c("APP","CD200","TGFB2",'APLN','IL19','CALR','DKK2','MSTN'))


