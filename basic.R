load("G:/sampleComparisonEnv.RData")

#change the clusters name from metadata
RenameIdents(samples.combined, '0'=	'Fibroblasts 1',
             '1'=	'Macrophage 1',
             '2'=	'Myelinating SC 1',
             '3'=	'Neurons 1',
             '4'=	'Fibroblasts 2',
             '5'=	'Non-myelinating SC',
             '6'=	'Macrophage 2',
             '7'=	'Neurons 1',
             '8'=	'Satellite cells',
             '9'=	'Macrophage 2',
             '10'=	'Macrophage 3',
             '11'=	'Myelinating SC 2',
             '12'=	'Perineurial cells',
             '13'=	'Endothelial cells',
             '14'=	'Pericytes',
             '15'=	'Neurons 2',
             '16'=	'Pericytes',
             '17'=	'Endoneurial fibroblasts 1',
             '18'=	'Endoneurial fibroblasts 2',
             '19'=	'Neurons 3',
             '20'=	'Neurons 4',
             '21'=	'SC precursors',
             '22'=	'T cells',
             '23'=	'Macrophage 4',
             '24'=	'Perineurial cells',
             '25'=	'Fibroblasts 1',
             '26'=	'MDSCs',
             '27'=	'Neurons 5',
             '28'=	'Neurons 6',
             '29'=	'Neurons 4',
             '30'=	'Proliferating cells'
)

DimPlot(samples.combined, label = T, repel = T, label.box = T, split.by = 'orig.ident') #no new populations from pretumor to tumor found
