source("analysis.R")

# Parameters
colorbars = 'green'
barsize = 2
fontsize = 12
pvx = 'left'
pvy = 'bottom'

## chlororaphis
dataset = 'chlororaphis'
folder = paste0('/tmp/rarefan_test/', dataset, '/out')
treefile = paste0(dataset, '.nwk')


plotREPINs(folder, treefile, 0, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 1, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 2, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 3, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 4, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 5, colorbars, barsize, fontsize)

plotCorrelationSingle(folder, 0, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 1, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 2, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 3, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 4, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 5, theme,  fontsize, pvx, pvy)

drawRAYTphylogeny(folder)

## neisseria
dataset = 'neisseria'
folder = paste0('/tmp/rarefan_test/', dataset, '/out')
treefile = paste0(dataset, '.nwk')


plotREPINs(folder, treefile, 0, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 1, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 2, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 3, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 4, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 5, colorbars, barsize, fontsize)

plotCorrelationSingle(folder, 0, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 1, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 2, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 3, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 4, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 5, theme,  fontsize, pvx, pvy)

drawRAYTphylogeny(folder)
## dokdonia
dataset = 'dokdonia'
folder = paste0('/tmp/rarefan_test/', dataset, '/out')
treefile = paste0(dataset, '.nwk')


plotREPINs(folder, treefile, 0, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 1, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 2, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 3, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 4, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 5, colorbars, barsize, fontsize)

plotCorrelationSingle(folder, 0, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 1, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 2, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 3, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 4, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 5, theme,  fontsize, pvx, pvy)

drawRAYTphylogeny(folder)



## ecoli
dataset = 'ecoli'
folder = paste0('/tmp/rarefan_test/', dataset, '/out')
treefile = paste0(dataset, '.nwk')


plotREPINs(folder, treefile, 0, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 1, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 2, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 3, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 4, colorbars, barsize, fontsize)
plotREPINs(folder, treefile, 5, colorbars, barsize, fontsize)

plotCorrelationSingle(folder, 0, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 1, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 2, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 3, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 4, theme,  fontsize, pvx, pvy)
plotCorrelationSingle(folder, 5, theme,  fontsize, pvx, pvy)

drawRAYTphylogeny(folder)
