####example code
## this code is designed to reproduce the figures shown in the associated paper 
## 'Morphodynamic styles: characterising the behaviour of gravel-bed rivers using a new, quantitative index'
## by W.H. Booker and B.C. Eaton
## Data provided includes pre-processed (volumetric) data with which to produce plots,
## raw data used to derive products are available upon request.

require(magicaxis)


dir = 'C:/Users/willt/Documents/phd_work/papers/methods'

##Figure 4

f_sed_out = read.csv(paste0(dir, '/fix_sed.csv'), header = T)
m_sed_out = read.csv(paste0(dir, '/mob_sed.csv'), header = T)


plot(f_sed_out[,1], f_sed_out[,2], type = 'o',
     axes = F, ann = F, col = 'red', ylim = c(0, 0.0015), xlim = c(0,960), pch = 0, cex = 0.8, log = '', lwd =1.5)
points(m_sed_out[,1], m_sed_out[,2], type = 'o',col = 'blue', pch = 0, cex = 0.8, lwd =1.5)
legend('topright', col = c('red','blue'),
       legend = c('Fixed', 'Mobile'),
       pch = c(0,0), bty = 'n',lty = c(rep(1,2)), lwd =1.5)
magaxis(frame.plot = T,  xlab = 'Time (min)', ylab = expression(paste(q[b], ' ( ',m^2, ' mi', n^{-1},')')))


##Figure 5
#time step now between successive DEMs

f_dem_out = read.csv(paste0(dir, '/fix_dem.csv'), header = T)
m_dem_out = read.csv(paste0(dir, '/mob_dem.csv'), header = T)

#delta v
plot(f_dem_out[,1], f_dem_out[,3], type = 'o',
     axes = F, ann = F, col = 'red', ylim = c(-0.0008, 0.0001), xlim = c(0,960), pch = 0, cex = 0.8, log = '', lwd =1.5)
points(cumsum(m_dem_out[,1]), m_dem_out[,3], type = 'o',col = 'blue', pch = 0, cex = 0.8, lwd =1.5)
legend('bottomright', col = c('red','blue'),
       legend = c('Fixed','Mobile'),
       pch = c(0,0), bty = 'n',lty = c(rep(1,2)), lwd =1.5)
magaxis(frame.plot = T,  xlab = 'Time (min)', ylab = expression(paste(Delta,'v', ' ( ',m^2, ' mi', n^{-1},')')))
#m
plot(f_dem_out[,1], f_dem_out[,4], type = 'o',
     axes = F, ann = F, col = 'red', ylim = c(0, 0.003), xlim = c(0,960), pch = 0, cex = 0.8, log = '', lwd =1.5)
points(m_dem_out[,1], m_dem_out[,4], type = 'o',col = 'blue', pch = 0, cex = 0.8, lwd =1.5)
legend('topright', col = c('red','blue'),
       legend = c('Fixed','Mobile'),
       pch = c(0,0), bty = 'n',lty = c(rep(1,2)), lwd =1.5)
magaxis(frame.plot = T,  xlab = 'Time (min)', ylab = expression(paste('m', ' ( ',m^2, ' mi', n^{-1},')')))


##Figure 6
# throughput ratio; volumetric qb/morphologic activity

f_zeta = f_dem_out[,2]/f_dem_out[,4]
m_zeta = m_dem_out[,2]/m_dem_out[,4]

plot(f_dem_out[,1], f_zeta, type='o', col = 'red',
     axes = F, ann = F, ylim = c(0, 10), xlim = c(0,960), pch = 0, cex = 0.8, log = '', lwd =1.5)
points(m_dem_out[,1], m_zeta, type='o', col = 'blue', pch = 0, cex = 0.8, lwd =1.5)
legend('topleft', col = c('red','blue'),
       legend = c('Fixed','Mobile'),
       pch = c(0,0), bty = 'n',lty = c( rep(1,2)), lwd =1.5)
magaxis(frame.plot = T,  xlab = 'Time (min)', ylab = expression(zeta))

##Figure 7
# violin plots of m
require(vioplot)
require(scales)

f_vio = read.csv(paste0(dir, '/fix_vio_data.csv'), header = T)
m_vio = read.csv(paste0(dir, '/mob_vio_data.csv'), header = T)

vioplot(f_vio[,-1], xlim =c(1,24), at = 1:24, col=alpha('red', 0.5), yaxt = 'n', ylim = c(0,7.967313e-07))
magaxis(side = c(2),ylab = expression(paste('m (', m^2,' ', min^{-1},')')), labels = c(T))
axis(1, labels= (f_dem_out[,1]),at = 1:24, las = 2)
mtext(text ='Time (min)', side = 1, line = 3)

vioplot(m_vio[,-1], xlim =c(1,32), at = 1:32, col=alpha('blue', 0.5), yaxt = 'n', ylim = c(0,7.967313e-07))
magaxis(side = c(2),ylab = expression(paste('m (', m^2,' ', min^{-1},')')), labels = c(T))
axis(1, labels= (m_dem_out[,1]),at = 1:32, las = 2)
mtext(text ='Time (min)', side = 1, line = 3)

##Figure 8-11
#Code necessary to spatialise a DEM is shown; dod4

fig_8_dem = raster(paste0(dir, '/fig_8.tif'))
fig_8_mat = as.matrix(fig_8_dem)
na.sum = sum(is.na(values(fig_8_dem)))
not.na = ncell(dod_list[[1]]) - na.sum
ares = not.na/1000000
width = not.na/ncol(dod_list[[1]])/1000

xsection_v = apply(fig_8_mat, 2, function(x){sum(((x)*0.001*0.001*(1-0.54)/width), na.rm = T)})
xsection_m = apply(fig_8_mat, 2, function(x){sum((abs(x)*0.001*0.001*(1-0.54)/width), na.rm = T)})

#calculate total morphologic activity
M = sum((abs(fig_8_mat)*0.001*0.001*(1-0.54)/width), na.rm = T)
#throughput ratio for this DoD
z = f_zeta[4]
#number of cross sections
n = dim(fig_8_mat)[2]
#calc. observed Qb per xsection using z/xsection_m; multiply by expected ratio of total M over number of observations
zeta_ratio = z/xsection_m*M/n
#calc. spatial version of throughput ratio
spatial_zeta = zeta_ratio/(mean(zeta_ratio,na.rm=T)/z)