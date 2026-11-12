Xcc <- function(im, cen) Xc(im) - width(im)/2 + cen
Ycc <- function(im, cen) Yc(im) - height(im)/2 + cen

##### go through the ROIs 4 degree acptance over the DRA array

source("/Users/bznwr/OneDrive - University of Bristol/Research/Research Projects/2.Visual Space/wrapper.R")
library(plyr)

#################################################################################

site <- "Gavarnie"

par <- paste0(site,"_par.txt")  #open the par and perp files
perp <- paste0(site,"_perp.txt")


### calcualtes the angles of the rhabdoms across the DRA from back to front

phi1 <- c(seq(-pi/4,0, by = (pi/4)/(22))) 
phi2 <- c(seq(0,pi/2, by = (pi/2)/(10)))

phi <- rep(c(phi1,phi2), 3)

########################################
#read in both images
im <- read.table(par, header=T, sep=",")
bpar <- im %>% as.cimg

imp <- read.table(perp, header=T, sep=",")
bperp <- imp %>% as.cimg

####################################################
#finds the diamter of the image to set the centre

#diam <- (2*(width(bo)/2)/3)  #if landscape
diam <- (width(bpar)/2)^2 # in portrait
########################

########################################
#set the coords of the  x and y centres of ecach of the 100 omatidia
## the ((sqrt(diam))/6) is for auto set the coords to go out to just past 30 degrees for the FOV
xcentres <- ((sqrt(diam))/6)*c(rep(seq(-1,3, by = 4/(33)), 3))  
ycentres <- c(rep(-10, 36), rep(-35, 34), rep(-60, 32))
#ycentres <- c(rep(-70, 36), rep(-95, 34), rep(-120, 32))

########################################
plot(bperp)
## rotate to allign with the band - SET THE ANGLE
#bpar <-  imrotate(bpar, 90)
#bperp <- imrotate(bperp, 90)


#plot(bpar)


###################################################
#Test code to plot and visualize the DRA on the image

cont <- im
## calculate receptor contrasts
cont$value <- log(im$value/imp$value)


p <- ggplot(cont,aes(x,y))+geom_raster(aes(fill=value))
p <- p + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())
p <- p + scale_fill_ecovis(palette = "Cont", discrete = FALSE, limits=c(-2.302585,2.302585)) + coord_fixed()
p

############# Plt FOVs on image
#cont <- imrotate(cont %>% as.cimg, -8)
plot(bpar)

#cf <- bpar

for (i in 1:101) {
ppar <- (Xcc(bpar, xcentres[i])^2+Ycc(bpar, ycentres[i])^2) < diam/200
highlight(ppar) }

####################################################
#loop to calcualte the receptro output of each of the parallel photoreceptors

fli <- c(1:102)
vc <- laply(fli, function(ff){
  
  ppar <- (Xcc(bpar, xcentres[ff])^2+Ycc(bpar, ycentres[ff])^2) < diam/200
  pperp <- (Xcc(bperp, xcentres[ff])^2+Ycc(bperp, ycentres[ff])^2) < diam/200
  
 
  ppa <- (bpar*ppar) %>% as.data.frame
  ppe <- (bperp*pperp) %>% as.data.frame
  
  pa <- subset(ppa, ppa$value != 0)
  pe <- subset(ppe, ppe$value != 0)


  pa$value <- pa$value*cos(phi[ff])^2+pe$value*sin(phi[ff])^2
  pe$value <- pa$value*sin(phi[ff])^2+pe$value*cos(phi[ff])^2

  
  ta <- sum(ppar)*1.818182
  
  poe <- sum(pa$value)/(ta)
  pee <- sum(pe$value)/(ta)
  prc <- log(poe/pee)
  
  ##entropy
  sumd <- (poe+pee)
  ent.sd <- poe/sumd
  ent.ld <- pee/sumd
  inf <- -((ent.sd * logb(ent.sd, 2))+(ent.ld * logb(ent.ld, 2)))
  return(c(prc, inf))
})


################################# 
## plot the values with the correct colour scheme

dx <- 6*cos(phi)
dy <- 6*sin(phi)



plo_hx <- data.frame(x = xcentres, y = ycentres, vc = vc[,1])

p <- ggplot(plo_hx, aes(x, y, fill=vc)) + geom_hex(stat="identity")
p <- p + scale_x_continuous(expand=c(-1.1,1.1))+scale_y_continuous(expand=c(-2.5,2.5),trans=scales::reverse_trans())
p <- p + scale_fill_ecovis(palette = "Cont", discrete = FALSE, limits=c(-2.302585,2.302585)) + coord_fixed()
p <- p +geom_segment(aes(x=x-dx/2, y=y-dy/2, xend=x+dx, yend=y+dy, color = "white")) 
p

ggsave(paste0(site,"_DRAcontrast_plot.png"))

fi_t <- file.path(paste0(site,"_DRAcontrast_Array.txt"))
write.table(plo_hx, file = fi_t, sep = ',', quote = F, row.names = F)


abs(max(vc[,1])-min(vc[,1]))

##############################################################################
####################################################

##Do the infromation

plo_in <- data.frame(x = xcentres, y = ycentres, vc = vc[,2])

p <- ggplot(plo_in, aes(x, y, fill=vc)) + geom_hex(stat="identity")
p <- p + scale_x_continuous(expand=c(-1.1,1.1))+scale_y_continuous(expand=c(-2.5,2.5),trans=scales::reverse_trans())
p <- p + scale_fill_viridis(option = "magma", name = "Entropy", alpha = 1.0, limits = c(0.92, 1),# make fill a bit brighter
                   discrete = F, # discrete classes, thus guide_legend instead of _colorbar
                   direction = -1, # dark is lowest, yellow is highest
                   guide = guide_legend(
                     keyheight = unit(5, units = "mm"),
                     title.position = "top",
                     reverse = T)) + coord_fixed()
p <- p +geom_segment(aes(x=x-dx/2, y=y-dy/2, xend=x+dx, yend=y+dy, color = "white")) 
p

ggsave(paste0(site,"_DRAinformation_plot.png"))

fi_t <- file.path(paste0(site,"_DRAinformation_Array.txt"))
write.table(plo_hx, file = fi_t, sep = ',', quote = F, row.names = F)



(max(vc[,2])-min(vc[,2]))/(max(vc[,2])+min(vc[,2]))












##############################################################################
####################################################

##########
###want to get the vectro field of the DRA correct


library(grid)

xcentres <- 100*c(rep(seq(-1,3, by = 4/(33)), 3))  
ycentres <- c(rep(-10, 36), rep(-35, 34), rep(-60, 32))


phi1 <- c(seq(-pi/4,0, by = (pi/4)/(22)))
phi2 <- c(seq(0,pi/2, by = (pi/2)/(10)))

phi <- rep(c(phi1,phi2), 3)

dx <-6*cos(phi)
dy <- 6*sin(phi)
df <- data.frame(x=xcentres,y=ycentres,dx = dx, dy = dy)

ggplot(data=df, aes(x=x-dx/2, y=y-dy/2)) + 
  geom_segment(aes(xend=x+dx, yend=y+dy)) + coord_fixed()














 #################################################################################

#################################################################################

#################################################################################

#################################################################################
#################################################################################

### signal

path <- getwd()
filenames <- list.files(path)
filenames <- filenames[-c(1,5)]



sigg <- c(0.274, 0.439, 0.419, 0.470, 0.165, 0.193, 0.149, 0.321, 0.192, 0.262,
          0.276, 0.455, 0.423, 0.157, 0.122, 0.123, 0.144, 0.236)

DRAc <- c(0.407, 0.603, 0.692, 0.239, 0.159, 0.093, 0.448, 0.206, 0.378, 0.539,
          0.620, 0.629, 0.142, 0.211, 0.071, 0.109, 0.216)

ENT <- c(0.015, 0.0323, 0.042, 0.0051, 0.002, 0.0007, 0.018, 0.004, 0.013, 0.026,
         0.034, 0.035, 0.002, 0.004, 0.0004, 0.001, 0.006)

Sig.data <- data.frame(name = filenames, signal = sigg)  # spider
Sig.data <- data.frame(name = filenames, signal = DRAc) # DRA PRC moth
Sig.data <- data.frame(name = filenames, signal = ENT) # DRA entropy moth



ppr = laply(filenames, function(filename){
  dap <- read.csv(filename)
  dar <- na.omit(dap)
  t.data <- ddply(dar, .(), summarise,
                  me = mean(value),
                  std = sd(value) 
  )
  dan <- c(t.data$me, t.data$std)
  return(dan)
}
)
 

final.data <- cbind(Sig.data, ppr)
colnames(final.data) <- c("name, total", "signal", "rad", "sd")


p<- ggplot(final.data, aes(x=(rad), y=signal)) + geom_point() +
  geom_errorbarh(aes(xmin=rad-sd, xmax=rad+sd))
p

o <- order(final.data$rad)
df_o <- final.data[o, ] # order it by increasing x


fit <- nls(signal ~ SSasymp(rad, yf, y0, log_alpha), data = df_o)
qplot(rad, signal, data = augment(fit)) + geom_line(aes(y = .fitted)) + 
  xlim(c(0,70)) + ylim(c(0,0.05))


ggplot(final.data, aes(x=(rad), y=signal)) + geom_point() +
  geom_errorbarh(aes(xmin=rad-sd, xmax=rad+sd)) + 
  xlim(c(0,70)) + ylim(c(0,0.05))




##################################################################

###


# HEX BIN

##########################################################


made <- data.frame(x = 1, y = 1, vc = 0.4386080)

p <- ggplot(made, aes(x, y, fill=vc)) + geom_hex(stat="identity")
p <- p + scale_x_continuous(expand=c(0,1))+scale_y_continuous(expand=c(0,1),trans=scales::reverse_trans())
p <- p + scale_fill_ecovis(palette = "Cont", discrete = FALSE, limits=c(-2.302585,2.302585)) + coord_fixed()
p

