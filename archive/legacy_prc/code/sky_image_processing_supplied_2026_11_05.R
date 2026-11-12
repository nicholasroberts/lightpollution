#################################################################################################


library(dplyr)
library(ggplot2)
library(reshape)
library(tiff)
library(fields)
library(imager)
library(magrittr)
library(mmand)
library(jpeg)
library(raster)
library(RNiftyReg)
library(RColorBrewer)

#################################################################################################

#################################################################################################

#################################################################################################




DoP <- function(z,f,n,o) {
  S.1 <- (z - n)/(z+n)
  S.2 <- (f-o)/(f+o)
  d <- sqrt(S.1^2+S.2^2)
  return(d)
}

AoP <- function(z,f,n,o) {
  S.1 <- (z - n)/(z+n)
  S.2 <- (f-o)/(f+o)
  da <- 0.5*atan2(S.2, S.1)
  return(da)
}


ecovis_colors <- c(
  `red`        = "#FF0000",
  `endred`        = "#550000",
  `darkgreen`   = "#00FF1B",
  `brightgreen` = "#80FF00",
  `brightblue`   = "#00FFFF",
  `darkblue`   = "#002EFF",
  `endblue` = "#002151",
  `orange`     = "#FF8000",
  `yellow`     = "#FFF300",
  `light grey` = "#cccccc",
  `dark grey`  = "#8c8c8c",
  `zero` = "#f6e653",
  `fifteen` = "#98d25e",
  `thrity` = "#5ab080",
  `fortyfive` = "#458b8c",
  `sixty` = "#3e6289",
  `seventyfive` = "#42347a",
  `ninety` = "#3e6289",
  `oneofive` = "#7e2f9f",
  `onetwoone` = "#a73c88",
  `onethreefive` = "#c65c6f",
  `onefifty` = "#e68b57",
  `onesixfive` = "#f5c14e",
  `oneeighty` = "#f6e653",
  `neg` = "#091be2",
  `less` = "#020a66",
  `ze` = "#1b1b1c",
  `more` = "#7a7517",
  `pos` = "#efe62f")


#' Function to extract ecovis colors as hex codes
#'
#' @param ... Character names of ecovisvis_colors 
#'
ecovis_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (ecovis_colors)
  
  ecovis_colors[cols]
}


ecovis_palettes <- list(
  `main`  = ecovis_cols("blue", "brightgreen", "red"),
  
  `Cont`  = ecovis_cols("brightgreen", "neg", "less", "ze", "more", "pos", "red"),
  
  `AoP`   = ecovis_cols("zero", "fifteen", "thirty", "fortyfive",
                        "sixty", "seventyfive", "ninety", "oneofive",
                        "onetwoone", "onethreefive", "onefifty", "onesixfive",
                        "oneeighty"),
  
  `DoP` = ecovis_cols("endblue", "darkblue", "brightblue", "yellow", "orange", "red", "endred"),
  
  `grey`  = ecovis_cols("light grey", "dark grey")
)

#' Return function to interpolate a ecovis color palette
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'
ecovis_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- ecovis_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}


#' Color scale constructor for ecovis colors
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_ecovis <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- ecovis_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("ecovis_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for ecovis colors
#'
#' @param palette Character name of palette in ecovis_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_ecovis <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- ecovis_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("ecovis_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}


#################################################################################################
# SET WORKING DIRECTORY


#################################################################################################
############## Images one
one <- readTIFF("1.tif")
two <- readTIFF("4.tif")
three <- readTIFF("7.tif")
four <- readTIFF("10.tif")

#three <- readTIFF("10.tif")
#two <- readTIFF("1.tif")
#one <- readTIFF("4.tif")
#four <- readTIFF("7.tif")



#extra values
extra2  <- 0
extra3 <- 0
extra4 <- 0

trans.tw <- applyTransform(buildAffine(angles=-(pi/4)+(extra2), anchor="centre", source=two[,,3]), two[,,3])
trans.th <- applyTransform(buildAffine(angles=-(pi/2)+(extra3), anchor="centre", source=three[,,3]), three[,,3])
trans.fo <- applyTransform(buildAffine(angles=-((3*pi)/4)+(extra4), anchor="centre", source=four[,,3]), four[,,3])

################################################################################################
## extra section to check the extra rotation in the initial allinment
#### its the extra 1, extra2, extra3 values that set this. 
#plot(one[,,3] %>% as.cimg)
#plot(S0.1 %>% as.cimg)
#plot(trans.three$image %>% as.cimg)
#plot(trans.four$image %>% as.cimg)
################################################################################################


trans.two <- niftyreg(trans.tw, one[,,3])
trans.three <- niftyreg(trans.th, one[,,3])
trans.four <- niftyreg(trans.fo, one[,,3])

forward.two <- (forward(trans.two))
forward.three <- (forward(trans.three))
forward.four <- (forward(trans.four))

trans.two <- applyTransform(forward.two, trans.tw)
trans.three <- applyTransform(forward.three, trans.th)
trans.four <- applyTransform(forward.four, trans.fo)

im.1.p <- DoP(one[,,3], trans.two, trans.three, trans.four)
im.1.a <- AoP(one[,,3], trans.two, trans.three, trans.four)

im.1.p[is.na(im.1.p)] <- 0
im.1.p[im.1.p > 0.95] <- 0

S0.1 <- (one[,,3] + trans.three)
S1.1 <- (one[,,3] - trans.three)

#im.1.a[is.na(im.1.a)] <- 0



###################

#second set 

############## 
one <- readTIFF("2.tif")
two <- readTIFF("5.tif")
three <- readTIFF("8.tif")
four <- readTIFF("11.tif")

#one <- readTIFF("5.tif")
#three <- readTIFF("11.tif")
#four <- readTIFF("8.tif")
#two <- readTIFF("2.tif")


trans.tw <- applyTransform(buildAffine(angles=-(pi/4)+(extra2), anchor="centre", source=two[,,3]), two[,,3])
trans.th <- applyTransform(buildAffine(angles=-(pi/2)+(extra3), anchor="centre", source=three[,,3]), three[,,3])
trans.fo <- applyTransform(buildAffine(angles=-((3*pi)/4)+(extra4), anchor="centre", source=four[,,3]), four[,,3])

trans.two <- applyTransform(forward.two, trans.tw)
trans.three <- applyTransform(forward.three, trans.th)
trans.four <- applyTransform(forward.four, trans.fo)


im.2.p <- DoP(one[,,3], trans.two, trans.three, trans.four)
im.2.a <- AoP(one[,,3], trans.two, trans.three, trans.four)

S0.2 <- (one[,,3] + trans.three)
S1.2 <- (one[,,3] - trans.three)

im.2.p[is.na(im.2.p)] <- 0
im.2.p[im.2.p > 0.95] <- 0

#im.2.a[is.na(im.2.a)] <- 0


###################

#third set 

############## 
one <- readTIFF("3.tif")
two <- readTIFF("6.tif")
three <- readTIFF("9.tif")
four <- readTIFF("12.tif")

#one <- readTIFF("6.tif")
#two <- readTIFF("3.tif")
#three <- readTIFF("12.tif")
#four <- readTIFF("9.tif")


trans.tw <- applyTransform(buildAffine(angles=-(pi/4)+(extra2), anchor="centre", source=two[,,3]), two[,,3])
trans.th <- applyTransform(buildAffine(angles=-(pi/2)+(extra3), anchor="centre", source=three[,,3]), three[,,3])
trans.fo <- applyTransform(buildAffine(angles=-((3*pi)/4)+(extra4), anchor="centre", source=four[,,3]), four[,,3])

trans.two <- applyTransform(forward.two, trans.tw)
trans.three <- applyTransform(forward.three, trans.th)
trans.four <- applyTransform(forward.four, trans.fo)

im.3.p <- DoP(one[,,3], trans.two, trans.three, trans.four)
im.3.a <- AoP(one[,,3], trans.two, trans.three, trans.four)

im.3.p[is.na(im.3.p)] <- 0
im.3.p[im.3.p > 0.95] <- 0

#im.3.a[is.na(im.3.a)] <- 0

S0.3 <- (one[,,3] + trans.three)
S1.3 <- (one[,,3] - trans.three)

#########
## Average all the images

im.p <- (im.1.p + im.2.p + im.3.p)/3.0
im.a <- (180/pi)*(im.1.a + im.2.a + im.3.a)/3.0

P0 <- (S0.1+S0.2+S0.3)/3
P1 <- (S1.1+S1.2+S1.3)/3


##### SET THE DR
a_par <- 0.4
a_perp <- 0.04

sent.r <- 10^(-a_par)
sent.q <- 10^(-a_perp)

perp <- P0-(0.5*(sent.q*(P0-P1)+sent.r*(P0+P1)))
par <- P0-(0.5*(sent.q*(P0+P1)+sent.r*(P0-P1)))

int.im <- par+perp
pol.im <- (par-perp)
#plot(int.im %>% as.cimg)
#plot(pol.im %>% as.cimg)
######

#PLot DoP
boats <- im.p %>% as.cimg
df <- boats %>% as.data.frame
p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=value))
p <- p + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())
p <- p + scale_fill_ecovis(palette = "DoP", discrete = FALSE)
p
ggsave("DoP_000.pdf", width=18, height=11.5, units = "cm")

#Plot AoP
boats <- im.a %>% as.cimg
df <- boats %>% as.data.frame
p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=value))
p <- p + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())
p <- p + scale_fill_ecovis(palette = "AoP", discrete = FALSE)
p
ggsave("angle_0000.pdf", width=18, height=11.5, units = "cm")

#PLot Receptor Contrast
boats <- pol.im %>% as.cimg
df <- boats %>% as.data.frame
p <- ggplot(df,aes(x,y))+geom_raster(aes(fill=value))
p <- p + scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())
p <- p + scale_fill_scico(palette = 'tofino', limits = c(-0.6,0.6))
p  
ggsave("Contrast_0000.pdf", width=18, height=11.5, units = "cm")


################################################################################################


