##### TRAITEMENT DATA COUCHE DE MELANGE KERGUELEN  // VINCENT DORIOT // avr. 2023
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_arg <- args[startsWith(args, file_arg)]
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_arg[1]))))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }

  getwd()
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))
setwd(project_dir) ### REPERTOIRE CONTENANT LES DATA + FILES ...
input_dir <- file.path(project_dir, "processed", "1_gridded_data", "r_input")
output_dir <- file.path(project_dir, "processed", "2_fPCA", "GLORYS_CL_masked")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data=read.table(file.path(input_dir, "GLORYS_CL_masked.txt"),header=TRUE, na.strings=c("NA", ""))

# data[4980,6]=0.7075137 #### Y'A UNE DONNEE MANQUANTE... MOYENNE DES AUTRES....
data[data == ''] <- NA
rownames(data) <- 1:nrow(data) 
#### ORGANISATION DES DATA
#### DATES PAR an.mois

tps=paste(rep(unique(data$year),each=12),c(paste(0,1:9,sep=""),10:12),sep=".")
nmonth=length(tps)

datam=by(data[,5],rep(tps,each=1521),list)
names(datam)=tps

### ON GARDE 2007->2023 SOIT 18x12=216 MOIS
# data=data[-c(1:12996),]
# datam=datam[-c(1:36)]

tps=paste(rep(unique(data$year),each=12),c(paste(0,1:9,sep=""),10:12),sep=".")
nmonth=length(tps)

#### ON ORGANISE LES DONNEES DANS UNE MATRICE 1521 LIGNES X (12X18) COLONNES
###SERIE TEMPORELLE PAR CELLULE
#x11()
st_cell=NULL
for(k in 1:1521){
  y=NULL
  for(j in 1:length(tps)){
    y=c(y,datam[[j]][k])
  }
  st_cell=rbind(st_cell,y)
  # obs=which(is.na(y))
  # if(length(obs)!=192){
  #   plot(dates,y,main=k,xlim=c(2004,2023))
  #   abline(v=2004:2023,lty=2)
  # }
  #text(locator(1),"")
}


############# PARTIE I : TRAVAIL SUR LA SERIE TEMPORELLE TOTALE / CELLULE
### DETECTION DES SERIE TEMPS VIDES = CELLULES VIDES
ligNA=NULL
for(k in 1:1521){
  ligNA=c(ligNA,sum(is.na(st_cell[k,])))
}
indNA=which(ligNA==nmonth)

#### ON CREE DEUX MATRICES AVEC DES NAs : LES DATES ET LES DATAS (QUI SONT DEJA DANS st_cell)
dates=seq(2007,2023,length=nmonth+1)
dates=dates[-(nmonth+1)]
datemat=matrix(dates,1521,nmonth,byrow=TRUE)

#### CASE indNA COMPLETEMENT VIDE
donmat=st_cell
if (length(indNA) > 0) {
  datemat=datemat[-indNA, , drop=FALSE]
  donmat=donmat[-indNA, , drop=FALSE]
}

library(fdapace)
### ACP NON PARAM
datemat=t(datemat)
donmat=t(donmat)

set.seed(1)

res=FPCA(data.frame(donmat),data.frame(datemat)
         ,optns=list(useBinnedData='OFF',
                     methodMuCovEst='smooth',
                     userBwCov=0.25,
                     userBwMu=0.05,
                     kernel="epan",
                     nRegGrid=204,
                     dataType="Sparse",
                     methodSelectK = "FVE",
                     FVEthreshold = 1,
                     methodXi='CE',
                     maxK=100,
                     plot=FALSE
         ))

#x11()
#plot(t(datemat),t(donmat),pch=39,col="grey")
#lines(res$workGrid,res$mu,type="l",col=2,lwd=2)
#abline(v=2007:2023,lty=2)

#### RECONSTITUTION SERIE TEMP /CELL
####A FAIRE TOURNER OBLIGATOIREMENT POUR AVOIR LA SUITE NOTAMMENT Xestot
nbcp=length(res$lambda)
rx=c(2007,2023)
ry=range(donmat,na.rm=TRUE)
nobs=ncol(donmat)
#x11()
Xestot=NULL
### CI-DESSOUS : DEGAGER LES # POUR OBTNIR LES PLOTS DE MANIERE DYNAMIQUE
for(numobs in 1:nobs){
  Xest=res$mu+apply(sweep(res$phi,2,res$xiEst[numobs,],"*")[,1:nbcp],1,sum)
  Xestot=cbind(Xestot,Xest)
  
  #
  #plot(datemat[,numobs],donmat[,numobs],main="Reconstitution temporelle d'un bin après l'utilisation du kriging"
   #    ,cex.main=2.5,cex.axis=2,cex.lab=2,cex=2.5,lwd=10,
    #   xlab="Années",ylab="MLD (m)",xlim=rx,ylim=c(50,200),col="darkgrey",pch=39)
  #matpoints(res$workGrid,Xconf,lty=c(2,1,2),type="l",col="dark blue")
  #mtext("",side = 3, line = -3, outer = TRUE)
  #text(locator(1), "")
}

#### RECONSTITUTION DES IMAGES PAR DATE

#x11()
for(tps in 1:204){
  XestotNA=matrix(-99.99,1521,1)
  XiestNA=matrix(-99.99,1521,1)
  cpt=0
  for(k in 1:1521){ ### 353 EST LE NB DE CELL NON-NA POUR UN TEMPS tps
    cpt=cpt+1
    dec=which(indNA==k)
    if(length(dec)!=0){
      XestotNA[k]=NA
      XiestNA[k]=NA
      cpt=cpt-1
    }
    if(length(dec)==0){
      XestotNA[k]=Xestot[tps,cpt]
      XiestNA[k]=res$xiEst[cpt,2]
    }
  }
  #pics=matrix(XestotNA,39,39,byrow=TRUE)
  #image(unique(data[,1]),unique(data[,2]),pics)
  #text(locator(1),"")
}


################################### MOVIE
graphics.off()

# saveVideo({
#   for(tps in 1:204){
#     XestotNA=matrix(-99.99,1521,1)
#     cpt=0
#     for(k in 1:1521){ ### 353 EST LE NB DE CELL NON-NA POUR UN TEMPS tps
#       cpt=cpt+1
#       dec=which(indNA==k)
#       if(length(dec)!=0){
#         XestotNA[k]=NA
#         cpt=cpt-1
#       }
#       if(length(dec)==0){
#         XestotNA[k]=Xestot[tps,cpt]
#       }
#     }
#     pics=matrix(XestotNA,39,39,byrow=TRUE)
#     image(unique(data[,1]),unique(data[,2]),pics,main=dates[tps],col=colo)
#   }
#   
# },video.name = "mld_ker.mp4", other.opts = "-pix_fmt yuv420p -b 204k")

#### RECONSTITUTION DES IMAGES PAR PCs
##x11()
quelpc=10
XiestNA=matrix(-99.99,1521,1)
cpt=0
for(k in 1:1521){ ### 353 EST LE NB DE CELL NON-NA POUR UN TEMPS tps
  cpt=cpt+1
  dec=which(indNA==k)
  if(length(dec)!=0){
    XiestNA[k]=NA
    cpt=cpt-1
  }
  if(length(dec)==0){
    XiestNA[k]=res$xiEst[cpt,quelpc]
  }
}
# pics=matrix(XiestNA,39,39,byrow=TRUE)
# image(unique(data[,1]),unique(data[,2]),pics)

# #x11()
# matplot(res$workGrid,res$phi[,1:k], type ="l",lty=c("solid","dashed","dotted")
#         ,lwd=c(3,3,3),xlab="Années",ylab="Valeurs de PC", col=c("black","blue","red"),
#         ,main="3 premiers modes de variabilité de la MLD"
#         ,cex.main=2.5,cex.axis=2,cex.lab=2) 
# legend("topright", legend=c("1","2","3"), col=c("black","blue","red"),
#        lty=c("solid","dashed","dotted"), horiz=F,cex=2,lwd=c(3,3,3))

write.csv(res$workGrid, file.path(output_dir, "PCA_GRID.csv"), row.names=TRUE)
write.csv(res$lambda, file.path(output_dir, "PCA_LAMBDA.csv"), row.names=TRUE)
write.csv(Xestot, file.path(output_dir, "PCA_MLD.csv"), row.names=TRUE)
write.csv(res$phi, file.path(output_dir, "PCA_PHI.csv"), row.names=TRUE)
write.csv(res$xiEst, file.path(output_dir, "PCA_XIEST.csv"), row.names=TRUE)
write.csv(res$mu, file.path(output_dir, "PCA_MU.csv"), row.names=TRUE)
