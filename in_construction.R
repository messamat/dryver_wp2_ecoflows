library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('R/packages.R')
source("R/functions.R")

hydromod_present_dir <- file.path('data', 'wp1', 'Results_present_period_final', 'data')
bio_dir <- file.path('data', 'wp2', '01_WP2 final data')
#datdir <- file.path('data', 'data_annika')
resdir <- 'results'

tar_option_set(format = "qs")

metacols <- c('campaign', 'site', 'running_id', 'country', 'date',
              'summa_sample', 'sample.id', 'sample_id', 'organism')
drn_dt <- data.table(
  country = c("Croatia", "Czech", "Finland", "France",  "Hungary", "Spain"),
  catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal")
)

hydro_combi <- expand.grid(
  in_country = drn_dt$country,
  in_varname =  c('isflowing', 'qsim'),
  stringsAsFactors = FALSE)

############### PARAMETERS

in_country <- 'France'
in_network_path = tar_read(network_directed_gpkg_list)[[in_country]]

############### Function
#At complex confluences (when more than two line segments flow into a node, 
#i.e. more than two inflows to an outflow), new confluences of only two streams 
#are created:
#1. complex confluences are found based on the fact that the outflow has more 
#than two previous streams 
#2. the inflow with the shortest cumulative length from its source is found; 
#the end of this segment will be moved 
#3. the inflow with the smallest angle to this inflow is found; 
#this segment will be cut into tow segments close to the junction 
#using the GRASS function v.edit(tool = break) creating a new confluence 
#4. the shortest inflow found in 2 is moved to the newly created confluence 
#using v.edit(tool = vertexmove) 
#5. all lengths are updated (segment length, cumulative length, 
#i.e. length of the stream from the source, distance to the outlet). 

#The distance the shortest confluence is moved depends on the number of inflows. 
#For three inflows, it is moved 1/12 time the DEM cellsize upstream, 
#for seven (the extremely rare maximum) 5/12 * cellsize. 

#Read input network
rivnet <- st_read(rivnet_path) %>%
  st_cast("LINESTRING") %>%
  #Make sure that the geometry column is equally named regardless 
  #of file format (see https://github.com/r-spatial/sf/issues/719)
  st_set_geometry('geometry') 



######################################################################################################
# FIXES MULTIPLE CONFLUENCES IN STREAM NETWORKS. #
##################################################
# Author: Miguel Porto
# All nodes which have >2 streams flowing to it are corrected. The outermost streams' end vertices
# are adjusted by "step" meters along the downgoing stream. New nodes are created at suitable places,
# and existing lines suitably split.
#*** Requires a line shapefile with the proper FROM_NODE and TO_NODE fields. The file is assumed to be
#*** correct in all other aspects, there is no error checking.
###### USAGE EXAMPLE
# rios=readOGR("streams_Pt.shp","streams_Pt")
# correctedshp=fix.streams(rios,step=10)
# writeOGR(correctedshp,"streams_corrected.shp","streams_corrected","ESRI Shapefile")
######################################################################################################
#https://github.com/miguel-porto/fix-streams/blob/master/fix-streams.r

library(rgdal)
fix.streams=function(shp,from="FROM_NODE",to="TO_NODE",step=10) {
  # step is the desired length (in map units) by which the river sinks are adjusted (separated) downstream.
  pieces = list()
  probrivers = list()
  removeindexes = integer(0)
  CRS = shp@proj4string
  
  # find multiple confluence nodes
  nv=table(shp@data[,to])
  mc=as.numeric(names(nv[nv>2]))	# these are the nodes with >2 rivers flowing to	them
  maxnode=max(c(shp@data[,to],shp@data[,from]))	# max node ID, for creating new nodes
  cat(length(mc),"multiple confluences found.\n")
  cat("Cutting lines and tweaking vertices...\n");flush.console()
  for(i in mc) {	# for each problematic node
    #msrc=shp[shp@data[,to] %in% i,]	# get source rivers flowing to it
    privers=which(shp@data[,to] %in% i)	# get source rivers flowing to it (problematic rivers)
    sinkriver=which(shp@data[,from] %in% i)	# get sink river
    msrc=shp[privers,]
    mto=shp[sinkriver,]
    delta=sum(shp@data[,to] %in% i)-2	# how many problematic rivers flow to there? leave only two, the others correct
    newnodes=c(mto@data[,from],(maxnode+1):(maxnode+delta),mto@data[,to])	# the IDs of the nodes that will be created (the first remains the same for the two "good" rivers)
    
    coo=coordinates(mto)[[1]][[1]]	# coordinates of the sink river
    # order the rivers by their angle, so that the outermost rivers are first adjusted (alternating the side)
    v1=coo[1,]-coo[2,]
    v2=matrix(nc=2,nr=length(msrc))
    for(j in 1:length(msrc)) {
      tmp=coordinates(msrc[j,])[[1]][[1]]
      v2[j,]=tmp[dim(tmp)[1]-1,]-tmp[dim(tmp)[1],]
    }
    ang=atan2(v2[,2],v2[,1])-atan2(v1[2],v1[1])
    ang[ang>0]=ang[ang>0] %% pi
    ang[ang<0]=-((-ang[ang<0]) %% pi)
    names(ang)=privers
    ang=ang[order(abs(ang),decreasing=T)]
    if(sum(ang>0)>sum(ang<0)) {
      angneg=c(as.numeric(names(ang[ang<0])),rep(NA,sum(ang>0)-sum(ang<0)))
      angpos=as.numeric(names(ang[ang>=0]))
    } else {
      angpos=c(as.numeric(names(ang[ang>=0])),rep(NA,sum(ang<0)-sum(ang>0)))
      angneg=as.numeric(names(ang[ang<0]))
    }
    angs=matrix(c(angpos,angneg),nc=2)
    privers=na.omit(as.vector(t(angs)))[1:delta]
    
    # split sink river in delta pieces plus the remainder		
    tmplines=split.line(coo,delta,step)
    if(length(tmplines)<=delta) stop("You must decrease step: river ",sinkriver)
    for(j in 1:length(tmplines)) {
      # cut sink river into pieces, as many as necessary
      rid=runif(1,10^6,10^7)	# random ID for the piece
      # create a new piece with the j'th (step+1) vertices of the sink river
      piece=SpatialLines(list(Lines(list(Line(tmplines[[j]])),rid)),proj4string=CRS)
      newdata=mto@data
      newdata[1,]=NA
      newdata[1,from]=newnodes[j]
      newdata[1,to]=newnodes[j+1]
      rownames(newdata)=rid
      pieces=c(pieces,list(SpatialLinesDataFrame(piece,newdata,match=F)))	# save pieces for later use
      if(j>1) {
        # now change coords of problematic rivers
        pri=privers[j-1]	# pick the j'th problematic river
        tmp=coordinates(shp[pri,])[[1]][[1]]
        tmp[dim(tmp)[1],]=tmplines[[j]][1,]	# change the coordinate of the last vertex of problematic river
        tmp1=SpatialLines(list(Lines(list(Line(tmp)),shp[pri,]@lines[[1]]@ID)),proj4string=CRS)	# keep same ID (original will be removed)
        tmp1=SpatialLinesDataFrame(tmp1,shp[pri,]@data)
        tmp1@data[1,to]=newnodes[j]
        probrivers=c(probrivers,list(tmp1))	# collect new rivers to replace old
      }
    }
    
    removeindexes=c(removeindexes,c(privers[1:delta],sinkriver))
    maxnode=maxnode+delta
  }
  cat("Now reassembling shape...\n");flush.console()
  newlines=pieces[[1]]
  for(j in 2:length(pieces)) {
    newlines=rbind(newlines,pieces[[j]])
  }
  for(j in 1:length(probrivers)) {
    newlines=rbind(newlines,probrivers[[j]])
  }
  
  # remove all problematic + sink rivers
  newshp=shp[-removeindexes,]
  newshp=rbind(newshp,newlines)
  return(newshp)
}

split.line=function(line,n,length,debug=F) {
  # splits a Lines object (or coordinate matrix) in n segments of length length (starting in the begining) plus the remaining segment (what is left)
  if(debug) plot(line)
  coo=coordinates(line)
  if(inherits(coo,"list")) {
    if(length(coo)>1) {
      stop("Multiple lines not allowed")
    } else {
      coo=coo[[1]]
      if(!inherits(coo,"matrix")) {
        if(!inherits(coo,"list")) stop("Invalid line object")				
        if(length(coo)>1) stop("Multiple lines not allowed")
        coo=coo[[1]]
      }
    }
  } else {
    if(!inherits(coo,"matrix")) stop("Invalid line object")
  }
  
  pieces=list()	
  accum=0
  i=1
  remainder=0
  newcoords=matrix(nc=2,nr=0)	
  repeat {
    newcoords=rbind(newcoords,coo[i,])
    v=c(coo[i+1,]-coo[i,])
    hyp=sqrt(sum(v^2))
    accum=accum+hyp
    if(accum>length && length(pieces)<n) {	# cut in this segment
      oript=coo[i,]
      #			accum=hyp
      repeat {
        newpt=oript+v/hyp*(length-remainder)
        newcoords=rbind(newcoords,newpt)
        pieces=c(pieces,list(newcoords))
        newcoords=matrix(newpt,nr=1)
        remainder=0
        if(debug) points(newpt[1],newpt[2],pch=19)
        accum=accum-length
        if(accum<length || length(pieces)>=n) break
        oript=newpt
      }
      remainder=accum
    } else remainder=accum
    
    i=i+1
    if(i>=dim(coo)[1]) {
      newcoords=rbind(newcoords,coo[dim(coo)[1],])
      pieces=c(pieces,list(newcoords))
      break
    }
  }
  if(debug) {
    for(i in 1:length(pieces)) {
      pieces[[i]][,1]=pieces[[i]][,1]-20*i
      lines(pieces[[i]],col="red")
      points(pieces[[i]][1,1],pieces[[i]][1,2])
    }
  }
  return(pieces)
}


