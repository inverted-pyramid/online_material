# The code below replicates the results presented in: 
# Cahen-Fourot, L., Campiglio, E., Dawkins, E., Godin, A., and Kemp-Benedict, E. (2019) 
# "Looking for the inverted pyramid: An application using input-output networks"
# For correspondence, please write to: Emanuele.campiglio@wu.ac.at

# The code is structured in four sections:
# 1. Introductory statements
# 2. Calculation of sectoral forward linkages
# 3. Definition and plotting of the inverted pyramids
# 4. Write-up of results

######################################################################################
# 1. Introductory statements
######################################################################################

# Load required libraries
library(iotables)
library(igraph)
library(plotrix)
library(writexl)
# Load the functions written in the library file
source('Inverted_pyramids_function_library.R')

# Countries to be analysed
countries<- c("BE","BG", "CZ", "DE", "EL", "ES", "FR", "HR", "IE", "IT", "CY", "LV", "HU", "AT", "PL", "SI", "SE", "UK")
# For reference: entire set of countries: countries<- c("BE","BG", "CZ", "DE", "EL", "ES", "FR", "HR", "IE", "IT", "CY", "LV", "HU", "AT", "PL", "SI", "SE", "UK")

# List of NACE sectors
sectors<-c(  "A01",    "A02",    "A03",    "B",      "C10-12", "C13-15", "C16",    "C17",    "C18",    "C19",    "C20",    "C21",    "C22",    "C23",    "C24",    "C25",    "C26",    "C27",    "C28",    "C29",    "C30",    "C31_32", "C33" ,   "D",      "E36",    "E37-39", "F",      "G45",    "G46",    "G47",    "H49",    "H50",    "H51",    "H52",    "H53",    "I",      "J58",    "J59_60", "J61",    "J62_63", "K64",    "K65",    "K66",    "L68A",   "L68B",   "M69_70", "M71",    "M72",   "M73",    "M74_75", "N77",    "N78",    "N79",    "N80-82", "O",      "P",      "Q86",    "Q87_88", "R90-92", "R93",    "S94",    "S95",    "S96")

# Choose sector at the bottom of the pyramid
# Note: if A01, A02 or A03 are selected, Ireland (IE) needs to be excluded from the country list (it doesn't report those sectors)
sector <- "B"

# We define the dataframes where to store forward linkage results
# Normalised forward linkages
int_fln<-data.frame(matrix(0, nrow=63, ncol=length(countries), dimnames=list(sectors,countries)))
# Sector rankings in terms of forward linkages
int_fln_rank<-data.frame(matrix(0, nrow=63, ncol=length(countries), dimnames=list(NULL, countries)))
# Sector positions in national ranking
int_sect_fln_pos<-data.frame(matrix(0, nrow=63, ncol=length(countries), dimnames=list(sectors, countries)))
# GDP values
int_GDP<-array(0, dim=length(countries), dimnames = list(countries))
# Sectors on the first layer of the inverted pyramid
neighbours<-data.frame(matrix(0, nrow=63, ncol=length(countries), dimnames=list(NULL, countries)))
# Sector positions in inverted pyramid (order of layer)
int_sect_layer<-data.frame(matrix(0, nrow=63, ncol=length(countries), dimnames=list(sectors, countries)))
sect_layer<-array(0, dim=length(sectors), dimnames = list(sectors))

# This is the start of the country for loops. It runs the same code for all the countries defined above. 
count<-1 
for (geo in countries) {

######################################################################################
# 2. Forward linkages
######################################################################################
  
# We use the 'iotables' package to download IO table from Eurostat
# Package ocumentation: https://cran.r-project.org/web/packages/iotables/iotables.pdf
# We define the following default settings: source: "naio_10_cp1700" (Symmetric input-output table at basic prices (product by product)); unit: "MIO_EUR" (million Euros); labelling: "short"
# Choose type of table (DOM or TOTAL)
stk_flow<-"TOTAL"
# Download IO table
io_dat<-as.data.frame(iotable_get(labelled_io_data = NULL, source = "naio_10_cp1700", geo = geo, year = 2010, unit = "MIO_EUR", stk_flow = stk_flow, labelling = "short", data_directory = NULL,
                      force_download = TRUE))

# Transform the first column into row names
io_dat[,1]<-as.vector(io_dat[,1])
rownames(io_dat)<-io_dat[,1]
io_dat<-io_dat[,-1]

#Calculate GDP as Total Value Added at Basic Prices (B1G, TOTAL) plus Net Taxes on Products (D21X31, TU) (see Eurostat Manual, p.305)
GDP<-io_dat["B1G","TOTAL"]+io_dat["D21X31","TU"]

#Some manipulation to make sure that the matrix is consistent (ie. Eurostat NA data is downloaded in the rows but not in the columns, creating symmetry issues)
ncol_ind<-which(colnames(io_dat)=="TOTAL")-1
nrow_ind<-which(rownames(io_dat)=="TOTAL")-1
if (ncol_ind!=nrow_ind) {
  redundant_rows=vector()
  for (i in 1:nrow_ind){
    if (rownames(io_dat)[i]!=colnames(io_dat)[i]){
      redundant_rows<-cbind(redundant_rows,i)
    }
  }
  io_dat<-io_dat[-redundant_rows,]    
}

#In case the T and U columns still exist, we eliminate them
for(i in 1:ncol_ind){
  for(removesect in c("CPA_T", "CPA_U")){
    if (colnames(io_dat)[i]==removesect){
      removencol<-i
      io_dat<-io_dat[,-removencol]
      io_dat<-io_dat[-removencol,]
    }
  }
}
#Redefine ncol_ind and nrow_ind
ncol_ind<-which(colnames(io_dat)=="TOTAL")-1
nrow_ind<-which(rownames(io_dat)=="TOTAL")-1

# In case there is some sector with all NAs, remove those sectors
removencol<-vector()
for(i in 1:ncol_ind){
  if(is.na(io_dat["P1",i])){
    removencol<-cbind(removencol,i)
  }
}
if (length(removencol)!=0){
  io_dat<-io_dat[,-removencol]
  io_dat<-io_dat[-removencol,]
}
#Redefine ncol_ind and nrow_in in case something was deleted in the previous passages
ncol_ind<-which(colnames(io_dat)=="TOTAL")-1
nrow_ind<-which(rownames(io_dat)=="TOTAL")-1

# Reformulation of codes in the inter-industry matrix to eliminate the 'CPA_' part and improve readability of pyramid charts
sect_codes<-gsub("CPA_", "", rownames(io_dat)[1:ncol_ind])
rownames(io_dat)[1:ncol_ind]<-sect_codes
colnames(io_dat)[1:ncol_ind]<-sect_codes

#Subset the interindustry matrix by selecting the first 63*63 elements (up to sector S96; T and U excluded)
io_dat_ind<-io_dat[1:ncol_ind,1:ncol_ind]

# Extract output to use in the linkage calculations
# Output is equal to "P1" (Output) when looking at the DOM table, to "TS_BP" (Total Supply at Basic Prices) when looking at the TOTAL table
if (stk_flow=="DOM"){
  output_dat<-io_dat['P1',1:ncol_ind]
} else {
  output_dat<-io_dat['TS_BP',1:ncol_ind]
}

#transform them in matrices/numeric vectors (otherwise no matrix multiplications possible for Ghosh)
io_mat_ind<-as.matrix(io_dat_ind)
output_mat<-as.numeric(output_dat)

#Find Ghosh matrix 
B<-solve(diag(output_mat))%*%io_mat_ind
G<-solve(diag(nrow(B))-B)

# Calculate forward linkages
# fl are forward linkages, fln are normalised forward linkages (fl divided by the mean of fl)
# The _rank vectors put the sectors in decreasing order of forward linkages
fl<-colSums(t(G))
fln<-fl/mean(fl)
names(fl)<-names(output_dat)
names(fln)<-names(fl)
fl_rank<-fl[order(fl, decreasing=T)]
fln_rank<-fln[order(fln, decreasing=T)]
fln_rank_codes<-names(fln_rank)
fln_rank_fin<-vector()
for (i in 1:length(fln_rank)){
  fln_rank_fin[i]<-paste(fln_rank_codes[i]," (", round(fln_rank[i], digits = 3),")", sep = "")
}

# Create vector with positions of sectors in the ranking
sect_fln_pos<-vector()
for (i in sectors){
  fln_pos<-which(names(fln_rank)==i)
  sect_fln_pos<-c(sect_fln_pos, fln_pos)
}
names(sect_fln_pos)<-names(fln)

#The lines below deal with the fact that certain countries do not report certain sectors. NA is inserted in their place.
for (i in 1:length(sectors)){
  if (sectors[i]!=names(fln)[i]){
    fln<-c(fln[1:i-1],NA,fln[i:length(fln)])
    names(fln)[i]<-sectors[i]
  }
  if (sectors[i]!=names(sect_fln_pos)[i]){
    sect_fln_pos<-c(sect_fln_pos[1:i-1],NA,sect_fln_pos[i:length(sect_fln_pos)])
    names(sect_fln_pos)[i]<-sectors[i]
  }
}

######################################################################################
# 3. Inverted pyramids
######################################################################################

#Create minimal fully connected matrix
io_mfc<-mfc(io_mat_ind)

# The lines below create a directed graph from the interindustry matrix; gets rid of loops and multiple edges (simplify); 
# and creates a table of the number of edges connected to each sector (treated as a vertex). This is done through the degree function. 
# We distinguish between 'in' and 'out' edges. 

# Select a method (either "mfc" or "q")
method<-"q"
if (method=="mfc"){
  ip.graph <- simplify(graph_from_adjacency_matrix(as.matrix(io_mfc), mode = "directed", weighted = TRUE))
  m.inout <- data.frame(cbind(degree(ip.graph, mode = "in"),degree(ip.graph, mode = "out")))
  names(m.inout) <- c("in", "out")
} else {
  colnames(G)<-rownames(G)
  G.graph <- simplify(graph_from_adjacency_matrix(as.matrix(G), mode = "directed", weighted = TRUE))
  ip.graph <- cascade_network(G.graph, sector, q = 0.1 , depth=NA)
}


# We define the layers ('shells') of the pyramid (see explanation in the ego_shells function declaration)
# In ego_shells we declare the graph we want to examine, the sector at the bottom of the pyramid, the type of connections to examine ("in", "out" or "all")
shells <- ego_shells(ip.graph, sector, "out")
# 
# # Now we want to assign each vertex (sector) to its shell
# # We define shndx as a vector of length equal to the number of vertices (gorder function returns the number of vertices of a graph)
# # Then assign to each element the number of the shell where the corresponding sector lies.
 shndx <- vector(mode = "integer", length = gorder(ip.graph))
 for (i in 1:length(shells)) {
   shndx[as.integer(V(ip.graph)[shells[[i]]])] <- i
 }
 for (i in 1:length(V(ip.graph))) {
   V(ip.graph)$layer[i]<-shndx[i]
 }

for (i in 1:length(sectors)){
  if (sectors[i] %in% names(V(ip.graph))){
    sect_layer[i]<- V(ip.graph)$layer[V(ip.graph)$name==sectors[i]]
  } else {
    sect_layer[i]<- NA
  }
}
 
## Register the ellipse shape with igraph
add_shape("ellipse", clip=shapes("circle")$clip,
          plot=myellipse)

# assign colors
#colbar <- rainbow(length(shells))
# create layout #V for Vertex; E for Edge
V(ip.graph)$size <- 3.3
V(ip.graph)$shape <- "ellipse"
V(ip.graph)$color <- "white"
V(ip.graph)$label.color <- "black"
V(ip.graph)$label.family <- "sans"
V(ip.graph)$label.cex <-1.3
V(ip.graph)$label.dist <- 0
#V(ip.graph)$frame.color <- "gray40"
#V(ip.graph)$label.degree <-pi
E(ip.graph)$arrow.size <- 0.0
E(ip.graph)$color <- "gray60"
E(ip.graph)$curved <- 0

ll <- ego_layout(ip.graph, sector, "out", jitter=0.075, inverted = T, nmax = 11)
plot(ego_directed(ip.graph, sector, "out"), layout=ll)
title(main = geo)


######################################################################################
# 4. Write-up of results
###################################################################################### 

# Assign results to the dataframes summarising international results
int_fln[count]<-fln 
int_fln_rank[count]<-c(fln_rank_fin,rep(NA,length(sectors)-length(fln_rank_fin)))
int_sect_fln_pos[count]<-sect_fln_pos
int_GDP[count]<-GDP
neighbours[count]<-c(names(shells[[2]]), rep(NA,length(sectors)-length(names(shells[[2]]))))
int_sect_layer[count]<-sect_layer

count<-count+1
}
# This is the end of the country for loop


#Create mean columns and attach them to int_fln
if (length(countries)>1){
  fwd_mean<-rowMeans(int_fln[,1:length(int_fln)],na.rm=TRUE)
  GDP_share<-int_GDP/sum(int_GDP)
  fwd_mean_weighted<-vector()
  for (i in 1:nrow(int_fln)){
    fwd_mean_weighted[i]<-weighted.mean(int_fln[i,],GDP_share, na.rm = TRUE)
  }
  names(fwd_mean_weighted)<-rownames(int_fln)

  fwd_mean_rank<-fwd_mean[order(fwd_mean, decreasing=T)]
  fwd_mean_weighted_rank<-fwd_mean_weighted[order(fwd_mean_weighted, decreasing=T)]

  fwd_mean_rank_fin<-vector()
  for (i in 1:length(fwd_mean_rank)){
    fwd_mean_rank_fin[i]<-paste(names(fwd_mean_rank)[i]," (", round(fwd_mean_rank[i], digits = 3),")", sep = "")
  }
  fwd_mean_weighted_rank_fin<-vector()
  for (i in 1:length(fwd_mean_rank)){
    fwd_mean_weighted_rank_fin[i]<-paste(names(fwd_mean_weighted_rank)[i]," (", round(fwd_mean_weighted_rank[i], digits = 3),")", sep = "")
  }

  int_fln<-cbind(int_fln, fwd_mean, fwd_mean_weighted)
  int_fln_rank<-cbind(int_fln_rank, fwd_mean_rank_fin, fwd_mean_weighted_rank_fin)
}

# Recreate a first column in int_fln with sector names for exporting to Excel
int_fln<-cbind(sectors, int_fln)
int_sect_fln_pos<-cbind(sectors, int_sect_fln_pos)


# Code to study layer positioning of sectors
# int_layer_res summarise sector positioning internationally
# int_mode_layers returns the mode of sector positioning
int_layer_res<-data.frame(matrix(0, nrow=63, ncol=8, dimnames=list(sectors, c("l1","l2","l3","l4","l5","l6","l7","NA"))))
for (sect in sectors){
  for (l in 1:(ncol(int_layer_res))){
    int_layer_res[sect,l]<-sum(int_sect_layer[sect,]==l,na.rm=TRUE)
  }
  int_layer_res[sect,8]<-sum(is.na(int_sect_layer[sect,]))
}
layer_mode<-array(0, dim=length(sectors), dimnames = list(sectors))
for (sect in sectors){
  layer_mode[sect]<- which.max(int_layer_res[sect,])
}
int_mode_layers<-layer_mode[order(layer_mode)]


#Export results (Uncomment as needed)
 # write.csv(int_fln,"Results/Int_fln.csv")
 # write.csv(int_fln_rank,"Results/Int_fln_rank.csv")
 # write.csv(int_sect_fln_pos,"Results/Int_sect_fln_pos.csv")
 # write.csv(int_sect_layer, paste0("Results/TESTInt_mode_layers_",sector,".csv"))
 # write.csv(int_layer_res, paste0("Results/TESTInt_mode_layers_",sector,".csv"))
 # write.csv(int_mode_layers, paste0("Results/TESTInt_mode_layers_",sector,".csv"))
