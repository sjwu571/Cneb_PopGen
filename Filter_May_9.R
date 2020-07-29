library(dartR)

gl <- gl.read.dart(filename="Report_DCyn19-4420_1_moreOrders_SNP_mapping_2_reorg.csv", 
                   covfilename="ind_metrics_locations_reorg.csv", probar = FALSE)

save(gl, file="gl_DCyn19_reorg.rdata")

# delete 1_77 & 1_78 due to mixed samples 

gl <- gl.drop.ind(gl, c("1_77", "1_78"))
gl <- gl.filter.repavg(gl,t = 0.9) # 59884
gl <- gl.filter.monomorphs(gl)
gl <- gl.filter.callrate(gl, method = "loc", threshold = 0.9, recalc = T) # 39896
gl <- gl.filter.secondaries(gl) #29457
gl <- gl.filter.callrate(gl, method = "ind", threshold = 0.9, recalc = T) # 29452
gl <- gl.filter.maf(gl, threshold = 0.01) # 17506


glout <- gl

if (is.null(glout@other$loc.metrics$rdepth)) {
  
  cat("Read depth calculated and added to @loc.metrics slot.\n")
  glout@other$loc.metrics$rdepth <- array(NA,nLoc(glout))
  for (i in 1:nLoc(glout)){
    called.ind <- round(nInd(glout)*glout@other$loc.metrics$CallRate[i],0)
    ref.count <- called.ind*glout@other$loc.metrics$OneRatioRef[i]
    alt.count <- called.ind*glout@other$loc.metrics$OneRatioSnp[i]
    sum.count.ref <- ref.count*glout@other$loc.metrics$AvgCountRef[i]
    sum.count.alt <- alt.count*glout@other$loc.metrics$AvgCountSnp[i]
    glout@other$loc.metrics$rdepth[i] <- round((sum.count.alt + sum.count.ref)/called.ind,1)
  }
  cat("All read in. Please check carefully the output above\n")  
  
} 

x <- glout
n0 <- nLoc(x); lower <- 5; upper <- 100

index <- (x@other$loc.metrics["rdepth"]>=lower & x@other$loc.metrics["rdepth"]<= upper)
x2 <- x[, index]
# Remove the corresponding records from the loci metadata
x2@other$loc.metrics <- x@other$loc.metrics[index,]

x2 # 15187

save(x2, file = "gl_May_9_filter.RData")
