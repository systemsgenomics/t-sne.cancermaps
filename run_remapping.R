# run_remapping: Wrapper function to run remapping algorithm
# originalData: gene expression data from original map
# map_coord: coordinates from originalData TSNE projection
# newData: gene expression data from new samples
# perplexity: "soft" number of neighbors assumed. Used in TSNE run
# theta: portion of approximation used when running Barnes-Hut algorithm

run_remapping = function(originalData, map_coord, newData, perplexity=30, theta=0.3, CORES=1){
  library(uuid)
  library(parallel)
  
  # Filter to common genes between data
  common_genes = intersect(rownames(originalData), rownames(newData))
  originalData=originalData[match(common_genes, rownames(originalData)),]
  newData=newData[match(common_genes, rownames(newData)),]
  
  # Runs remapping on given expression matrix against premapped maps.
  remapped_coord = fast_tsne_remapping(originalData,
                                       perplexity=perplexity,
                                       theta=theta,
                                       map_coord,
                                       newData,
                                       CORES)
  return(remapped_coord)
}

# fast_tsne_remapping: actual function to run remapping algorithm
fast_tsne_remapping = function(expr_samples_in_map,
                               perplexity=30,
                               theta=0.3,
                               lockedMappings,
                               expr_new_samples,
                               CORES=1){
  
  write_data = function(expr, theta, perplexity,
                        lockedMappings, file_){
    data_file = file(file_, 'wb')
    writeBin(nrow(expr), data_file, size = 4)
    writeBin(ncol(expr), data_file, size = 4)
    writeBin(theta, data_file)
    writeBin(perplexity, data_file)
    writeBin(as.vector(t(expr)), data_file)
    writeBin(as.vector(as.matrix(t(lockedMappings))), data_file)
    close(data_file)
  }
  
  read_results = function(res_file){
    conn = file(res_file, "rb")
    nrows = readBin(conn, "integer", 1, 4)
    ncols = readBin(conn, "integer", 1, 4)
    print("nrows and ncols")
    print(nrows)
    print(ncols)
    coord = readBin(conn, "double", nrows*ncols)
    close(conn)
    
    # Reshape to matrix
    dim(coord) = c(ncols, nrows)
    return(t(coord))
    
  }
  
  tmp_dir = paste0("tmp_data-",  UUIDgenerate())
  dir.create(tmp_dir, showWarnings = FALSE)
  
  # Transpose and center gene expressions
  nLocked = ncol(expr_samples_in_map)
  nNew = ncol(expr_new_samples)
  combined_samples = rbind( t(expr_samples_in_map), t(expr_new_samples))
  combined_samples.ctr = scale(combined_samples, scale=F)
  expr_samples_in_map = combined_samples.ctr[1:nLocked,]
  expr_new_samples = combined_samples.ctr[(nLocked+1):(nLocked+nNew), ,drop=F]
  
  remappedCoord = matrix(0, nrow = nrow(expr_new_samples), ncol = 2)
  
  # TODO: Parallellize, done
  # for (i in 1:nNew){
  #   tmp_expr = rbind(expr_samples_in_map,
  #                    expr_new_samples[i,])
  #   tmp_file = file.path(tmp_dir,
  #                        paste0("data", as.character(i), ".dat"))
  #   res_file = file.path(tmp_dir,
  #                        paste0("result", as.character(i), ".dat"))
  #   print(tmp_file)
  #   write_data(tmp_expr, theta, perplexity, lockedMappings,
  #              file_ = tmp_file)
  #   cmd = paste("./useCase2/bh_tsne_remapping", tmp_file, res_file, sep=" ")
  #   print(paste0("Running command ", cmd))
  #   system(cmd)
  #   resCoord = read_results(res_file)
  #   remappedCoord[i,] = resCoord[nrow(resCoord),]
  #   #unlink(tmp_file)
  #   #unlink(res_file)
  # }
  # 
  # parallelization
  remappedCoord=do.call(rbind, mclapply(seq(nNew), function(i){
    tmp_expr = rbind(expr_samples_in_map,
                     expr_new_samples[i,])
    tmp_file = file.path(tmp_dir,
                         paste0("data", as.character(i), ".dat"))
    res_file = file.path(tmp_dir,
                         paste0("result", as.character(i), ".dat"))
    print(tmp_file)
    write_data(tmp_expr, theta, perplexity, lockedMappings,
               file_ = tmp_file)
    cmd = paste("./useCase2/bh_tsne_remapping", tmp_file, res_file, sep=" ")
    print(paste0("Running command ", cmd))
    system(cmd)
    resCoord = read_results(res_file)
    unlink(tmp_file)
    unlink(res_file)
    return(resCoord[nrow(resCoord),])
  }, mc.cores=CORES))
  
  return(remappedCoord)
  
}
