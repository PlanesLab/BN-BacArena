simEnvBN <- function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, diffusion=TRUE, diff_par=FALSE, cl_size=2, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=TRUE, verbose=TRUE, bacCoeff=NULL, nutCoeff=NULL){
  
  # Check coefficient matrix
  if(length(bacCoeff) != 0 && (nrow(bacCoeff) != length(arena@specs) || ncol(bacCoeff) != length(arena@specs))){
    stop('Cell coefficient matrix should have the same number of rows and columns as the number of cell types')
  }
  if(length(nutCoeff) != 0 && (nrow(nutCoeff) != length(arena@specs) || ncol(nutCoeff) != sum(is.element(colnames(nutCoeff), names(arena@media))))){
    stop('Nutrient coefficient matrix should have the same number of rows as the number of cell types and columns as the number of metabolites in the media')
  }
  
  if(length(object@media)==0) stop("No media present in Arena!")
  switch(class(object),
         "Arena"={arena <- object; evaluation <- Eval(arena)},
         "Eval"={arena <- getArena(object); evaluation <- object},
         stop("Please supply an object of class Arena or Eval."))
  if( sec_obj=="none" & any(sapply(object@specs, function(s){s@limit_growth})))
    warning("If growth is limited by maximum weight for an organism (max_weight=TRUE) it is recommended to use minimize total flux (sec_obj='mtf').")
  
  if(is.null(lrw)){lrw=estimate_lrw(arena@n,arena@m)}
  for(i in names(arena@specs)){
    phensel <- arena@phenotypes[which(names(arena@phenotypes)==i)]
    if(length(phensel)==0){
      test = getPhenotype(arena@specs[[i]], cutoff=pcut, fbasol=arena@specs[[i]]@fbasol)
      pvec = rep(0,length(arena@mediac))
      names(pvec) = arena@mediac
      pvec[names(test)] = test
      pvec <- paste(pvec,collapse='')
      names(pvec) = i
      arena@phenotypes <- c(arena@phenotypes,pvec)
    }
  }
  if(class(object)!="Eval"){addEval(evaluation, arena)}
  arena@sublb <- getSublb(arena)
  arena@exchanges <- data.frame() # remember exchanges
  diff_t=0
  if(arena@stir){ #create all possible positions on arena
    allxy = expand.grid(1:arena@n,1:arena@m)
    colnames(allxy) = c("x","y")
  }
  if(length(arena@specs) > 0) biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
  for(i in 1:time){
    init_t <- proc.time()[3]
    sublb <- arena@sublb
    if(nrow(arena@orgdat) > 1){
      new_ind = sample(1:nrow(arena@orgdat),nrow(arena@orgdat)) #shuffle through all bacteria to increase randomness
      arena@orgdat = arena@orgdat[new_ind,]
      sublb = sublb[new_ind,] #apply shuffeling also to sublb to ensure same index as orgdat
    }
    #if(verbose) cat("\niteration-start:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
    #if(i==1){ 
    #org_stat <- sapply(seq_along(arena@specs), function(x){dim(arena@orgdat[which(arena@orgdat$type==x),])[1]})
    #if(length(arena@specs) > 0){
    #  old_biomass<-biomass_stat; biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
    #  org_stat <- cbind(org_stat, biomass_stat); rownames(org_stat) <- names(arena@specs); colnames(org_stat) <- c("count", "biomass")
    #  if(verbose) print(as.data.frame(org_stat))}}
    #if(verbose & i!=1) print(as.data.frame(org_stat)[,1:2])
    arena@mflux <- lapply(arena@mflux, function(x){numeric(length(x))}) # empty mflux pool
    arena@shadow <-lapply(arena@shadow, function(x){numeric(length(x))}) # empty shadow pool
    if(nrow(arena@orgdat) > 0){ # if there are organisms left
      org.count <- nrow(arena@orgdat)
      for(j in 1:org.count){ # for each organism in arena
        if(verbose) cat("\rOrganims",j,"/",org.count)
        org <- arena@specs[[arena@orgdat[j,'type']]]
        bacnum = round((arena@scale/(org@cellarea*10^(-8)))) #calculate the number of bacteria individuals per gridcell
        switch(class(org),
               "Bac"= {arena = simBacBN(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow, bacCoeff = bacCoeff, nutCoeff = nutCoeff)}, #the sublb matrix will be modified within this function
               "Human"= {arena = simHum(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow)}, #the sublb matrix will be modified within this function
               stop("Simulation function for Organism object not defined yet."))
      }
      test <- is.na(arena@orgdat$biomass)
      if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
      rm("test")
    }
    if(verbose) cat("\r")
    if(diffusion && !arena@stir){
      if(diff_par){
        diff_t <- system.time(arena <- diffuse_par(arena, cluster_size=cl_size, lrw=lrw, sublb=sublb) )[3]
      }else diff_t <- system.time(arena <- diffuse(arena, lrw=lrw, sublb=sublb, verbose=verbose) )[3]
    }
    if(!diffusion){
      if(nrow(sublb)>0){
        for(k in 1:length(arena@media)){
          for(l in 1:nrow(sublb)){
            arena@media[[k]]@diffmat[sublb[l,"y"],sublb[l,"x"]] = sublb[l,k+2] # first two columns are coordinates
          }
        }
      }
      arena@sublb <- getSublb(arena)
    }
    if(arena@stir){ #stir environment -> random movement of bacteria + perfect diffusion
      sublb_tmp = arena@orgdat[,c("x","y")]
      for(sub in names(arena@media)){ #go through each metabolite in medium
        sumc = sum(arena@media[[sub]]@diffmat) #sum of all concentrations
        meanc = sumc/(arena@n*arena@m) #mean per grid cell
        conc = ((sumc-(meanc*nrow(sublb)))+sum(sublb[,sub]))/(arena@n*arena@m) #remove concentrations where bacteria are sitting + add the current concentration in their position
        arena@media[[sub]]@diffmat = Matrix::Matrix(conc,nrow=arena@m,ncol=arena@n,sparse=TRUE) #create matrix with homogen concentration
        sublb_tmp[,sub] = conc #create a new sublb matrix
      }
      newpos = allxy[sample(1:nrow(allxy),nrow(arena@orgdat)),]
      arena@orgdat[,c('x','y')] = newpos
      sublb_tmp[,c('x','y')] = newpos
      arena@sublb = as.matrix(sublb_tmp)
    }
    idx.rm <- which(arena@removeM > 0, arr.ind=TRUE) # remove organisms given removal matrix
    if( nrow(idx.rm) > 0 ){
      idx.rm.str <- apply(idx.rm, 1, function(r){paste0(r,collapse=",")})
      idx.orgdat.str <- apply(arena@orgdat[,c('x','y')], 1, function(r){paste0(r,collapse=",")})
      rm.rows <- which(!is.na(match(idx.orgdat.str, idx.rm.str)))
      if( length(rm.rows)> 0 ){
        arena@orgdat <- arena@orgdat[-rm.rows,]
        if(verbose) cat("removed", length(rm.rows), "organisms\n")
      }
    }
    
    addEval(evaluation, arena)
    if(reduce && i<time){evaluation = redEval(evaluation)}
    if(nrow(arena@orgdat)==0 && !continue){
      if(verbose) print("All organisms died!")
      break
    }
    step_t <- proc.time()[3] - init_t
    if(verbose) cat("\niteration:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
    if(verbose) cat("\r")
    org_stat <- sapply(seq_along(arena@specs), function(x){dim(arena@orgdat[which(arena@orgdat$type==x),])[1]})
    if(length(arena@specs) > 0){
      old_biomass<-biomass_stat; biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
      org_stat <- cbind(org_stat, biomass_stat, 100*(biomass_stat-old_biomass)/old_biomass); rownames(org_stat) <- names(arena@specs); colnames(org_stat) <- c("count", "biomass", "%")
      if(verbose) print(as.data.frame(org_stat))}
    if(verbose) cat("\r")
    if(verbose) cat("\ttime total: ", round(step_t,3), "\tdiffusion: ", round(diff_t,3), " (", 100*round(diff_t/step_t,3),"%)\n\n" )
    if(verbose) cat("--------------------------------------------------------------------\n")
  }
  return(evaluation)
}

simBacBN <- function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE, bacCoeff=NULL, nutCoeff=NULL){
  predator_found <- FALSE
  if( object@predator != ""){
    pos  <- arena@orgdat[,c('x','y')]
    nb <- emptyHood(object, pos, arena@n, arena@m, pos[j,1], pos[j,2], arena@occupyM, inverse=T)  
    unlist(strsplit(nb,'_'))
    nb_types <- sapply(strsplit(nb,'_'), function(coord){
      arena@orgdat[arena@orgdat$x==coord[1] & arena@orgdat$y==coord[2], "type"]
    })
    nb_names <- unique(names(arena@specs)[unlist(nb_types)])
    if( object@predator %in% nb_names){
      predator_found <- TRUE
    }
  }
  if( predator_found ){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,], factor=arena@orgdat[j,"biomass"])))
    dead <- TRUE
    arena@orgdat[j, "biomass"] <- NA
    #print("died because of predator")
  }else{
    const <- constrain(object, object@medium, lb=-sublb[j,object@medium]/bacnum, #scale to population size
                       dryweight=arena@orgdat[j,"biomass"], tstep=arena@tstep, scale=arena@scale, j)
    lobnd <- const[[1]]; upbnd <- const[[2]]
    optimization <- optimizeLP(object, lb=lobnd, ub=upbnd, j=j, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)
    fbasol <- optimization[[1]]
    
    ############################## BIOMASS UPDATE ##############################
    if (length(bacCoeff)!=0 && length(nutCoeff)!=0){
      cellAbundance <- (unname(table(arena@orgdat[,'type'])) * 100) / nrow(arena@orgdat)
      nutAbundance <- sublb[j, match(colnames(nutCoeff), colnames(sublb))] / (10^12 *0.01 * arena@scale)
      factor <- (sum(bacCoeff[arena@orgdat[j,'type'],] * cellAbundance) + sum(nutCoeff[arena@orgdat[j,'type'],] * nutAbundance) + cellAbundance[arena@orgdat[j,'type']]) / cellAbundance[arena@orgdat[j,'type']]
      arena@orgdat[j,'biomass'] <- factor * arena@orgdat[j,'biomass']
    } else if (length(bacCoeff)!=0 && length(nutCoeff)==0){
      cellAbundance <- (unname(table(arena@orgdat[,'type'])) * 100) / nrow(arena@orgdat)
      factor <- (sum(bacCoeff[arena@orgdat[j,'type'],] * cellAbundance) + cellAbundance[arena@orgdat[j,'type']]) / cellAbundance[arena@orgdat[j,'type']]
      arena@orgdat[j,'biomass'] <- factor * arena@orgdat[j,'biomass']
    } else if (length(bacCoeff)==0 && length(nutCoeff)!=0){
      nutAbundance <- sublb[j, match(colnames(nutCoeff), colnames(sublb))] / (10^12 *0.01 * arena@scale)
      factor <- (sum(nutCoeff[arena@orgdat[j,'type'],] * nutAbundance) + cellAbundance[arena@orgdat[j,'type']]) / cellAbundance[arena@orgdat[j,'type']]
      arena@orgdat[j,'biomass'] <- factor * arena@orgdat[j,'biomass']
    }
    ############################################################################
    
    eval.parent(substitute(sublb[j,] <- consume(object, sublb[j,], bacnum=bacnum, fbasol=fbasol, cutoff) )) #scale consumption to the number of cells?
    
    dead <- growth(object, arena, j, arena@occupyM, fbasol=fbasol, tstep=arena@tstep)
    arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, org=object, fbasol=fbasol, cutoff=pcut))
    
    type <- object@type
    arena@mflux[[type]]  <- arena@mflux[[type]] + fbasol$fluxes # remember active fluxes
    arena@shadow[[type]] <- arena@shadow[[type]]+ optimization[[2]]
    idx <- match(arena@mediac, names(fbasol$fluxes))
    exchanges <- data.frame(type,t(fbasol$fluxes[idx]))
    colnames(exchanges) <- c("species", unname(arena@mediac))
    arena@exchanges <- rbind(arena@exchanges, exchanges) # remember exchanges
  }
  
  
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,])))
  }
  if(!dead && !arena@stir && object@speed != 0){
    if(object@chem[1] == ''){
      pos <- arena@orgdat[,c('x','y')]
      mov_pos <- move(object, pos, arena@n, arena@m, j, arena@occupyM)
      arena@orgdat[,c('x','y')] <- mov_pos
    }else{
      for (v in seq_along(object@chem)){
        chemo <- object@chem[[v]]
        chemo_pos <- chemotaxis(object, arena, j, chemo, arena@occupyM)
        if(!is.null(chemo_pos)){arena@orgdat[j,c('x','y')] <- chemo_pos}
      }
    }
  }
  return(arena)
}

estimate_lrw <- function(grid_n, grid_m){
  x=c(10*10, 25*25, 51*51, 61*61, 71*71, 81*81, 91*91, 101*101)
  y=c(3901, 29911, 160000, 230000, 330000, 430000, 580000, 710000)
  lm <- lm(y~x)
  #summary(lm)
  #plot(x,y)
  #abline(coef(lm))
  #abline(coef=c(0, lm$coefficients[2]))
  lrw <- as.numeric(lm$coefficients[2]*grid_n*grid_m + grid_n*grid_m*100)
  #lrw <- ((grid_n*grid_m)*18.5 + 20)*10 -> alternative function
  return(lrw)
}

loadMAT <- function (file){
  
  print(system.time(data <- R.matlab::readMat(file)))
  dat.mat <- data[[1]]
  mod.var <- dimnames(dat.mat)[[1]]
  if ("modelID" %in% mod.var) {
    mod.id <- as.character(dat.mat[[which(mod.var == "modelID")]])
  }
  else {
    mod.id <- names(data)[1]
  }
  if ("modelName" %in% mod.var) {
    mod.name <- as.character(dat.mat[[which(mod.var == "modelName")]])
  }
  else {
    mod.name <- mod.id
  }
  if ("description" %in% mod.var) {
    mod.desc <- as.character(dat.mat[[which(mod.var == "description")]])
  }
  else {
    mod.desc <- mod.id
  }
  mod.S <- Matrix(dat.mat[[which(mod.var == "S")]], sparse = T)
  mod.react_id <- unlist(dat.mat[[which(mod.var == "rxns")]])
  mod.react_name <- unlist(dat.mat[[which(mod.var == "rxnNames")]])
  if ("rev" %in% mod.var) {
    mod.react_rev <- as.vector(dat.mat[[which(mod.var == 
                                                "rev")]]) == TRUE
  }
  else {
    mod.react_rev <- as.vector(dat.mat[[which(mod.var == 
                                                "lb")]]) < 0
  }
  mod.met_id <- unlist(dat.mat[[which(mod.var == "mets")]])
  mod.met_name <- unlist(dat.mat[[which(mod.var == "metNames")]])
  
  mod.lb <- as.vector(dat.mat[[which(mod.var == "lb")]])
  mod.ub <- as.vector(dat.mat[[which(mod.var == "ub")]])
  met_comp <- stringr::str_extract_all(mod.met_id, "(?<=\\[)[a-z](?=\\])")
  if (all(sapply(met_comp, length) == 0)) {
    met_comp <- stringr::str_extract_all(mod.met_id, "(?<=_)[a-z][0-9]?(?=$)")
  }
  mod.mod_compart <- unique(unlist(met_comp))
  mod.met_comp <- match(met_comp, mod.mod_compart)
  sub <- sapply(dat.mat[[which(mod.var == "subSystems")]], 
                unlist)
  sub.unique <- unique(sub)
  mod.subSys <- Matrix(FALSE, nrow = length(sub), ncol = length(sub.unique), 
                       sparse = T)
  for (i in 1:length(sub)) {
    j <- match(sub[i], sub.unique)
    mod.subSys[i, j] <- TRUE
  }
  colnames(mod.subSys) <- sub.unique
  model <- sybil::modelorg(id = mod.id, name = mod.name)
  model@mod_desc <- mod.desc
  model@S <- mod.S
  model@lowbnd <- mod.lb
  model@uppbnd <- mod.ub
  model@met_id <- mod.met_id
  model@met_name <- mod.met_name
  model@met_num <- length(mod.met_id)
  model@react_id <- mod.react_id
  model@react_name <- mod.react_name
  model@react_num <- length(mod.react_id)
  model@react_rev <- mod.react_rev
  model@mod_compart <- mod.mod_compart
  model@met_comp <- mod.met_comp
  model@subSys <- mod.subSys
  obj.idx <- which(dat.mat[[which(mod.var == "c")]] != 0)
  if (length(obj.idx) > 0) {
    model <- sybil::changeObjFunc(model, react = obj.idx)
    print(sybil::optimizeProb(model))
  }
  
  return(model)
}