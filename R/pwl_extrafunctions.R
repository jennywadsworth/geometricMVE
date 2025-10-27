# A clone of the geometry::delaunayn function, without the recent
# update (Feb 2025) with typos

qhull.options <- function(options, output.options, supported_output.options, full=FALSE) {
  if (full) {
    if (!is.null(output.options)) {
      stop("full and output.options should not be specified together")
    }
    output.options = TRUE
    ## Enable message in 0.4.1
    ## Turn to warning in 0.4.7
    message("delaunayn: \"full\" option is deprecated; adding \"Fa\" and \"Fn\" to options.
      Please update your code to use \"output.options=TRUE\" or set \"output.options\" to a
      string containing desired QHull options.")
  }

  if (is.null(output.options)) {
    output.options <- ""
  }
  if (is.logical(output.options)) {
    if (output.options) {
      output.options <- paste(supported_output.options, collapse=" ")
    } else {
      output.options  <- ""
    }
  }
  if (!is.character(output.options)) {
    stop("output.options must be a string, logical or NULL")
  }

  ## Input sanitisation
  options <- paste(options, output.options, collapse=" ")
  return(options)
}

delaunayn <-
  function(p, options=NULL, output.options=NULL, full=FALSE) {
    tmp_stdout <- tempfile("Rf")
    tmp_stderr <- tempfile("Rf")
    on.exit(unlink(c(tmp_stdout, tmp_stderr)))

    ## Coerce the input to be matrix
    if (is.data.frame(p)) {
      p <- as.matrix(p)
    }

    ## Make sure we have real-valued input
    storage.mode(p) <- "double"

    ## We need to check for NAs in the input, as these will crash the C
    ## code.
    if (any(is.na(p))) {
      stop("The first argument should not contain any NAs")
    }

    ## Default options
    default.options <- "Qt Qc Qx"
    if (ncol(p) < 4) {
      default.options <- "Qt Qc Qz"
    }
    if (is.null(options)) {
      options <- default.options
    }

    ## Combine and check options
    options <- tryCatch(PWLExtremes::qhull.options(options, output.options, supported_output.options  <- c("Fa", "Fn"), full=full), error=function(e) {stop(e)})

    ## It is essential that delaunayn is called with either the QJ or Qt
    ## option. Otherwise it may return a non-triangulated structure, i.e
    ## one with more than dim+1 points per structure, where dim is the
    ## dimension in which the points p reside.
    if (!grepl("Qt", options) & !grepl("QJ", options)) {
      options <- paste(options, "Qt")
    }

    out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout, tmp_stderr, PACKAGE="geometry")

    ## Check for points missing from triangulation, but not in the case
    ## of a degenerate trianguation (zero rows in output)
    if (nrow(out$tri) > 0) {
      missing.points <- length(setdiff(seq(1,nrow(p)), unique(as.vector(out$tri))))
      if (missing.points > 0) {
        warning(paste0(missing.points, " points missing from triangulation.
It is possible that setting the 'options' argument of delaunayn may help.
For example:
options = \"", default.options, " Qbb\"
options = \"", default.options, " QbB\"
If these options do not work, try shifting the centre of the points
to the origin by subtracting the mean coordinates from every point."))
      }
    }

    # Remove NULL elements
    out[which(sapply(out, is.null))] <- NULL
    if (is.null(out$areas) & is.null(out$neighbours)) {
      attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
      return(out$tri)
    }
    class(out) <- "delaunayn"
    out$p <- p
    return(out)
  }



###############################################################################
################################  d=2  ########################################
###############################################################################

# Given angles and reference angles, return a list containing the angles and
# the indices of reference angles
which.adj.angles.2d = function(angles,locs){
  # w         -> a vector or 2-column matrix of angles in the d=1 unit simplex
  # locs      -> vector of length Nmesh, the reference angles in the d=1 unit simplex

  n.pars = length(locs)

  lst.val = lapply(angles,function(w){
    for(j in 1:(n.pars-1)){
      if(w>=locs[j] & w<=locs[j+1]){
        break
      }
    }
    return(list(w=w,
                idx.locs=c(j,j+1)))
  })
  return(lst.val)
}

G.vol.2d = function(gauge.pars,par.locs){
    val = 0.5*sum(sapply(c(1:(length(par.locs)-1)),
                   function(j){
                       par.locs.1.x = par.locs[j]
                       par.locs.2.x = par.locs[j+1]
                       par.locs.1.y = 1-par.locs[j]
                       par.locs.2.y = 1-par.locs[j+1]
                       return(gauge.pars[j]*gauge.pars[j+1]*
                                abs(par.locs.1.x*par.locs.2.y - par.locs.2.x*par.locs.1.y))
                   }))
    return(val)
}

###############################################################################
################################  d>2  ########################################
###############################################################################

which.adj.angles = function(angles,locs){
  # For each angle, where does it live in the partition of the simplex?

  require(geometry)

  # angles are d-dimensional

  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  if(is.null(dim(angles))){
    angles=t(as.matrix(angles))
  }

  num.cols = dim(locs)[2]

  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  tsearchn.output = tsearchn(x=locs[,-num.cols],#rbind(locs,0)[,-num.cols],
                          t=del.tri$tri,
                          xi=matrix(angles[,-num.cols],ncol=num.cols-1))
  locs.idx = tsearchn.output$idx
  w.adj.angles = lapply(c(1:nrow(angles)),function(i){
    return(list(w=angles[i,],                                # the angle of interest, w
                tri.idx=locs.idx[i],                         # the index of the triangle where w belongs
                loc.idx=del.tri$tri[locs.idx[i],],           # the indices of the enclosing vertices
                vertices=locs[del.tri$tri[locs.idx[i],],]))  # the enclosing vertices (a matrix of size d x d)
  })
  return(w.adj.angles)
}

adj.DT = function(par.locs){
  # given reference angles, return the indices of the neighboring
  # triangles in the Delaunay triangulation

  num.cols = dim(par.locs)[2]
  del.tri = PWLExtremes::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)
  tri = del.tri$tri
  num.tris = nrow(tri)
  adj.lst = list(list(NULL))
  for(k in 1:num.tris){
    tri.k = tri[k,]
    adj.idx = NULL
    for(j in 1:num.tris){
      if(sum(table(c(tri.k,tri[j,]))==2)==(num.cols-1)){
        # if this is satisfied, tirangles i and j are neighbours
        adj.idx = c(adj.idx,j)
      }
    }
    adj.lst[[k]] = sort(adj.idx)
  }
  return(adj.lst)
}

ij.couples = function(par.locs){
  # given N par locs (nodes), for each node, find the (i,j) pairings in the penalty
  require(dplyr)

  num.cols = dim(par.locs)[2]
  tri = PWLExtremes::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)$tri
  adj.list = adj.DT(par.locs)

  ij.couples.lst = list(list(NULL))  # eventually will be of length N
  for(l in 1:nrow(par.locs)){
    #i=20
    # which triangle has location "i" as a vertex?
    which.has.l = apply(tri,1,function(ll){
      l %in% ll
    })
    if(sum(which.has.l)==1){
      # only 1 face, can't compare gradient to another
      ij.couples.lst[[l]] = NA
      next
    } else if(sum(which.has.l)>1){
      # which.grad.gs = (1:length(grad.g))[which.has.i]
      idx.w.l = c(1:length(adj.list))[which.has.l]
      adj.list.l = adj.list[idx.w.l]
      names(adj.list.l) = idx.w.l
      adj.list.l = lapply(adj.list.l, function(vec){
        return(vec[vec %in% idx.w.l])
      })
      ij.couples = lapply(1:length(idx.w.l), function(idx){
        idx.w.l.idx = idx.w.l[idx]
        lapply(adj.list.l[[idx]], function(vec){
          lapply(vec,function(vec.i){
            return(sort(c(idx.w.l.idx,vec.i)))
          })
        })
      })
      ij.couples = matrix(unlist(ij.couples),ncol=2,byrow=T)
      ij.couples = unname(as.matrix(distinct(as.data.frame(ij.couples,col.names=NULL))))
      ij.couples.lst[[l]] = ij.couples
    }
  }
  return(ij.couples.lst)
}

G.vol = function(gauge.pars,par.locs){
  num.cols = dim(par.locs)[2]

  # Method 2: coordinate geometry
  del.tri = PWLExtremes::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)
  nodes = del.tri$tri
  vol = (1/factorial(num.cols))*sum(apply(nodes,1,function(loc){
    pp = gauge.pars[loc]
    ll = par.locs[loc,]
    return(abs(det(ll*pp)))
  }))
  return(vol)
}

# Given angles and reference angles, return a list containing the reference angles,
# their indices, and how many angles are in each region
n.DT.region = function(which.adj.angles.res,locs){
  # which.adj.angles.res  -> output of which.adj.angles()
  # locs                  -> the reference angles in the d=2 unit simplex

  all.regions = lapply(which.adj.angles.res, function(lst) sort(lst$loc.idx))
  idx.locs = unique(all.regions)

  res = lapply(idx.locs,function(loc){
    freq = sum(sapply(which.adj.angles.res, function(lst) all(sort(lst$loc.idx) == loc)))
    return(list(reg = locs[loc,],  # the region is defined by these vertices
                loc = loc,
                freq=freq))        # there are this many angles in the region
  })

  return(res)
}


###################################################################

# other functions

# sampling from f=e^-g

f.mcmc.g.2d<-function(niter,nburn,theta,alpha,thin,g){
  # initialize
  w<-rbeta(1,alpha,alpha)
  x<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)
  draws<-matrix(ncol=2,nrow=niter)
  it<-1
  while(it<=niter){
    w<-rbeta(1,alpha,alpha)
    xcan<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)

    accn<-g(xy=x,par=theta)*dbeta(x[1]/(x[1]+x[2]),alpha,alpha)*(xcan[1]+xcan[2])^2
    accd<-g(xy=xcan,par=theta)*dbeta(xcan[1]/(xcan[1]+xcan[2]),alpha,alpha)*(x[1]+x[2])^2

    if(runif(1)<accn/accd){x<-xcan}

    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

#########################################################


f.mcmc.g.3d<-function(niter,nburn,alpha=rep(1,3),thin,g){
  # initialize
  w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
  x<-rexp(1)*w/g(w)
  draws<-matrix(ncol=3,nrow=niter)
  it<-1
  while(it<=niter){
    w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
    xcan<-rexp(1)*w/g(w)

    accn<-g(x)*LaplacesDemon::ddirichlet(x/sum(x),alpha)*(sum(xcan))^3
    accd<-g(xcan)*LaplacesDemon::ddirichlet(xcan/sum(xcan),alpha)*(sum(x))^3

    if(runif(1)<accn/accd){x<-xcan}

    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

###########################################################

