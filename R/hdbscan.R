#'  Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x  matrix of shape (n_samples, n_features)
#' @param algorithm specify which algorithm to use
#' @param alpha A distance scaling parameter as used in robust single linkage.
#' @param approx_min_span_tree Whether to accept an only approximate minimum spanning tree.
#' @param gen_min_span_tree Whether to generate the minimum spanning tree with regard to mutual reachability distance for later analysis.
#' @param leaf_size If using a space tree algorithm (kdtree, or balltree) the number of points ina leaf node of the tree.
#' @param metric The metric to use when calculating distance between instances in a feature array.
#' @param min_cluster_size The minimum size of clusters
#' @param min_samples The number of samples in a neighbourhood for a point to be considered a core point.
#' @param cluster_selection_epsilon A distance threshold. Clusters below this value will be merged.
#' @param cluster_selection_method The method used to select clusters from the condensed tree.
#' @param nThreads number of parallel threads
#' @param prediction_data  not sure what this is for. Will update later.
#'
#' @return A list object (length=5) with the cluster labels for each point (labels [num]), the strength of a sample's membership to its assigned cluster (probabilities [num]), the stability of the cluster (cluster_persistance [num]),  list of exemplar points for clusters, and outlier scores for clusterd points (outlier_scores [num]))
#' @export
#'
#' @import reticulate
#'
#' @examples
#'\dontrun{
#' # load Seurat object
#' data('pbmc_small', package = 'SeuratObject')
#'
#' # tsne embedding matrix
#' emb <- Seurat::Embeddings(object = pbmc_small, 'tsne')
#'
#' # HDBSCAN.matrix
#' res.matrix <- dbsinglecell::HDBSCAN(emb, min_cluster_size = 5)
#' res.matrix
#'
#' # HDBSCAN.Seurat
#' res.seurat <- dbsinglecell::HDBSCAN(pbmc_small,reduction = 'tsne', min_cluster_size = 5)
#' slot(res.seurat, 'misc')$hdbscan
#'
#'}
#'
HDBSCAN <- function(object, ...){
  UseMethod(generic = 'HDBSCAN', object = object)
}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x  matrix of shape (n_samples, n_features)
#' @param algorithm specify which algorithm to use
#' @param alpha A distance scaling parameter as used in robust single linkage.
#' @param approx_min_span_tree Whether to accept an only approximate minimum spanning tree.
#' @param gen_min_span_tree Whether to generate the minimum spanning tree with regard to mutual reachability distance for later analysis.
#' @param leaf_size If using a space tree algorithm (kdtree, or balltree) the number of points ina leaf node of the tree.
#' @param metric The metric to use when calculating distance between instances in a feature array.
#' @param min_cluster_size The minimum size of clusters
#' @param min_samples The number of samples in a neighbourhood for a point to be considered a core point.
#' @param cluster_selection_epsilon A distance threshold. Clusters below this value will be merged.
#' @param cluster_selection_method The method used to select clusters from the condensed tree.
#' @param nThreads number of parallel threads
#' @param prediction_data  not sure what this is for. Will update later.
#'
#' @return A list object (length=5) with the cluster labels for each point (labels [num]), the strength of a sample's membership to its assigned cluster (probabilities [num]), the stability of the cluster (cluster_persistance [num]),  list of exemplar points for clusters, and outlier scores for clusterd points (outlier_scores [num]))
#' @export
#'
#' @examples
#'\dontrun{
#' # load Seurat object
#' data('pbmc_small', package = 'SeuratObject')
#'
#' # tsne embedding matrix
#' emb <- Seurat::Embeddings(object = pbmc_small, 'tsne')
#'
#' # HDBSCAN.matrix
#' res.matrix <- dbsinglecell::HDBSCAN(emb, min_cluster_size = 5)
#' res.matrix
#'
#' # HDBSCAN.Seurat
#' res.seurat <- dbsinglecell::HDBSCAN(pbmc_small,reduction = 'tsne', min_cluster_size = 5)
#' slot(res.seurat, 'misc')$hdbscan
#'
#'}
#'
HDBSCAN.default <- function(x,
                            algorithm='best',
                            alpha=1.0,
                            approx_min_span_tree = TRUE,
                            gen_min_span_tree=FALSE,
                            leaf_size=30,
                            metric='euclidean',
                            prediction_data=TRUE,
                            min_cluster_size =1000L,
                            min_samples = 1L,
                            cluster_selection_epsilon = 0.1,
                            cluster_selection_method = 'leaf',
                            nThreads = parallel::detectCores()-1,
                            return_full = FALSE
){

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm = algorithm,
                               allow_single_cluster = FALSE,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )



  reticulate::py_main_thread_func(clusterer$fit(x))


  result <- list(
    labels = clusterer$labels_,
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_)



  if(return_full==TRUE){
    return(clusterer)
  } else {
    return(result)
  }

}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param x  matrix of shape (n_samples, n_features)
#' @param algorithm specify which algorithm to use
#' @param alpha A distance scaling parameter as used in robust single linkage.
#' @param approx_min_span_tree Whether to accept an only approximate minimum spanning tree.
#' @param gen_min_span_tree Whether to generate the minimum spanning tree with regard to mutual reachability distance for later analysis.
#' @param leaf_size If using a space tree algorithm (kdtree, or balltree) the number of points ina leaf node of the tree.
#' @param metric The metric to use when calculating distance between instances in a feature array.
#' @param min_cluster_size The minimum size of clusters
#' @param min_samples The number of samples in a neighbourhood for a point to be considered a core point.
#' @param cluster_selection_epsilon A distance threshold. Clusters below this value will be merged.
#' @param cluster_selection_method The method used to select clusters from the condensed tree.
#' @param nThreads number of parallel threads
#' @param prediction_data  not sure what this is for. Will update later.
#'
#' @return A list object (length=5) with the cluster labels for each point (labels [num]), the strength of a sample's membership to its assigned cluster (probabilities [num]), the stability of the cluster (cluster_persistance [num]),  list of exemplar points for clusters, and outlier scores for clusterd points (outlier_scores [num]))
#' @export
#'
#' @examples
#'\dontrun{
#' # load Seurat object
#' data('pbmc_small', package = 'SeuratObject')
#'
#' # tsne embedding matrix
#' emb <- Seurat::Embeddings(object = pbmc_small, 'tsne')
#'
#' # HDBSCAN.matrix
#' res.matrix <- dbsinglecell::HDBSCAN(emb, min_cluster_size = 5)
#' res.matrix
#'
#' # HDBSCAN.Seurat
#' res.seurat <- dbsinglecell::HDBSCAN(pbmc_small,reduction = 'tsne', min_cluster_size = 5)
#' slot(res.seurat, 'misc')$hdbscan
#'
#'}
#'
HDBSCAN.matrix <- function(x,
                           algorithm='best',
                           alpha=1.0,
                           approx_min_span_tree = TRUE,
                           gen_min_span_tree=FALSE,
                           leaf_size=30,
                           metric='euclidean',
                           prediction_data=TRUE,
                           min_cluster_size =1000L,
                           min_samples = 1L,
                           cluster_selection_epsilon = 0.1,
                           cluster_selection_method = 'leaf',
                           nThreads = parallel::detectCores()-1,
                           return_full = FALSE
){

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm = algorithm,
                               allow_single_cluster = FALSE,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )



  reticulate::py_main_thread_func(clusterer$fit(x))


  result <- list(
    labels = clusterer$labels_,
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_)

  if(return_full==TRUE){
    return(clusterer)
  } else {
    return(result)
  }

  reticulate::py_de
  reticulate::py_run_string('del clusterer')
  reticulate::py_gc <- import("gc")
  reticulate::py_gc$collect()
}



#'  Hierarchical Density-Based Spatial Clustering of Applications with Noise
#'
#' @param object Seurat object
#' @param reduction name of embedding to use (default='umap')
#' @param dims number of dimensions to use from the reduction embedding. If not set (default=NULL), uses all available dimensions.
#' @param algorithm specify which algorithm to use
#' @param alpha A distance scaling parameter as used in robust single linkage.
#' @param approx_min_span_tree Whether to accept an only approximate minimum spanning tree.
#' @param gen_min_span_tree Whether to generate the minimum spanning tree with regard to mutual reachability distance for later analysis.
#' @param leaf_size If using a space tree algorithm (kdtree, or balltree) the number of points ina leaf node of the tree.
#' @param metric The metric to use when calculating distance between instances in a feature array.
#' @param min_cluster_size The minimum size of clusters
#' @param min_samples The number of samples in a neighbourhood for a point to be considered a core point.
#' @param cluster_selection_epsilon A distance threshold. Clusters below this value will be merged.
#' @param cluster_selection_method The method used to select clusters from the condensed tree.
#' @param nThreads number of parallel threads
#' @param prediction_data  not sure what this is for. Will update later.
#' @param return_seurat  logical to return the result within the orignal object or as the raw HDBSCAN result
#'
#' @return returns the result from HDBSCAN.matrix or the original Seurat object with the result from HDBSCAN.matrix stored in the misc slot.
#' @export
#'
#' @examples
#'\dontrun{
#' # load Seurat object
#' data('pbmc_small', package = 'SeuratObject')
#'
#' # tsne embedding matrix
#' emb <- Seurat::Embeddings(object = pbmc_small, 'tsne')
#'
#' # HDBSCAN.matrix
#' res.matrix <- dbsinglecell::HDBSCAN(emb, min_cluster_size = 5)
#' res.matrix
#'
#' # HDBSCAN.Seurat
#' res.seurat <- dbsinglecell::HDBSCAN(pbmc_small,reduction = 'tsne', min_cluster_size = 5)
#' slot(res.seurat, 'misc')$hdbscan
#'
#'}
#'
HDBSCAN.Seurat <- function(object,
                           reduction = 'umap',
                           dims = NULL,
                           algorithm='best',
                           alpha=1.0,
                           prediction_data = TRUE,
                           approx_min_span_tree = TRUE,
                           gen_min_span_tree=FALSE,
                           leaf_size=30,
                           metric='euclidean',
                           min_cluster_size =1000L,
                           min_samples = 1L,
                           cluster_selection_epsilon = 0.1,
                           cluster_selection_method = 'leaf',
                           nThreads = parallel::detectCores()-1,
                           return_seurat = TRUE,
                           return_full = FALSE
){

  if(is.null(dims)){
    x <- Seurat::Embeddings(object, reduction = reduction)
  } else {
    x <- Seurat::Embeddings(object, reduction = reduction)[,dims]
  }

  hdbscan <- reticulate::import('hdbscan', delay_load = TRUE)



  clusterer <- hdbscan$HDBSCAN(algorithm=algorithm,
                               alpha = alpha,
                               prediction_data = prediction_data,
                               approx_min_span_tree = approx_min_span_tree,
                               gen_min_span_tree = gen_min_span_tree,
                               leaf_size = leaf_size,
                               core_dist_n_jobs = nThreads,
                               metric = metric,
                               min_cluster_size = as.integer(min_cluster_size),
                               min_samples = as.integer(min_samples),
                               cluster_selection_epsilon =  cluster_selection_epsilon,
                               cluster_selection_method = cluster_selection_method
  )
  reticulate::py_main_thread_func(clusterer$fit(x))



  result <- list(
    labels = factor(clusterer$labels_),
    probabilities = clusterer$probabilities_,
    cluster_persistance = clusterer$cluster_persistence_,
    exemplars = clusterer$exemplars_,
    outlier_scores = clusterer$outlier_scores_)

  if(return_seurat){
    if(return_full==TRUE){
      object@misc$hdbscan <- clusterer
    } else {
      object@misc$hdbscan <- result
    }

    object$hdbscan_clusters <- result$labels
    object$outlier_scores <- result$outlier_scores
    return(object)
  } else {

    if(return_full==TRUE){
      return(clusterer)
    } else {
      return(result)
    }



  }
  reticulate::py_de
  reticulate::py_run_string('del clusterer')
  reticulate::py_gc <- import("gc")
  reticulate::py_gc$collect()
}
