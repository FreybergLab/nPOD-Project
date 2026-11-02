# irlba_small_matrix_fix.R
#
# Lets monocle3's preprocess_cds() run on the pseudotime analysis's 5-feature
# matrix (cell counts per islet).
#
# preprocess_cds() computes its PCA with irlba, which finds a few singular
# values of a large matrix quickly. Here it is asked for 4 of the 5 components
# (num_dim = 5, capped at one fewer than the number of features), too large a
# share for irlba: it warns "You're computing too large a percentage of total
# singular values, use a standard svd instead" and then fails with
# "object 'dv' not found".
#
# This replaces irlba::irlba, for the rest of the R session, with a version that
# passes such small problems to base R's exact svd(). For this matrix the exact
# SVD gives the same components irlba would, apart from arbitrary signs, which
# change no distances, so the UMAP, clusters and pseudotime are unaffected.
# Larger problems still go to the original irlba. The analysis code itself is
# unchanged.
#
# Source after loading monocle3 and before preprocess_cds():
#   source("Pseudotime/R/irlba_small_matrix_fix.R")

local({
  if (isTRUE(attr(irlba::irlba, "small_matrix_fix"))) return(invisible())   # already applied this session
  original <- irlba::irlba
  patched <- function(A, nv = 5, nu = nv, ..., center = NULL, scale = NULL) {
    if (nv < min(dim(A)) - 1) return(original(A, nv = nv, nu = nu, ..., center = center, scale = scale))
    X <- as.matrix(A)
    if (!is.null(center)) X <- sweep(X, 2, center)
    if (!is.null(scale))  X <- sweep(X, 2, scale, "/")
    s <- svd(X, nu = nu, nv = nv)
    message(sprintf("irlba_small_matrix_fix: exact SVD for a %d x %d matrix (%d components)",
                    nrow(X), ncol(X), nv))
    list(d = s$d[seq_len(nv)], u = s$u, v = s$v, iter = 0, mprod = 0)
  }
  attr(patched, "small_matrix_fix") <- TRUE
  utils::assignInNamespace("irlba", patched, ns = "irlba")
})
