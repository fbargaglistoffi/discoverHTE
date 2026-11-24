test_that("process_cluster drops unused factor levels before counting clusters", {
  # Regression test: qF() keeps unused levels, and the small-sample adjustment counts
  # levels, so a cluster variable carrying unused levels used to overstate the number of
  # clusters and mis-scale every standard error.
  full <- factor(paste0("c", 1:10))
  kept <- full[1:6]                      # 4 unused levels retained
  dropped <- droplevels(kept)

  d_kept <- data.frame(cl = kept)
  d_drop <- data.frame(cl = dropped)

  pc_kept <- process_cluster(~cl, d_kept)
  pc_drop <- process_cluster(~cl, d_drop)

  expect_identical(fnlevels(pc_kept$cluster[[1L]]), 6L)
  expect_equal(pc_kept$adj, pc_drop$adj)
  expect_equal(unname(pc_kept$adj), 6 / 5)
})

test_that("unused cluster levels do not change the cluster-robust variance", {
  d <- sim_erf_data(n = 300, clsize = 10)

  d_kept <- d
  d_kept$cl <- factor(d_kept$cl, levels = c(levels(d_kept$cl), "unused1", "unused2"))

  e_drop <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, cluster = ~cl)
  e_kept <- estimate_erf(Y ~ A | X1 + X2, data = d_kept, n = 21L, cluster = ~cl)

  expect_equal(e_drop$erf(e_drop$values, vcov = TRUE)$vcov,
               e_kept$erf(e_kept$values, vcov = TRUE)$vcov)
})

test_that("more than one clustering variable is rejected", {
  d <- sim_erf_data(n = 200, clsize = 10)
  d$cl2 <- factor(rep(seq_len(4), length.out = nrow(d)))

  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, cluster = ~ cl + cl2),
               regexp = "one clustering variable", ignore.case = TRUE)
  expect_error(process_cluster(~ cl + cl2, d),
               regexp = "one clustering variable", ignore.case = TRUE)
})

test_that("cross-fitting folds keep whole clusters together", {
  cl <- factor(rep(seq_len(20), each = 5))

  for (K in c(2L, 3L, 5L)) {
    f <- make_folds(K, length(cl), cl)

    expect_length(f, length(cl))
    expect_setequal(unique(f), seq_len(K))

    spread <- tapply(f, cl, function(z) length(unique(z)))
    expect_true(all(spread == 1L))
  }
})

test_that("make_folds falls back to row-level folds when unclustered", {
  f <- make_folds(4L, 100L, NULL)

  expect_length(f, 100L)
  expect_setequal(unique(f), 1:4)
  expect_equal(as.integer(table(f)), rep(25L, 4L))
})

test_that("estimate_erf stores the clustering aligned to the analysis rows", {
  d <- sim_erf_data(n = 300, clsize = 10)

  e_none <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)
  expect_null(e_none$cluster)

  e_cl <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, cluster = ~cl, crossfit = 3)
  expect_false(is.null(e_cl$cluster))
  expect_identical(ncol(e_cl$cluster), 1L)
  expect_identical(nrow(e_cl$cluster), length(e_cl$treat))
})

test_that("singleton clusters reproduce the iid variance", {
  d <- sim_erf_data(n = 200, clsize = 10)
  d$one_each <- factor(seq_len(nrow(d)))

  e_iid <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)
  e_cl <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, cluster = ~one_each)

  v_iid <- e_iid$erf(e_iid$values, vcov = TRUE)$vcov
  v_cl <- e_cl$erf(e_cl$values, vcov = TRUE)$vcov

  # One observation per cluster: the cluster-robust sandwich collapses to the iid one.
  expect_equal(v_cl, v_iid)
})
