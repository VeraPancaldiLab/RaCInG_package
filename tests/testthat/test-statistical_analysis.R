test_that("wilcox_group_test() reports fold change and runs all pairwise comparisons", {
  set.seed(1)
  mat <- matrix(c(rnorm(10, 5), rnorm(10, 2)), nrow = 20, ncol = 1,
                dimnames = list(paste0("P", 1:20), "feat1"))
  groups <- rep(c("A", "B", "C"), length.out = 20)

  res <- wilcox_group_test(mat, groups)
  expect_true(all(c("Comparison", "Feature", "Wilcox_statistic", "Fold_change",
                     "Log2FC", "P_value", "Adjusted_P_value") %in% colnames(res)))
  # 3 groups -> 3 pairwise comparisons, 1 feature each
  expect_equal(sort(unique(res$Comparison)), sort(c("A_vs_B", "A_vs_C", "B_vs_C")))
  expect_equal(nrow(res), 3)
})

test_that("top_features_plot() requires a Log2FC column and ranks top up/down features", {
  bad_input <- data.frame(Feature = "f1", Adjusted_P_value = 0.01, Wilcox_statistic = 5)
  expect_error(top_features_plot(bad_input), "Log2FC")

  good_input <- data.frame(
    Feature = c("f1", "f2", "f3"),
    Adjusted_P_value = c(0.01, 0.01, 0.01),
    Log2FC = c(1.5, -2, 0.5)
  )
  p <- top_features_plot(good_input, top_n = 5)
  expect_s3_class(p, "ggplot")
  # a feature can never be labeled both Up and Down
  expect_false(any(duplicated(p$data$Feature)))
})

test_that("correlate_features_with_score() aligns patients by name and supports groups", {
  set.seed(1)
  feat <- matrix(1:20, nrow = 5, ncol = 4, dimnames = list(paste0("P", 1:5), paste0("F", 1:4)))
  score <- setNames(c(5, 4, 3, 2, 1), paste0("P", 1:5))
  group <- c("g1", "g1", "g2", "g2", "g2")

  res <- correlate_features_with_score(feat, score, group = group)
  expect_true(all(c("All", "g1", "g2") %in% res$Group))
  # Each column is perfectly monotonic with the (reversed) score -> rho == -1 pooled
  pooled <- res[res$Group == "All" & res$Feature == "F1", ]
  expect_equal(pooled$Rho, -1, tolerance = 1e-8)
})
