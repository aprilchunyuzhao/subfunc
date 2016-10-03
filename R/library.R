#### adonis_repeated_measures
# vegan::adonis permutational multivariate analysis of variance using distance matrices
# permanova with repeated measures
adonis_repeated_measures <- function (uu, s) {
  set.seed(10)
  # unrestricted permutations
  a_ixn <- adonis(uu ~ study_group * study_day, s, perm=999)

  # Interaction
  f_ixn <- a_ixn$aov.tab[3, 4] # F staprint.adonist for the interaction
  fs_permuted <- replicate(999, {
    s_permuted <- within(s, {
      study_group <- shuffle_between_groups(study_group, SampleID)
      study_day <- shuffle_within_groups(study_day, SampleID)
    })
    a_permuted <- adonis(uu ~ study_group * study_day, s_permuted, permutations = 4)
    a_permuted$aov.tab[3, 4]
  })
  # WHY there is NA ?? anyway, remove the NA
  #fs_permuted <- fs_permuted[!is.na(fs_permuted)]

  p_ixn <- sum(c(f_ixn, fs_permuted) >= f_ixn) / (length(fs_permuted) + 1)
  a_ixn$aov.tab[3,6] <- p_ixn

  # First factor: study_group
  f_ixn1 <- a_ixn$aov.tab[1, 4]
  fs_permuted1 <- replicate(999, {
    s_permuted <- within(s, {
      study_group <- shuffle_between_groups(study_group, SampleID)
    })
    a_permuted1 <- adonis(uu ~ study_group * study_day, s_permuted, permutations = 4)
    a_permuted1$aov.tab[1, 4]
  })
  p_ixn1 <- sum(c(f_ixn1, fs_permuted1) >= f_ixn1) / (length(fs_permuted1) + 1)
  a_ixn$aov.tab[1,6] <- p_ixn1

  # Second factor: has to be time points (study_day)
  a_ixn_day <- adonis(uu ~ study_day, data = s, strata = s$SubjectID)
  a_ixn$aov.tab[2,6] <- a_ixn_day$aov.tab[1,6]

  return(a_ixn)
}


#### tidy_anova_repeated_measure_posthoc seris
###### sub functions yay ######
## 1way or 2way anova with repeated measures
tidy_anova_repeated_measures <- function(anov){
  temp <- data.frame(summary(anov)[[2]][[1]]) # Error: SubjectID; Error: Within
  data.frame(term=gsub(" ","",rownames(temp),perl = TRUE), temp, row.names = NULL) %>%
    filter(!grepl("Residuals", term)) %>%
    rename(p.value = Pr..F.)
}

## 1way or 2way anova
tidy_anova <- function(anov) {
  temp <- data.frame(summary(anov)[[1]])
  data.frame(term = gsub(" ","", rownames(temp), perl = TRUE), temp, row.names = NULL) %>%
    filter(!grepl("Residuals", term)) %>%
    rename(p.value = Pr..F.)
}

## posthoc test for anova with repeated measures
tidy_anova_repeated_measures_posthoc <- function (s_sub, form1, group_label, p_cutoff = 0.05) {
  # one-way / two-way anova with repeated measures
  anov <- data.frame(comparison="all", tidy_anova_repeated_measures(aov(as.formula(form1), data=s_sub)))
  combs <- combn(unlist(unique(s_sub[,group_label])), 2)
  num_tests <- dim(combs)[2]

  # pairwise posthoc test, based on anova with repeated measures
  # anov$p.value < p_cutoff <- Ask Kyle: why the p-value are different...!!!
  if (anov$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(comparison = paste(combs[,x], collapse='-'),
                                               tidy_anova_repeated_measures(aov(as.formula(form1),
                                                                                data = s_sub[is.element(s_sub[[group_label]], combs[,x]),]))))

    anov <- rbind(anov, do.call(rbind, post_hocs))
  }
  anov[,!(colnames(anov) %in% "term")]
}

## posthoc test for anova
tidy_anova_posthoc <- function (s_sub, form1, group_label, p_cutoff = 0.05) {
  # one-way / two-way anova with repeated measures
  anov <- data.frame(comparison=group_label, tidy_anova(aov(as.formula(form1), data=s_sub)))
  combs <- combn(unlist(unique(s_sub[,group_label])), 2)
  num_tests <- dim(combs)[2]

  # pairwise posthoc test, based on anova with repeated measures
  # anov$p.value < p_cutoff <- Ask Kyle: why the p-value are different...!!!
  if (anov$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(comparison = paste(combs[,x], collapse='-'),
                                               tidy_anova(aov(as.formula(form1),
                                                              data = s_sub[is.element(s_sub[[group_label]], combs[,x]),]))))

    anov <- rbind(anov, do.call(rbind, post_hocs))
  }
  anov[,!(colnames(anov) %in% "term")]
}


#### define_functions: make_pcoa_plot
make_pcoa_plot <- function(uu, s, shape_by, color_by, title) {
  uu_pcoa <- pcoa(uu) ## uu: distance matrix
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="SampleID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)

  g_uu = ggplot(uu_df, aes(x=Axis.1, y=Axis.2)) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_colour_discrete(name=sub("_", " ", color_by)) +
    ggtitle(title)

  if (is.null(shape_by) & !is.null(color_by)) {
    g_uu <- g_uu + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    g_uu <- g_uu + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else if (is.null(shape_by) & !is.null(color_by)) {
    g_uu <- g_uu + geom_point(aes(colour=factor(get(color_by))))
  } else {
    g_uu <- g_uu + geom_point()
  }
  return(g_uu)
}
