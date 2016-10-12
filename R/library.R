#### adonis_repeated_measures ####
adonis_repeated_measures <- function (uu, s) {
  set.seed(10)
  # unrestricted permutations
  a_ixn <- adonis(uu ~ study_group * study_day, s, perm=999)

  # Interaction
  f_ixn <- a_ixn$aov.tab[3, 4] # F staprint.adonist for the interaction
  fs_permuted <- replicate(999, {
    s_permuted <- within(s, {
      study_group <- shuffle_between_groups(study_group, SubjectID)
      study_day <- shuffle_within_groups(study_day, SubjectID)
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
      study_group <- shuffle_between_groups(study_group, SubjectID)
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

#### test for difference in group centroid position ####
#### tidy permanova repeated measures posthoc series ####
dataframe_permanova <- function(a_ixn) {
  ## vegan::adonis permutational multivariate analysis of variance using distance matrices
  ## one factor / two factor
  data.frame(Term = rownames(a_ixn$aov.tab), a_ixn$aov.tab, row.names = NULL) %>%
    rename(p.value = Pr..F.)
}

tidy_permanova_one_way <- function(dist_mat, s_toTest, grp1, perm=99) {
  dist_toTest <- dist_subset(dist_mat, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", grp1)
  set.seed(11)
  dataframe_permanova(adonis(as.formula(form1), data=s_toTest, permutations=perm))
}

tidy_permanova_one_way_posthoc <- function(dist_mat, s_toTest, grp1, perm=99, p_cutoff=0.05){
  set.seed(12)
  ## one way permanova post hoc
  a_ixn <- data.frame(comparison=grp1, tidy_permanova_one_way(dist_mat, s_toTest, grp1, perm=perm)[1,])
  combs <- combn(unique(s_toTest[[grp1]]), 2)
  num_tests <- dim(combs)[2]

  # do post hoc tests
  if (a_ixn$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(comparison = paste(combs[,x], collapse=' - '),
                                               tidy_permanova_one_way(dist_mat,
                                                                      s_toTest[is.element(s_toTest[[grp1]], combs[,x]),], grp1,perm=perm)[1,]) )
    a_ixn <- rbind(a_ixn, do.call(rbind, post_hocs))
  }
  a_ixn[,!(colnames(a_ixn) %in% "Term")]
}

tidy_permanova_one_way_strata <- function(dist_mat, s_toTest, grp1, reps, perm=99) {
  dist_toTest <- dist_subset(dist_mat, s_toTest$SampleID)
  form1 <- paste("dist_toTest", "~", grp1)
  set.seed(20)
  a_ixn_day <- adonis(as.formula(form1), data = s_toTest, strata = s_toTest[[reps]], permutations=perm)
  dataframe_permanova(a_ixn_day)
}

tidy_permanova_one_way_strata_posthoc <- function(dist_mat, s_toTest, grp1, reps, perm=99, p_cutoff=0.05){
  ## one way permanova post hoc
  set.seed(30)
  a_ixn <- data.frame(comparison=grp1, tidy_permanova_one_way_strata(dist_mat, s_toTest, grp1, reps, perm=perm)[1,])
  combs <- combn(unique(s_toTest[[grp1]]), 2)
  num_tests <- dim(combs)[2]

  # do post hoc tests
  if (a_ixn$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(comparison = paste(combs[,x], collapse=' - '),
                                               tidy_permanova_one_way_strata(dist_mat,
                                                                             s_toTest[is.element(s_toTest[[grp1]], combs[,x]),], grp1, reps, perm=perm)[1,]) )
    a_ixn <- rbind(a_ixn, do.call(rbind, post_hocs))
  }
  a_ixn[,!(colnames(a_ixn) %in% "Term")]
}

tidy_permanova_repeated_measures <- function (dist_mat, s_toTest, grp1, grp2=NULL, reps, perm=99) {
  ## unrestricted permutations adnois
  set.seed(100)
  dist_toTest <- dist_subset(dist_mat, s_toTest$SampleID)

  if(is.null(grp2)){
    form1 <- paste("dist_toTest", "~", grp1)
    title = "one way permanova with repeated measures result"
  } else {
    form1 <- paste("dist_toTest", "~", grp1, "*", grp2)
    title = "two way permanova with repeated measures result"
  }

  a_ixn <- adonis(as.formula(form1), data=s_toTest, permutations = perm)

  #### update the p value for the first factor grp1 ####
  f_ixn1 <- a_ixn$aov.tab[1, 4]
  fs_permuted1 <- replicate(perm, {
    s_permuted <- s_toTest
    s_permuted[,grp1] <- shuffle_between_groups(s_permuted[[grp1]], s_permuted[[reps]])
    a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
    a_permuted$aov.tab[1, 4]
  })
  p_ixn1 <- sum(c(f_ixn1, fs_permuted1) >= f_ixn1) / (length(fs_permuted1) + 1)
  ## we want to keep the observed statistics with the perm ones => otherwise the estimate is baised
  a_ixn$aov.tab[1,6] <- p_ixn1

  if(!is.null(grp2)){
    #### upadte the p value for interaction term ####
    ## shuffle_between_groups: same SubjectID didn't get different study groups
    ## shuffle_within_groups: same SubjectID didn't get same days
    f_ixn <- a_ixn$aov.tab[3, 4]
    fs_permuted <- replicate(perm, {
      s_permuted <- s_toTest
      s_permuted[,grp1] <- shuffle_between_groups(s_permuted[[grp1]], s_permuted[[reps]])
      s_permuted[,grp2] <- shuffle_within_groups(s_permuted[[grp2]],s_permuted[[reps]])
      a_permuted <- adonis(as.formula(form1), s_permuted, permutations = 4)
      a_permuted$aov.tab[3, 4]
    })
    ## update the p value calculated from permutations
    p_ixn <- sum(c(f_ixn, fs_permuted) >= f_ixn) / (length(fs_permuted) + 1)
    a_ixn$aov.tab[3,6] <- p_ixn

    #### update the p value for the second factor: has to be time points (grp2) ####
    form2 <- paste("dist_toTest", "~", grp2)
    a_ixn_day <- adonis(as.formula(form2), data = s_toTest, strata =s_toTest[[reps]])
    a_ixn$aov.tab[2,6] <- a_ixn_day$aov.tab[1,6]
  }

  #### show the result
  dataframe_permanova(a_ixn)
}

tidy_permanova_one_way_repeated_measures_posthoc <- function(dist_mat, s_toTest, grp1, reps, perm=99, p_cutoff=0.05) {
  ## one way permanova post hoc with repeated measures
  set.seed(300)
  a_ixn <- data.frame(comparison=grp1,
                      tidy_permanova_repeated_measures(dist_mat, s_toTest, grp1, grp2=NULL,reps,perm=perm)[1,])
  combs <- combn(unique(s_toTest[[grp1]]), 2)
  num_tests <- dim(combs)[2]

  # do post hoc tests
  if (a_ixn$p.value < p_cutoff) {
    post_hocs <- lapply(1:num_tests,
                        function(x) data.frame(
                          comparison = paste(combs[,x], collapse=' - '),
                          tidy_permanova_repeated_measures(dist_mat,
                                                           s_toTest[is.element(s_toTest[[grp1]], combs[,x]),], grp1,grp2=NULL,reps, perm=perm)[1,]) )
    a_ixn <- rbind(a_ixn, do.call(rbind, post_hocs))
  }
  a_ixn[,!(colnames(a_ixn) %in% "Term")]
}

#### tidy permanova repeated measures posthoc series ####


#### tidy_anova_repeated_measure_posthoc series ####
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
#### tidy_anova_repeated_measure_posthoc series ####

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

filter_contam <- function(adf){
  # this function is used to filter out chloroplast, mitochondiral, archaea and unassigned sequences
  is_mitochondrial <- grepl("mitochondria", adf$Family)
  is_chloroplast <- grepl("Chloroplast", adf$Class)
  is_unassigned <- grepl("Unassigned", adf$Kingdom)
  is_archaea <- grepl("Archaea", adf$Kingdom)
  is_contam <- is_mitochondrial | is_chloroplast | is_unassigned | is_archaea
  return(is_contam)
}

#### useful functions from ceylan ####
filter_low_coverage <- function(props, perc_cutoff){
  frac_nonzero <- function (x) sum(x > 0) / length(x)
  apply(props, 1, frac_nonzero) >= perc_cutoff
}
rearrange_props <- function(props, sampleIDs) {
  props[,match(sampleIDs, colnames(props))] %>%
    subset(rowSums(.) != 0)
}
rearrange_sample_sheet <- function(s, sampleIDs) {
  s[match(sampleIDs, s$SampleID),]
}
#### useful functions from ceylan ####
