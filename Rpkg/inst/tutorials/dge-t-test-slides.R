# Code to support the "differences between groups / t-test" slides
# library(mbl2018)
theme_set(theme_bw(base_size = 18))

# First tail test expt: n = 3
set.seed(885)
test1 <- run_tail_length_experiment(num_mice = 4, actual_increase = 1,
                                  variance = 1)
mean1 <-  group_by(test1, group) %>%
  summarize(avg = mean(tail_length))
t.test(tail_length ~ group, test1)

ggplot(test1, aes(x = group, tail_length)) +
  geom_point(size = 4) +
  geom_point(aes(x = group, avg),
               data = mean1,
               color = "red", size = 4) +
  ylab("Tail Length (cm)")

test2 <- run_tail_length_experiment(num_mice = 4, actual_increase = 1,
                                    variance = 1)
mean2 <-  group_by(test2, group) %>%
  summarize(avg = mean(tail_length))
t.test(tail_length ~ group, test2)

ggplot(test2, aes(x = group, tail_length)) +
  geom_point(size = 4) +
  geom_point(aes(x = group, avg),
             data = mean2,
             color = "red", size = 4)

# n = 10
set.seed(885)
test1 <- run_tail_length_experiment(num_mice = 20, actual_increase = 1,
                                    variance = 1)
mean1 <-  group_by(test1, group) %>%
  summarize(avg = mean(tail_length))
t.test(tail_length ~ group, test1)

ggplot(test1, aes(x = group, tail_length)) +
  geom_boxplot() +
  geom_point(size = 4) +
  geom_point(aes(x = group, avg),
             data = mean1,
             color = "red", size = 4) +
  ylab("Tail Length (cm)")

test2 <- run_tail_length_experiment(num_mice = 20, actual_increase = 1,
                                    variance = 1)
mean2 <-  group_by(test2, group) %>%
  summarize(avg = mean(tail_length))
t.test(tail_length ~ group, test2)

ggplot(test2, aes(x = group, tail_length)) +
  geom_boxplot() +
  geom_point(size = 4) +
  geom_point(aes(x = group, avg),
             data = mean2,
             color = "red", size = 4)


pvals <- sapply(1:10000, function(i) {
  dat <- run_tail_length_experiment(20, 0)
  tt <- t.test(tail_length ~ group, dat)
  tt$p.value
})
sum(pvals < 0.05)
mean(pvals < 0.05)

library(edgeR)
ym <- mbl_load_rnaseq("mouse")
ym$design <- model.matrix(~ 0 + group, data = ym$samples)
ymf = ym[filterByExpr(ym, ym$design),]
cm <- makeContrasts(tissue = groupcheek_wildtype - grouppalate_wildtype,
                    levels = ymf$design)
vm <- voom(ymf, ymf$design)
res <- lmFit(vm, vm$design) %>%
  contrasts.fit(cm) %>%
  eBayes() %>%
  topTable(n = Inf)

vmt = mbl_tidy(vm)
dat <- filter(vmt, ens_gene %in% res$ens_gene[1:5],
              group %in% c("groupcheek_wildtype", "grouppalate_wildtype"))

