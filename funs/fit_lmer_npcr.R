#make a function which can run over each imputed dataset and run the mixed model for npcr
fit_lmer_npcr <- function(data) {
    lmer(npcr ~ visit_time*group + baseline_npcr + (1|id),
         data = data)
}
