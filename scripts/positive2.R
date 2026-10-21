#!/usr/bin/env Rscript
# Does a consensus of published inflamed signatures beat the individual signatures,
# and does supervised learning beat the consensus? Everything evaluated out-of-cohort.
suppressMessages({library(glmnet)}); set.seed(42)
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
CO <- readRDS(file.path(D,"positive_scores.rds"))
INFLAMED <- c("Chemokine_Messina","STAT1_Care","ImmuneScore_Roh","TIS_Damotte","APM_Thompson",
              "Inflammatory_Thompson","IRG_Ayers","APM_Wang","IFN_Ayers","Teff-IFNG_Fehrenbacher",
              "Tcell-inflamed_Ayers","CTL_Jiang","CYT_Davoli","CD39-CD8Tcell_Chow")
auc <- function(y,v){ ok<-is.finite(v); y<-y[ok]; v<-v[ok]
  if(length(unique(y))<2) return(NA_real_); r<-rank(v); n1<-sum(y==1); n0<-sum(y==0)
  (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
zs <- function(v){ s<-sd(v,na.rm=TRUE); if(!is.finite(s)||s==0) s<-1; (v-mean(v,na.rm=TRUE))/s }
common <- Reduce(intersect, lapply(CO,function(z) colnames(z$S)))
INF <- intersect(INFLAMED, common)
Zl <- lapply(CO, function(z){ Z<-as.data.frame(lapply(z$S, zs)); colnames(Z)<-colnames(z$S); Z[,common,drop=FALSE] })
Y  <- lapply(CO, function(z) z$y); NM <- names(CO)

## ---- 1. consensus (no fitting) vs every individual published signature ----
cons <- lapply(NM, function(n) rowMeans(Zl[[n]][,INF,drop=FALSE])); names(cons)<-NM
a_cons <- sapply(NM, function(n) auc(Y[[n]], cons[[n]]))
tab <- data.frame(signature=common,
  auc_med = sapply(common, function(s) median(sapply(NM, function(n){ a<-auc(Y[[n]],Zl[[n]][[s]]); max(a,1-a) }),na.rm=TRUE)),
  win_vs_consensus = sapply(common, function(s) sum(sapply(NM, function(n){ a<-auc(Y[[n]],Zl[[n]][[s]]); max(a,1-a) }) > a_cons, na.rm=TRUE)),
  p = sapply(common, function(s){ v<-sapply(NM, function(n){ a<-auc(Y[[n]],Zl[[n]][[s]]); max(a,1-a) })
                                  suppressWarnings(wilcox.test(v, a_cons, paired=TRUE)$p.value) }))
tab <- tab[order(-tab$auc_med),]
cat(sprintf("=== CONSENSUS of %d inflamed signatures: median out-of-cohort AUC = %.3f  (%d cohorts) ===\n",
            length(INF), median(a_cons), length(NM)))
cat("\nper-cohort:\n"); print(round(setNames(a_cons,NM),3))
cat("\n=== every published signature vs the consensus (paired across cohorts) ===\n")
print(head(tab,12), row.names=FALSE, digits=3)
cat(sprintf("\n  signatures beating the consensus in >6/%d cohorts: %d / %d\n",
            length(NM), sum(tab$win_vs_consensus>length(NM)/2), nrow(tab)))
cat(sprintf("  signatures significantly better than consensus (p<0.05): %d\n", sum(tab$p<.05 & tab$auc_med>median(a_cons), na.rm=TRUE)))

## ---- 2. supervised, strictly leave-one-cohort-out ----
X <- do.call(rbind, lapply(NM, function(n) as.matrix(Zl[[n]])))
yy <- unlist(lapply(NM, function(n) Y[[n]]))
grp <- rep(NM, sapply(NM, function(n) length(Y[[n]])))
sup <- setNames(rep(NA_real_,length(NM)), NM)
for(te in NM){
  tr <- grp!=te
  cv <- try(cv.glmnet(X[tr,],yy[tr],family="binomial",alpha=.5,nfolds=8,standardize=FALSE),silent=TRUE)
  if(inherits(cv,"try-error")) next
  p <- as.numeric(predict(cv, X[grp==te,,drop=FALSE], s="lambda.min"))
  sup[te] <- auc(Y[[te]], p)
}
cat("\n=== supervised elastic net on the 49 signature scores, leave-one-cohort-out ===\n")
print(round(sup,3))
cat(sprintf("  median = %.3f   vs consensus %.3f\n", median(sup,na.rm=TRUE), median(a_cons)))
cat("  paired Wilcoxon supervised vs consensus: ",
    format.pval(suppressWarnings(wilcox.test(sup,a_cons,paired=TRUE)$p.value),3),"\n")

out <- data.frame(cohort=NM, n=sapply(NM,function(n) length(Y[[n]])),
                  consensus=a_cons, supervised=sup,
                  TIS=sapply(NM,function(n) auc(Y[[n]],Zl[[n]][["TIS_Damotte"]])),
                  IFN_Ayers=sapply(NM,function(n) auc(Y[[n]],Zl[[n]][["IFN_Ayers"]])),
                  GEP_like=sapply(NM,function(n) auc(Y[[n]],Zl[[n]][["Tcell-inflamed_Ayers"]])))
write.csv(out, file.path(D,"positive_consensus.csv"), row.names=FALSE)
write.csv(tab, file.path(D,"positive_vs_published.csv"), row.names=FALSE)
cat("\n=== head-to-head, median out-of-cohort AUC ===\n")
for(c in c("consensus","supervised","TIS","IFN_Ayers","GEP_like"))
  cat(sprintf("  %-12s %.3f\n", c, median(out[[c]],na.rm=TRUE)))
for(c in c("TIS","IFN_Ayers","GEP_like"))
  cat(sprintf("  consensus vs %-10s gain=%+.3f  p=%s  wins %d/%d cohorts\n", c,
      median(out$consensus-out[[c]],na.rm=TRUE),
      format.pval(suppressWarnings(wilcox.test(out$consensus,out[[c]],paired=TRUE)$p.value),3),
      sum(out$consensus>out[[c]],na.rm=TRUE), nrow(out)))
