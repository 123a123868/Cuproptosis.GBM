lasso_cox <- function(xdata, ydata){
# 1. LASSO----   
library(glmSparseNet)
## 1.1 just see---
model_lasso <- glmnet(xdata,
                      Surv(ydata$time, ydata$event),
                      family = "cox",
                      alpha = 1)
plot(model_lasso, label = T)
plot(model_lasso, xvar = "lambda", label = T)  # lasso1
 
## 1.2 creat model---- 
fit <- cv.glmHub(xdata, 
                 Surv(ydata$time, ydata$event),
                 family  = 'cox',
                 lambda = buildLambda(1),
                 network = 'correlation',
                 network.options = networkOptions(cutoff = .3,
                                                  min.degree = .1))
print(fit)
plot(fit)

## 1.3 verity----
coefs.v <- coef(fit, s = 'lambda.min')[,1] %>% { .[. != 0]} 
coefs.v %>% {
  data.frame(gene.name   = names(.),
             coefficient = .,
             stringsAsFactors = FALSE)
} %>%
  arrange(gene.name) %>%
  knitr::kable()
coefs.v
ydata$status = ydata$event
separate2GroupsCox(as.vector(coefs.v),
                   xdata[, names(coefs.v)],
                   ydata,
                   plot.title = 'Full dataset', 
                   legend.outside = FALSE)
save(coefs.v, fit, ydata, xdata,
     file = 'Figure6/out_data//coefs.v_lasso_model.Rdata')
write.table(sort(names(coefs.v)),
            file = 'Figure6/out_data//model_lasso_cg_genes.txt',
            row.names = F,col.names = F,quote = F)

}




