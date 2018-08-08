# Daniel Alfonsetti (daniel.alfonsetti)
# 8/10/18
# A script for Ron's interesting (kidney?) gene (Kcnn4)
library(tidyverse)

################################
# Helper function

## Summarizes data. (This function was taken from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/)
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
################################



annot_row <- annot.mrna[annot.mrna$symbol == "Kcnn4", ]
current_id <- annot_row$gene_id
symbol <- annot_row$symbol

expr_data <- expr.mrna[, current_id]  # Use ID as key to get expression data from the expr.mrna matrix
df <- cbind(annot.sample$Sex, annot.sample$Age, expr_data)
colnames(df) <- c("Sex", "Age", "expr")

df <- as.data.frame(df)
df <- mutate(df, sex_factor = factor(df$Sex))
levels(df$sex_factor) <- c("Female", "Male")

# Perform t-test on sex
# Seperate groups
males <- df %>% filter(df$Sex == 1)
females <- df %>% filter(df$Sex == 2)
res.t <- t.test(x = males$expr, y = females$expr)

sex.p.value <- res.t$p.value


# Perform ANOVA on age
res.aov <- aov(df$expr ~ df$Age, data = df)
sum_obj <- summary(res.aov)
age.p.value <- sum_obj[[1]][, 5][1]


# Perform ANOVA on age:sex
df$Age_Sex <- paste0(as.character(df$Age), "_", as.character(df$sex_factor))
res.aov <- aov(df$expr ~ df$Age_Sex, data = df)
sum_obj <- summary(res.aov)
age_sex.p.value <- sum_obj[[1]][, 5][1]



# Make plots

setwd("~/do_heart/results/")
pdf("Kcnn4_plots_part2.pdf")


summary_df <- summarySE(df, measurevar="expr", groupvars=c("Age"))
plot <- ggplot(summary_df, aes(x=as.factor(Age), y=expr, group = 1)) +
  geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3) +
  ggtitle(label = paste(symbol, " (", current_id, ") ", "Expression", sep = "")) +
  ylab(label = "Mean Normalized Expression Level") +
  xlab(label = "Age (months)") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = paste0("Age ANOVA p-value: ", round(age.p.value, 5)), hjust = -.1, vjust = 1.5)
print(plot)


summary_df <- summarySE(df, measurevar="expr", groupvars=c("sex_factor"))
plot <- ggplot(summary_df, aes(x=sex_factor, y=expr, color = sex_factor)) +
  theme_minimal() +
  geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3) +
  ggtitle(label = paste(symbol, " (", current_id, ") ", "Expression", sep = "")) +
  ylab(label = "Mean Normalized Expression Level") +
  xlab(label = "Sex") +
  labs(color = "Sex") +
  annotate("text", x = -Inf, y = Inf, label = paste0("Sex ANOVA p-value: ", round(sex.p.value, 5)), hjust = -.1, vjust = 1.5)
print(plot)


  
summary_df <- summarySE(df, measurevar="expr", groupvars=c("sex_factor", "Age"))
plot <- ggplot(summary_df, aes(x=as.factor(Age), y=expr, color = sex_factor, group = sex_factor)) +
  theme_minimal() +
  geom_errorbar(aes(ymin=expr-se, ymax=expr+se), width=.1, position=position_dodge(0.1)) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1), size=3) +
  ggtitle(label = paste(symbol, " (", current_id, ") ", "Expression", sep = "")) +
  ylab(label = "Mean Normalized Expression Level") +
  xlab(label = "Age (months)") +
  labs(color = "Sex") +
  annotate("text", x = -Inf, y = Inf, label = paste0("Age:Sex ANOVA p-value: ", round(age_sex.p.value, 5)), hjust = -.1, vjust = 1.5)
print(plot)
  
  
dev.off()
