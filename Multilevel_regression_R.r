# Databricks notebook source
install.packages(c('broom.mixed','margins','lmerTest','purrr'))

# COMMAND ----------

# library(sapbio.ds.ml)
library(dplyr)
library(broom.mixed)
library(data.table)
library(lme4)
library(lmerTest)
library(margins)
library(purrr)

# COMMAND ----------

# MAGIC %md # load data

# COMMAND ----------

data=readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/data/std_impute_raw_data_w_pheno_091922.rdata')
data=setDT(data)

# COMMAND ----------

colnames(data)[startsWith(colnames(data),'coo')]

# COMMAND ----------

# MAGIC %md # filter

# COMMAND ----------

cooled_before_collection <- function(numbers) {
  sapply(
    numbers,
    function(number) {
      if(!is.na(number)){
        if (number>=0) 'Cooled' else 'Not_cooled'  }
        else NA
         
    }
  )
}
    
data$Collection_time=as.numeric(data$Collection_time)
data$Age_cooling=as.numeric(data$Age_cooling)
data$Time_cooled = data$Collection_time - data$Age_cooling 
data[, cooled_before_collection := cooled_before_collection(Time_cooled)]

data[,c('Age_cooling','Collection_time','Time_cooled','cooled_before_collection')]

# COMMAND ----------

T1_filtered_data= data[Sarnat_T1>=0,]
T1_filtered_data$Sarnat=T1_filtered_data$Sarnat_T1
T2_filtered_data= data[Sarnat_T2>=0,]
T2_filtered_data$Sarnat=T2_filtered_data$Sarnat_T2
lmm_data <- rbindlist(list(T1_filtered_data, T2_filtered_data))    # Rbind data.tables

# COMMAND ----------

lmm_data=lmm_data[!is.na(cooled_before_collection)]

# COMMAND ----------

# MAGIC %md # unit test

# COMMAND ----------

lmm_one_mtb=lmm_data[mtb_id %in% 'rLC_neg_mtb_6552530',]
lmm_one_mtb[,.N,by=Sarnat]

# COMMAND ----------

practice_cox_model_dt <- as.data.table(
      tidy(lmer(level ~ Sarnat+Collection_time+factor(cooled_before_collection)+(1|Patient_ID),data=lmm_one_mtb), conf.int = TRUE)
)
practice_cox_model_dt|> filter(effect == "fixed" & term != "(Intercept)")

# COMMAND ----------

summ=lmer(level ~ Sarnat+Collection_time+factor(cooled_before_collection)+(1|Patient_ID),data=lmm_one_mtb)|> summary()
summ$optinfo

# COMMAND ----------

# MAGIC %md ## margin effect

# COMMAND ----------

# lmm_one_mtb$Timepoint=as.factor(lmm_one_mtb$Timepoint)
# lmer(level ~ Sarnat*(Timepoint)+(1|Patient_ID),data=lmm_one_mtb) |> margins() |> tidy(conf.int = TRUE)

# COMMAND ----------

# run_margin_effct_lmm = function(mtb_data){
#   lmer(level ~ Sarnat*(Timepoint)+(1|Patient_ID),data=mtb_data) |> margins() |> tidy(conf.int = TRUE)
# }
# run_margin_effct_lmm(lmm_one_mtb)

# COMMAND ----------

# MAGIC %md ## capture warning

# COMMAND ----------

diagnostic_messages <- lmer(level ~ Sarnat+Collection_time+factor(cooled_before_collection)+(1|Patient_ID),data=lmm_one_mtb) |> summary()
# diagnostic_messages$optinfo$conv$lme4$messages[1]


# COMMAND ----------

d_messages <- diagnostic_messages$optinfo$message[1]
conv_messages <- diagnostic_messages$optinfo$conv$lme4$messages[1]
warnings <- if(length(diagnostic_messages$optinfo$warnings)==0) 'No' else diagnostic_messages$optinfo$warnings[1]
gradient=diagnostic_messages$optinfo$derivs$gradient[1]
Hessian=diagnostic_messages$optinfo$derivs$Hessian[1]

# COMMAND ----------

print(d_messages[1])
print(warnings)
print(derivs_messages)
print(gradient)
print(Hessian)

# COMMAND ----------

gradient=diagnostic_messages$optinfo$derivs$gradient
Hessian=diagnostic_messages$optinfo$derivs$Hessian

# COMMAND ----------

# MAGIC %md # test multiple metabolites

# COMMAND ----------

lmm_4_mtb=lmm_data[mtb_id %in% c('rLC_neg_mtb_6552530',"rLC_pos_mtb_4195645" , "rLC_pos_mtb_3406376" , "rLC_pos_mtb_173202")]
# lmm_4_mtb[,.N,by=Sarnat]
# lmm_4_mtb

# COMMAND ----------

run_lmm <- function(data) {
  formula <- level ~ Sarnat+Collection_time+factor(cooled_before_collection)+(1|Patient_ID)
  model <- lmer(formula, data)
  model_df <-  model |> tidy(conf.int = TRUE)
  model_df$d_messages <- diagnostic_messages$optinfo$message[1]
  model_df$conv_messages <- diagnostic_messages$optinfo$conv$lme4$messages[1]
  model_df$warnings <- if(length(diagnostic_messages$optinfo$warnings)==0) 'No' else diagnostic_messages$optinfo$warnings[1]
  model_df$gradient <- diagnostic_messages$optinfo$derivs$gradient[1]
  model_df$Hessian <- diagnostic_messages$optinfo$derivs$Hessian[1]
  return(model_df)
}

# COMMAND ----------

Res4_lmm_data <- lmm_4_mtb[
  , run_lmm(.SD) |> suppressMessages()
  , by = mtb_id
]

# COMMAND ----------

# Res4_lmm_data[,.N,by=Hessian]
Res4_lmm_data[,.N,by=gradient]
Res4_lmm_data[,.N,by=conv_messages]
# Res4_lmm_data$warnings

# COMMAND ----------

# MAGIC %md # run everything

# COMMAND ----------

run_lmm <- function(data) {
#   formula <- level ~ Sarnat+Collection_time+factor(cooled_before_collection)+(1|Patient_ID)
  formula <- level ~ Sarnat+Collection_time+(1|Patient_ID)
  model <- lmer(formula, data)
  model_df <-  model |> tidy(conf.int = TRUE)
  model_df$d_messages <- diagnostic_messages$optinfo$message[1]
  model_df$conv_messages <- diagnostic_messages$optinfo$conv$lme4$messages[1]
  model_df$warnings <- if(length(diagnostic_messages$optinfo$warnings)==0) 'No' else diagnostic_messages$optinfo$warnings[1]
  model_df$gradient <- diagnostic_messages$optinfo$derivs$gradient[1]
  model_df$Hessian <- diagnostic_messages$optinfo$derivs$Hessian[1]
  return(model_df)
}

# COMMAND ----------

# lmm_data=lmm_data[!is.na(Sarnat),]
# lmm_data=lmm_data[Sarnat>=0,]
lmm_data[,.N,by=Sarnat]

# COMMAND ----------

Res_lmm_data <- lmm_data[
  , run_lmm(.SD) |> suppressMessages()
  , by = mtb_id
]

# COMMAND ----------

Res_lmm_data[,.N,by=conv_messages]

# COMMAND ----------

# Res4_lmm_data[,.N,by=Hessian]
Res_lmm_data[,.N,by=gradient]

# COMMAND ----------

# saveRDS(Res_lmm_data,'/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/MElmm_std_impute_raw_data_w_pheno_092122_withcooled.rdata')

# COMMAND ----------

saveRDS(Res_lmm_data,'/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/MElmm_std_impute_raw_data_w_pheno_092122_without_cooled.rdata')

# COMMAND ----------


