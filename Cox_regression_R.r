# Databricks notebook source
install.packages('broom')

# COMMAND ----------

library(ggplot2)
library(data.table)
# library(lmerTest, include.only = "lmer")
library(broom, include.only = "tidy")
# library(broom.mixed, include.only = "tidy")

# COMMAND ----------

data=readRDS('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/data/std_impute_raw_data_w_pheno_091922.rdata')
data=setDT(data)

# COMMAND ----------

colnames(data)[!startsWith(colnames(data),'rLC')]

# COMMAND ----------

# MAGIC %md # filtering and cleaning

# COMMAND ----------

data = setDT(data)
data[, .N, by=WorstGradeTime]

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

convert_time_to_numeric <- function(strings) {
  sapply(
    strings,
    function(string) {
      switch(
        string,
        "<6H" = 3,
        "48H" = 48,
        "6H" = 6,
        "D7" = 24*7,
        "Data not available on source" = NA,
        "D6" = 24*6,
        "72H" = 72,
        "24H" = 24,
        "96H" = 96,
        "D10" = 10*24,
        "D8" = 8*24,
        NA
      )
    }
  )
}
convert_sarnat_text_to_numeric <- function(strings) {
  sapply(
    strings,
    function(string) {
      switch(
        string,
        "Moderate (2)" = 2,
        "Moderate" = 2,
        "Severe (3)" = 3,   
        'Severe' = 3,
        "Mild (1)" = 1,
        'Mild' = 1,
        "Normal (0)" = 0,
        "Normal" = 0,
        NA
      )
    }
  )
}

# COMMAND ----------

data$WorstSarnat |> unique()

# COMMAND ----------

data[
  , worst_Time:= convert_time_to_numeric(WorstGradeTime)
]
data$worst_Time |> unique()

# COMMAND ----------

filtered_data <- data[!is.na(WorstGradeTime),]
filtered_data <- filtered_data[WorstSarnat>=0,]
filtered_data <- filtered_data[cooled_before_collection=='Cooled',]

print(filtered_data[,.N,by=WorstSarnat])
print(filtered_data[,.N,by=Sarnat_T1])

# COMMAND ----------

# MAGIC %md # select for T1
# MAGIC T1 is taken from 0-6 hr, so we assume T1 is at 3 hr

# COMMAND ----------

T1_filtered_data= filtered_data[Sarnat_T1>=0,]
T1_filtered_data= T1_filtered_data[WorstGradeTime>=Collection_time,]
print(T1_filtered_data$'Patient ID' |> unique() |> length())
T1_filtered_data$'Sarnat_T1' |> length()
T1_filtered_data$'WorstSarnat' |> length()

print(T1_filtered_data[,.N,by=WorstGradeTime])
print(T1_filtered_data[,.N,by=WorstSarnat])
print(T1_filtered_data[,.N,by=Sarnat_T1])

# COMMAND ----------

T1_filtered_data <- T1_filtered_data[,
  .(
    Sarnat_baseline=Sarnat_T1,
    WorstSarnat,
    Collection_time,
    WorstGradeTime,
    Sarnat_Change =  as.numeric(WorstSarnat)-as.numeric(Sarnat_T1),
    Time_Taken = as.numeric(WorstGradeTime)-Collection_time,
    Patient_ID,
    mtb_id,
    level
  )
]
T1_filtered_data[,.N,by=Sarnat_Change]

# COMMAND ----------

T1_filtered_data

# COMMAND ----------

T1_filtered_data[Sarnat_Change==-1,.N, by=Patient_ID]

# COMMAND ----------

# MAGIC %md # select for T2

# COMMAND ----------

T2_filtered_data= filtered_data[Sarnat_T1<0 & Sarnat_T2 >=0 ,]
T2_filtered_data= T2_filtered_data[WorstGradeTime>=Collection_time,]
print(T2_filtered_data$'Patient ID' |> unique() |> length())


print(T2_filtered_data[,.N,by=WorstGradeTime])
print(T2_filtered_data[,.N,by=WorstSarnat])
print(T2_filtered_data[,.N,by=Sarnat_T2])

# COMMAND ----------

T2_filtered_data <- T2_filtered_data[,
  .(
    Sarnat_baseline=Sarnat_T2,
    WorstSarnat,
    WorstGradeTime,
    Sarnat_Change =  as.numeric(WorstSarnat)-as.numeric(Sarnat_T2),
    Time_Taken = as.numeric(WorstGradeTime)-as.numeric(Collection_time),
    Patient_ID,
    mtb_id,
    level
  )
]
T2_filtered_data[,.N,by=Sarnat_Change]

# COMMAND ----------

T2_filtered_data

# COMMAND ----------

# MAGIC %md # combine T1 and T2

# COMMAND ----------

cox_data <- rbindlist(list(T1_filtered_data, T2_filtered_data))    # Rbind data.tables
cox_data[,.N,by=Sarnat_Change]

# COMMAND ----------

convert_sarnat_cox <- function(Sarnat_Changes) {
  sapply(
    Sarnat_Changes,
    function(Sarnat_Change) {
      if(Sarnat_Change>0) 1 else 0
    }
  )
}

# COMMAND ----------

cox_data[, cox_score := convert_sarnat_cox(Sarnat_Change)]
cox_data

# COMMAND ----------

# MAGIC %md # cox regression for one

# COMMAND ----------

library("survival")
args(coxph)

# COMMAND ----------

library(broom, include.only = "tidy")
# library(broom.mixed, include.only = "tidy")

# COMMAND ----------

cox_one_mtb=cox_data[mtb_id %in% 'rLC_neg_mtb_6552530',]
cox_one_mtb[,.N,by=cox_score]

# COMMAND ----------

cox_for_one_mtb <- coxph(Surv(Time_Taken, cox_score) ~ level,data=cox_one_mtb)
cox_for_one_mtb

# COMMAND ----------

practice_cox_model_dt <- as.data.table(
  tidy(cox_for_one_mtb, conf.int = TRUE)
)
practice_cox_model_dt

# COMMAND ----------

# MAGIC %md # cox regression for all

# COMMAND ----------

run_cox_model <- function(data) {
  formula <- Surv(Time_Taken, cox_score) ~ level
  model <- coxph(formula, data)
  model_df <- tidy(model, conf.int = TRUE)
#   diagnostic_messages <- summary(model)$optinfo$message
#   model_df$messages <- diagnostic_messages
  return(model_df)
}

# COMMAND ----------

cox_data
cox_data[, .N, by = cox_score]

# COMMAND ----------

metabolite_cox_data_dt <- cox_data[
  , run_cox_model(.SD) |> suppressMessages()
  , by = mtb_id
]
metabolite_cox_data_dt

# COMMAND ----------

# saveRDS(metabolite_cox_data_dt,'/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/results/cox_std_impute_raw_data_w_pheno_091222.rdata')

# COMMAND ----------

# MAGIC %md # result visualization

# COMMAND ----------

data=as.data.frame(metabolite_cox_data_dt)
# data=as.data.frame(output)
# data=data[data$warnings_full=='None',]
data['logp']=-log10(data$p.value)
data['Type']='Below threshold'
data[data$logp>3 & data$estimate_Sarnat > 0 ,'Type']='Elevated above threshold'
data[data$logp>3 & data$estimate_Sarnat < 0 ,'Type']='Reduced above threshold'

print(length(data[data$logp>3 & data$estimate < 0 ,'mtb_id']))
print(length(data[data$logp>3 & data$estimate > 0 ,'mtb_id']))
reduced_cox=data[data$logp>3 & data$estimate < 0 ,'mtb_id']
elevated_cox=data[data$logp>3 & data$estimate > 0 ,'mtb_id']
cox=data[data$logp>3,'mtb_id']

myplot=ggplot(data, aes(x=estimate, y=logp,color=Type))+geom_point(alpha = 0.4 )
myplot
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlim(-2, 2) + 
scale_color_manual(values=c('grey','#039fbe','#cf1578'))+ 
theme(legend.text = element_text(size = 15))+
 theme(axis.title.y = element_text(size = 20))  + theme(axis.title.x = element_text(size = 20)) +
 theme(axis.text.y = element_text(size = 15)) + theme(axis.text.x = element_text(size = 15)) 

# COMMAND ----------

# MAGIC %md # look at the top molecules

# COMMAND ----------

data[order(p.value)]

# COMMAND ----------

cox_data

# COMMAND ----------

mtb='rLC_neg_mtb_3620056'
plot_mtb=cox_data[mtb_id %in% mtb,]
plot_mtb

# COMMAND ----------

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))

# COMMAND ----------

ggplot(plot_mtb, aes(x=factor(Sarnat_Change), y=Time_Taken, colour=level)) + 
  geom_point(size=6)+sc+
 ggtitle(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # # #   xlim(-1000, 1000) + 
#  theme(legend.position="none") +

# theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------

myplot=ggplot(plot_mtb, aes(y=Time_Taken, x=factor(Sarnat_Change), fill=level)) +
# myplot=ggplot(df1, aes(x=Sarnat, y=level, fill=Timepoint)) +
  geom_boxplot()
myplot
# # Add dots
myplot + geom_dotplot(binaxis='y', 
                      stackdir='center', dotsize = 0.1,
                 position=position_dodge(1),f)+

sc+
ggtitle(mtb)+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) +
# # # # # #   xlim(-1000, 1000) + 
 theme(legend.position="none") +

# theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 25)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =20)) + theme(
plot.title = element_text(color="black", size=25, face="bold"))


# COMMAND ----------


