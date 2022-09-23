# Databricks notebook source
library(ggplot2)
library(data.table)

# COMMAND ----------

# MAGIC %md # load phenotype data

# COMMAND ----------

pheno_df=fread('/dbfs/mnt/client-002sap21p026-pepper-neshie/04_data_analysis/data/pheno_raw_data_091922.csv')
pheno_df[1:3,1:4]

# COMMAND ----------

# MAGIC %md # convert phenotype data to numeric

# COMMAND ----------

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

# COMMAND ----------

pheno_df$Collection_time=as.numeric(pheno_df$Collection_time)
pheno_df$Age_cooling=as.numeric(pheno_df$Age_cooling)
pheno_df$Time_cooled = pheno_df$Collection_time - pheno_df$Age_cooling 
pheno_df[, cooled_before_collection := cooled_before_collection(Time_cooled)]

# COMMAND ----------

pheno_df$Sarnat = paste(pheno_df$Sarnat_T1,pheno_df$Sarnat_T2,sep='')
pheno_df$Sarnat |> unique()

# COMMAND ----------

pheno_df$worst_Sarnat=pheno_df$'Overall worst Sarnat grade'
pheno_df$worst_Time=pheno_df$'First time-point worst grade recorded'

# COMMAND ----------

pheno_df[
  , worst_Sarnat:= convert_sarnat_text_to_numeric(worst_Sarnat)
]
pheno_df$worst_Sarnat |> unique()
pheno_df[
  , worst_Time:= convert_time_to_numeric(worst_Time)
]
pheno_df$worst_Time |> unique()


# COMMAND ----------

pheno_df[
  , Sarnat := convert_sarnat_text_to_numeric(Sarnat)
]
pheno_df[,.N,by=Sarnat]

# COMMAND ----------

# MAGIC %md # plots

# COMMAND ----------

pheno_df <- pheno_df[,
  .(Timepoint,
    Collection_time,
    Age_cooling,
    Sarnat,
    worst_Sarnat,
    worst_Time,
    Patient_ID,
    Time_cooled,
    cooled_before_collection)
]
T1_df=pheno_df[Timepoint=='T1']
T2_df=pheno_df[Timepoint=='T2']
peak_df=pheno_df[!duplicated(Patient_ID)]

# COMMAND ----------

# MAGIC %md ## scatter plot (patient x DBS collection time)

# COMMAND ----------

myplot=ggplot(pheno_df, aes(x=Patient_ID, y=factor(Timepoint) ,color=factor(Sarnat), fill=factor(Sarnat))) +
  geom_point(size=0.8) + theme_classic() 
myplot=myplot+scale_color_brewer(palette="Dark2")+ggtitle('Distribution of infants')+ 
theme(legend.text = element_text(size = 20))+
 theme(axis.title.y = element_text(size = 25))  + theme(axis.title.x = element_text(size = 15)) +
 theme(axis.text.y = element_text(size = 20)) + theme(axis.text.x = element_text(size =5)) + theme(
plot.title = element_text(color="black", size=25, face="bold")) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

myplot

# COMMAND ----------

# MAGIC %md ## cooling time

# COMMAND ----------

scatterPlot <- ggplot(T1_df,aes(Collection_time, as.numeric(Age_cooling), color=cooled_before_collection)) + 
  geom_point(size=1) + theme_classic() + scale_y_continuous(name ="Age cooling",breaks=seq(0,7,1))+
  theme(legend.position=c(0,1), legend.justification=c(0,1))

scatterPlot

# COMMAND ----------

# MAGIC %md ## peak sarnat

# COMMAND ----------

scatterPlot <- ggplot(peak_df,aes(worst_Sarnat,worst_Time, color=cooled_before_collection)) + 
  geom_point(size=2) + theme_classic() 

scatterPlot

# COMMAND ----------


