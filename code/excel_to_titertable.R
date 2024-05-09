# format titer table from excel sheet
rm(list = ls())
library(tidyverse)
library(stringr)


set_threshold <- function(tab, thresh = 20) {
  tab[as.numeric(tab) < as.numeric(thresh)] <- paste0("<", thresh)
  tab[is.na(tab)] <- "*"
  
  return(tab)
}


# =============================== Omicron I sheet, published in NEJM
tab <- readxl::read_excel("./data/titer_data/03152024_Update 5.xlsx", sheet = 1)
colnames(tab) <- tab[1,]
tab <- tab[!is.na(tab$`NHP ID`),]
tab <- tab[2:nrow(tab),]

tab %>%
  fill(`Study#`, `Virus strain`, `number of samples`, `Virus dose/note`, `Time point post challenge`) -> tab

# make second sheet with vax cohort
vax <- readxl::read_excel("./data/titer_data/03152024_Update 5.xlsx", sheet = 2)
colnames(vax) <- vax[1,]
vax <- vax[!is.na(vax$`NHP ID`),]
vax <- vax[2:nrow(vax),]

# rename columns of vax such that they fit the infected names
vax_colnames <- c("Time point post challenge" = "Virus strain",
                  "Vaccination" = "Virus dose/note",
                  "Days post last vax dose" = "Time point post challenge"
                  )
colnames(vax)[colnames(vax) %in% names(vax_colnames)] <- vax_colnames[colnames(vax)[colnames(vax) %in% names(vax_colnames)]]

vax %>%
  select(!Antigen) %>%
  fill(`Study#`, `Virus strain`, `number of samples`, `Virus dose/note`, `Time point post challenge`) -> vax

ag_cols <- colnames(tab)[7:ncol(tab)]
# set sr info
rbind(tab, vax) %>%
  mutate(`Virus strain` = gsub("Delta", "delta", `Virus strain`),
          sr_group = ifelse(grepl("vax",`Virus strain`), `Virus strain`, paste0(`Virus strain`, " conv.")),
         sr_group_time = paste0(sr_group, ":",`Time point post challenge`),
         sr_group_dose = paste0(sr_group, ":",`Virus dose/note`),
         sr_info = paste(sr_group, `NHP ID`, `Time point post challenge`, `Virus dose/note`, `Study#`, sep = "_")) -> tab

tab %>%
  pivot_longer(cols = all_of(ag_cols), names_to = "ag_name", values_to = "Titer") %>%
  mutate("Titer_thresh1" = ifelse(as.numeric(Titer) <= 1, "<1", Titer),
         "Titer_thresh20" = ifelse(as.numeric(Titer) < 20, "<20", Titer)) -> tab_long

tab_long[is.na(tab_long$Titer), c("Titer","Titer_thresh1", "Titer_thresh20")] <- "*"


write.csv(tab_long, "./data/titer_data/titer_data_long.csv")
