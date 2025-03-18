###############################################################################
# This script examines the data downloaded to:
## Identify any duplicate participants over the two downloads for Removal
## Ensure proper representation on variables if interest
#

library(dplyr)


n517 <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\imaging data download\\DL_ID_n517.csv",header=FALSE)

n527 <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\imaging data download\\DL_ID_n527.csv",header=FALSE)

NROW(n517) # 374
NROW(n527) # 343

NROW(intersect(n517,n527)) #112 
NROW(setdiff(n517,n527)) # 262 are different and are in the n517 dataset
NROW(setdiff(n527,n517)) # 231 are different and are in the n527 dataset
NROW(intersect(n517,n527)) + NROW(setdiff(n517,n527)) + NROW(setdiff(n527,n517)) # total of 605 unique vars


# write.csv(intersect(n517,n527)$V1,"C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\imaging data download\\intersect_IDs.csv")


# as we can see we have 112 participants that overlap but some that do not 
#> This means we can combine these two downloads together. 
#> so the next step would be to see how all these participant numbers fall on CU traits, gender, race, etc

# stacking all IDs we downloaded. 
all_dl_ids <- rbind(n517,n527)
all_dl_ids <- all_dl_ids[which(!duplicated(data.frame(all_dl_ids))),]
# any duplicated now?
any(duplicated(all_dl_ids)) # no 

# how many?
NROW(all_dl_ids) # 605
#check form above
NROW(all_dl_ids) == (NROW(intersect(n517,n527)) + NROW(setdiff(n517,n527)) + NROW(setdiff(n527,n517)))
# yes we have the same # of participants 

# write.csv(all_dl_ids,"C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\imaging data download\\3165DL_IDs.csv")



# reading data 
cbcl_p <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Measure reconstruction measures\\abcd_cbcl01.csv", header=TRUE)
cbcl_p <- cbcl_p[which(cbcl_p$eventname == 'baseline_year_1_arm_1'),c("src_subject_id",
                                                                      "sex",
                                                                      "cbcl_q94_p",
                                                                      "cbcl_q50_p",
                                                                      "cbcl_q52_p")]

psb_p <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Measure reconstruction measures\\psb01.csv", header=TRUE)
psb_p <- psb_p[which(psb_p$eventname == 'baseline_year_1_arm_1'),c("src_subject_id",
                                                                   "sex",
                                                                   "prosocial_q1_p",
                                                                   "prosocial_q2_p",
                                                                   "prosocial_q3_p")]
psb_p$src_subject_id<-as.character(psb_p$src_subject_id)

# test if row lengths are same
NROW(psb_p) == NROW(cbcl_p) ## TRUE
NROW(cbcl_p)



## joining 
library(dplyr)
CUdata <- left_join(cbcl_p,psb_p,by=c("src_subject_id", "sex"))
CUdata$sex_male <- car::recode(CUdata$sex, "'M'=1; 'F' = 0") 
str(CUdata)

# test of row lengths are the same
NROW(CUdata) == NROW(cbcl_p) # TRUE

#### proper coding of the CU items _______________________________________________

# reverse coding
CUdata$cbcl_q94_p2 = car::recode(CUdata$cbcl_q94_p, '0=2; 1=1; 2=0')
CUdata$cbcl_q50_p2 = car::recode(CUdata$cbcl_q50_p, '0=2; 1=1; 2=0')
CUdata$cbcl_q52_p2 = car::recode(CUdata$cbcl_q52_p, '0=2; 1=1; 2=0')
CUdata$prosocial_q1_p2 = car::recode(CUdata$prosocial_q1_p, '0=2; 1=1; 2=0')
CUdata$prosocial_q2_p2 = car::recode(CUdata$prosocial_q2_p, '0=2; 1=1; 2=0')
CUdata$prosocial_q3_p2 = car::recode(CUdata$prosocial_q3_p, '0=2; 1=1; 2=0')


# collapsing variables
CUdata$cbcl_q26_p3 = car::recode(CUdata$cbcl_q26_p, 'c(1,2)=1; else=0')
CUdata$cbcl_q94_p3 = car::recode(CUdata$cbcl_q94_p2, 'c(1,2)=1; else=0')
CUdata$cbcl_q50_p3 = car::recode(CUdata$cbcl_q50_p2, 'c(1,2)=1; else=0')
CUdata$cbcl_q52_p3 = car::recode(CUdata$cbcl_q52_p2, 'c(1,2)=1; else=0')
CUdata$prosocial_q1_p3 = car::recode(CUdata$prosocial_q1_p2, 'c(1,2)=1; else=0')
CUdata$prosocial_q2_p3 = car::recode(CUdata$prosocial_q2_p2, 'c(1,2)=1; else=0')
CUdata$prosocial_q3_p3 = car::recode(CUdata$prosocial_q3_p2, 'c(1,2)=1; else=0')

# creating a score of CU traits.
  # reconstructing the CU traits measure according to Howes et al. 2020
CUdata$CU_sum <- CUdata$cbcl_q26_p3  + CUdata$prosocial_q1_p3 + CUdata$prosocial_q2_p3 + CUdata$prosocial_q3_p3


# keep downladed data IDs
CUdata<-CUdata[which(CUdata$src_subject_id %in% all_dl_ids),]

NROW(CUdata) # 605
  # same # of rows? 
NROW(CUdata) == NROW(all_dl_ids)# True

# distribution of values 
  ## participant n for each CU level
a1<-as.data.frame(rbind(as.character(table(CUdata$CU_sum)), round(table(CUdata$CU_sum)/sum(table(CUdata$CU_sum)),3)))
rownames(a1) <- c("n =", "% =")
a1

  ## CU level by sex
table(CUdata$sex,CUdata$CU_sum)

  # n of males and females overall
a2<-as.data.frame(rbind(as.character(table(CUdata$sex)), round(table(CUdata$sex)/sum(table(CUdata$sex)),3)))
rownames(a2) <- c("n =", "% =")
a2

# examine distribution of CU
  # all CU
psych::describe(CUdata$CU_sum)
hist(CUdata$CU_sum, breaks = 5)

  # only 1-3 bc we have more 0
psych::describe(CUdata$CU_sum[CUdata$CU_sum>0])
hist(CUdata$CU_sum[CUdata$CU_sum>0], breaks = 5)
  # the distribution is slightly kurtotic and more toward the 0 value...
    # I may have to randomly sample participants from the 0 group 



              #### check pubertal stage and age ####

pub <- read.delim("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Data_Behavioral102020\\abcd_ppdms01.txt",header=TRUE)
key_pub <- t(pub[1,])
pub<-pub[-1,]
# names(pub)
pub <- pub[which(pub$eventname == 'baseline_year_1_arm_1'),]
pub <- pub[which(pub$src_subject_id %in% all_dl_ids),
           c("src_subject_id",
             "sex",
             "interview_age",
             "pds_1_p",
             "pds_2_p",
             "pds_3_p",
             "pds_m4_p",
             "pds_m5_p",
             "pds_f4_p",
             "pds_f5b_p")]

NROW(pub)

# any duplicates
any(duplicated(pub$src_subject_id)) # FALSE: 0 duplicated 

pub <- na_if(pub,999)
pub <- na_if(pub,'')




  ## creating a male and female puberty var
male<-as.numeric(pub$pds_1_p) + as.numeric(pub$pds_2_p) + as.numeric(pub$pds_3_p) + as.numeric(pub$pds_m4_p) + as.numeric(pub$pds_m5_p)
sum(!is.na(male))

female <- as.numeric(pub$pds_1_p) + as.numeric(pub$pds_2_p) + as.numeric(pub$pds_3_p) + as.numeric(pub$pds_f4_p) + as.numeric(pub$pds_f5b_p)
sum(!is.na(female))

  # creating one column for both male and female puberty
pub$puberty<-rep(NA,nrow(pub))
pub$puberty<-ifelse(pub$sex =="F", female, male)

  # checking NAs
    # how many are not NA
length(which(!is.na(pub$puberty))) # 554 are not NA
    # how many are NA
length(which(is.na(pub$puberty))) # 51 are NA 

  # describing puberty
table(pub$puberty)
psych::describe(pub$puberty)

    # descrbing puberty by CU traits 
table(pub$puberty,CUdata$CU_sum)
psych::describeBy(pub$puberty,CUdata$CU_sum)

    #> Impressions
        #> Very few participants were in the higher categories of pubertal stage
        #> There are only 12 participants >= 14 and collapsing on those categories doesn't represent the data
        #> We could 
          #> remove the 12 participants.
          #> Or collapse categories above 14 in puberty 

  # I could just remove OR collapse these 
which(pub$puberty >=14)
    # I think that makes sense

psych::describe(CUdata$CU_sum[which(pub$puberty >=14)])
table(CUdata$CU_sum[which(pub$puberty >=14)])
  ## oddly an equal amount of CU traits across pubertal categories. 
    # CU 0:3 all have 3 participants 



# Describing age
table(pub$interview_age)
    # by CU traits 
table(pub$interview_age,CUdata$CU_sum)
psych::describeBy(as.numeric(pub$interview_age),CUdata$CU_sum)
    ## seems fairly evenly distributed.




              #### Checking race ####

#### Race ####

race <- read.delim("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Data_Behavioral102020\\pdem02.txt",header=TRUE)

key_race <- t(race[1,]) # keeping first row to look at key for each variable
race<-race[-1,] # removing key from dataframe
#names(race)

# Wrangling dataframe
# keeping only varibles we need
race_cols <- append(which(colnames(race)==c("src_subject_id", "eventname")), 
                    which(colnames(race)=="demo_race_a_p___10"):which(colnames(race)=="demo_race_a_p___25"))

race <- race[,race_cols]

# keepoing only baseline collection 
race <- race[which(race$eventname == 'baseline_year_1_arm_1'),]

# testing for subplicates
any(duplicated(race$src_subject_id))  # no duplicates 


# keeping we downloaded 
race<-race[which(race$src_subject_id %in% all_dl_ids),]
NROW(race)  # 605


# creating one race column

race$race_all<- rep(0,nrow(race))
#race_all[which(race$demo_race_a_p___10==1)] <-0 # White
race$race_all[which(race$demo_race_a_p___11==1)] <-1 # Black/African American
race$race_all[which(race$demo_race_a_p___12==1)] <-2 # native american 
race$race_all[which(race$demo_race_a_p___13==1)] <-7 # Alaska native
race$race_all[which(race$demo_race_a_p___14==1)] <-7 # native hawaiian
race$race_all[which(race$demo_race_a_p___15==1)] <-7 # Guamanian 
race$race_all[which(race$demo_race_a_p___16==1)] <-7 # Samoan 
race$race_all[which(race$demo_race_a_p___17==1)] <-7 # Other Pacific Islander
race$race_all[which(race$demo_race_a_p___18==1)] <-8 # Asian Indian
race$race_all[which(race$demo_race_a_p___19==1)] <-9 # Chinese 
race$race_all[which(race$demo_race_a_p___20==1)] <-10 # Filipino
race$race_all[which(race$demo_race_a_p___21==1)] <-14 # Japanese
race$race_all[which(race$demo_race_a_p___22==1)] <-14 # Korean
race$race_all[which(race$demo_race_a_p___23==1)] <-14 # Vietnamese 
race$race_all[which(race$demo_race_a_p___24==1)] <-14 # Other Asian
race$race_all[which(race$demo_race_a_p___25==1)] <-15 # other race

    # because japan, Korean, and vietminese are < 100 we place all of these in "other Asian" 
    # because alaskan, hawaiian, guamanian, Samoan are all < 100 we put all in the category pacific islander

table(race$race_all)

  ## we have a fairly equal distribution across races - mostly white 
    # but we also have distribution across races that is not disproportionate

      # race:     0   1   2   7   8   9  10  14  15 
      # number: 110  98  67  27  43  56  53  70  81 







                  #### Examining removal of items for imaging data ####


#CU traits 
  ## original distribution of CU traits
table(CUdata$CU_sum)
  ## removed for missing resting data 
table(CUdata[-c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$CU_sum)
  ## distribution of those removed
table(CUdata[c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$CU_sum)
    # dist of those removed
    #  0  1  2  3 
    # 11  7  2  7 



# Sex
table(pub$sex)
  ## removed for missing resting data 
table(pub[c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$sex)

table(pub[-c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$sex)




# Puberty
table(pub$puberty)
  ## removed for missing resting data 
table(pub[c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$puberty)

table(pub[-c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$puberty)




# age
table(pub$interview_age)
  ## removed for missing resting data 
table(pub[c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$interview_age)

table(pub[-c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$interview_age)




# race
table(race$race_all)
  ## removed for missing resting data 
table(race[c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$race_all)

table(race[-c(21, 27, 39, 98, 178, 224, 227, 247, 309, 316, 338, 339, 364, 370, 372, 400, 463, 467, 474, 484, 514, 529, 544, 570, 575, 600, 603),]$race_all)























