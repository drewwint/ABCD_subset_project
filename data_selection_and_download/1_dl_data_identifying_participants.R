
################################################################################
#                                                                              #
#   Code for matching participants to download from ABCD 3165 collection       #
#                                                                              #
#   Drew E. Winters, PhD.                                                      #
#                                                                              #
################################################################################


# This script is intended to identify participant IDs to download - the identified participants will be used in the downloader script

# Logic for participants selection
  #> We wanted to balance the sample we use for 
      #> CU traits
      #> sex 
  #> We also wanted to ensure a distribution of race
  #> Becasue we use collection 3165 we only have participants with acceptable qc scores
      #> Collection 3165 also took care of scanner concerns from different sites 
  #> We considered but did not end up balancing on
    #> Pubertal stage and age
        #> Given we are only using baseline data - Participants were relatively condensed on these measures
        #> As such, these appear more appropriate as controls
    #> Scanner
        #> It was difficult to balance all scanners and still obtain balanced demographics that were most important
        #> This also seems more appropriate as a control. 
        #> Concerns with scanner collection and preprocessing was addresssed by 3165 collection processing 
          #> see <https://collection3165.readthedocs.io/en/stable/>



# Packages
library(dplyr)
library(car) # to recode variables 




                      #### Behavioral data Set up #### 
# reading data 
cbcl_p <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Measure reconstruction measures\\abcd_cbcl01.csv", header=TRUE)
psb_p <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Measure reconstruction measures\\psb01.csv", header=TRUE)
psb_p$src_subject_id<-as.character(psb_p$src_subject_id)

# 3165 participant IDs
all_3165 <- read.csv("C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\imaging data download\\all_3165_IDs.csv", header = FALSE, fileEncoding="UTF-8-BOM") # file encoding used to prevent errors


## only retaning variables that we need 
cbcl_p <- cbcl_p[,c("src_subject_id", 
                    "eventname",
                    "sex",
                    "cbcl_q26_p", 
                    "cbcl_q50_p", 
                    "cbcl_q52_p", 
                    "cbcl_q94_p")]

psb_p <- psb_p[,c("src_subject_id", 
                  "eventname", 
                  "sex",
                  "prosocial_q1_p", 
                  "prosocial_q2_p", 
                  "prosocial_q3_p")]


## joining 
CUdata <- left_join(cbcl_p, 
                    psb_p, 
                    by= c("src_subject_id", "eventname", "sex"))
str(CUdata)

## recoding sex
CUdata$sex_male <- car::recode(CUdata$sex, "'M'=1; 'F' = 0") 



# keeping only files in the 3165 collection 
  # how many participant IDsintersect 
    # we will keep thsi for later as a test to make sure code is working correctly
(match1 <- NROW(intersect(CUdata$src_subject_id,all_3165$V1))) # 10030
  # keep only baseline b/c that's what in the 3165 data
CUdata <- CUdata[which(CUdata$eventname == 'baseline_year_1_arm_1'),]

CUdata <- CUdata[which(CUdata$src_subject_id %in% all_3165[,1]),]


# test to make sure the number matched is the same as we ID'ed before
if (NROW(CUdata) == match1){
  print("YES the sample size is CORRECT")
} else {
  print("NO the sample size is INCORRECT")
}
  # YES the code is working correctly



#### Proper coding of the items 

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

  # descriptives for CU measure
psych::describe(CUdata$CU_sum)
  # count of values by level of CU
table(CUdata$CU_sum)






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


  # keeping those in 3165 participants
race<-race[which(race$src_subject_id %in% all_3165[,1]),]
NROW(race)  # 10029


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





            #### Joining data and Selecting participants ####
# Joining data
cu_qc <- right_join(CUdata,race, by="src_subject_id")


# Function to randomly select participants 
  #> first we create an empty dataframe to plac values into
  #> Then we select participats randomly from each level of CU (0-4) for each sex and race
  #> We are aiming for 500 participants total
    #> so we use the specifier 500/5 (for 5 levels of CU traits)
    #> and we use 500/5/2 for equal numbers in each sex
    #> as we obtain a random sample of race at each of these levels
  #> This is done for each level of CU traits which is appended to the selected participant data frame

# note
  #> because the first download only produce 324 of the 527 participants identified 
  #> we ran this a second time and downloaded 317 of the 517 IDs identified the second time 
  #> We removed duplicated participant IDs and ended up with a final sample of 605
  

selected_part <- data.frame()
for (i in unique(cu_qc$race_all)){
  set.seed(1617)
  CU0f<-cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==0 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==0 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==0 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==0 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU0f<-data.frame(CU0f)
  selected_part <- rbind(selected_part,CU0f)
  CU0m<-cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==1 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==1 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==1 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==0 & cu_qc$sex_male==1 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU0m<-data.frame(CU0m)
  selected_part <- rbind(selected_part,CU0m)
  
  set.seed(1617)
  CU1f<-cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==0 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==0 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==0 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==0 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU1f<-data.frame(CU1f)
  selected_part <- rbind(selected_part,CU1f)
  CU1m<-cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==1 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==1 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==1 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==1& cu_qc$sex_male==1 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU1m<-data.frame(CU1m)
  selected_part <- rbind(selected_part,CU1m)
  
  set.seed(1617)
  CU2f<-cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==0 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==0 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==0 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==0 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU2f<-data.frame(CU2f)
  selected_part <- rbind(selected_part,CU2f)
  CU2m<-cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==1 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==1 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==1 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==2& cu_qc$sex_male==1 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU2m<-data.frame(CU2m)
  selected_part <- rbind(selected_part,CU2m)
  
  set.seed(1617)
  CU3f<-cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==0 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==0 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==0 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==0 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU3f<-data.frame(CU3f)
  selected_part <- rbind(selected_part,CU3f)
  CU3m<-cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==1 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==1 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==1 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==3& cu_qc$sex_male==1 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU3m<-data.frame(CU3m)
  selected_part <- rbind(selected_part,CU3m)
  
  set.seed(1617)
  CU4f<-cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==0 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==0 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==0 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==0 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU4f<-data.frame(CU4f)
  selected_part <- rbind(selected_part,CU4f)
  CU4m<-cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==1 & cu_qc$race_all==i),][sample(nrow(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==1 & cu_qc$race_all==i),]), size = ifelse(NROW(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==1 & cu_qc$race_all==i),]$race_all==i)>=round(((500/5)/2)/length(unique(cu_qc$race_all)),0), round(((500/5)/2)/length(unique(cu_qc$race_all)),0),max(NROW(cu_qc[which(cu_qc$CU_sum==4& cu_qc$sex_male==1 & cu_qc$race_all==i),]))), replace = FALSE), ]
  CU4m<-data.frame(CU4m)
  selected_part <- rbind(selected_part,CU4m)
  
}

  # list of participant IDs
(IDs <- cu_pc$src_subject_id)

# write.csv(IDs, "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\1 Publications\\ABCD DATA\\Data_Behavioral102020\\3165_IDs.csv")

