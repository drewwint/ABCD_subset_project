################################################################################
#                                                                              #
#    Parcelating ABCD resting state data                                       #
#                                                                              #
#    2022/12/31                  Drew E. Winters, PhD.                         #
#                                                                              #
################################################################################


                    ## Prep python and packages ##

Sys.setenv(RETICULATE_PYTHON = "C:\\Program Files\\Python39")


import pandas as pd
import numpy as np
import os, glob, pathlib
import re # to manupulate varaibles 




                          ## Getting IDs ##
        
# extracting all ID
aa = pd.read_csv(r"C:\Users\wintersd\OneDrive - The University of Colorado Denver\1 Publications\ABCD DATA\DL_3165_Code\3165DL_sub-IDs.csv", header= None)
len(aa) # 605

# placing into a list
ab =  []
for i in aa.values.tolist():
  ab.append(str(i[0]))

# Checking data is correct
type(ab) ## type now  list
len(ab) ## 605





                      ## Loops extracting resting state files ##

# These two loop identify resting state files 1 and 2
  #> identify files for each subject ID from above
  #> Idetify those that are missing files "None"

# Resting state T1

run1_files = []
for subnum in range(len(ab)):
  if glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
  run1_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
elif glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
  run1_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
else:
  run1_files.append('None')

# Testing data
len(run1_files) == 605 # True
#Identify which are missing 
np.where(pd.Series(run1_files) == 'None')
# length of those with missing values
len(np.where(pd.Series(run1_files) == 'None')[0]) # 22 participants without functional files
# The list of those ID's that are missing
r1_miss = np.array(aa.iloc[np.where(pd.Series(run1_files) == 'None')].T)[0]


# Resting state T2

run2_files = []
for subnum in range(len(ab)):
  #print(glob.glob("E:\\RKLND_data\\%s\\*%s\\REST_1400_00*\\swurest_1400*.nii" % (subjlist[subnum,0], subjlist[subnum,1])))
  # print(ab[subnum])
  if glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
  run2_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
elif glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
  run2_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-rest_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
else:
  run2_files.append('None')

# Testing data
len(run2_files) == 605 #True
# Identify which are missing 
np.where(pd.Series(run2_files) == 'None')
# Length of those that are missing values
len(np.where(pd.Series(run2_files) == 'None')[0]) # 27 participants without functional files
# The list of those ID's that are missing
r2_miss = np.array(aa.iloc[np.where(pd.Series(run2_files) == 'None')].T)[0]


# seeing which IDs are not in both (should = 27)
## same
len(set(r2_miss).intersection(set(r1_miss))) # 22
## different
len(set(r2_miss).difference(set(r1_miss))) # 5


# making list of IDs to remove from both runs
id_rm = list(set(r2_miss).intersection(set(r1_miss))) + list(set(r2_miss).difference(set(r1_miss)))

ind_rm = ([i for i, e in enumerate(ab) if e in id_rm])


# Removing filepaths not in either list
run1_files_rm = [i for j, i in enumerate(run1_files) if j not in ind_rm]
run2_files_rm = [i for j, i in enumerate(run2_files) if j not in ind_rm]

## Test of length
len(run1_files_rm) == len(run2_files_rm) ## TRUE

len(run1_files_rm) # 578 ## this is 605-27 = 578

## keeping only IDs we didnt remove
id_rm = [i for j, i in enumerate(ab) if j not in ind_rm]
  # Testing data
len(id_rm) == 578 # True


              #### Function for parcelating rest files and concatenating runs ####

# This function does the following for both rest T1 and T2
  #> extract high variance confounds
  #> parcelates into schefer cortical atlas
  #> parcelates into pauli subcortical adles
  #> horiconally stacks the parcelations of the same individual
  #> vertially stacks the timeseries of the same individual
  

cat_files = []
for subnum in range(len(run1_files_rm)): #a[90:]  range(len(run1_files_rm))
  cf1 = high_variance_confounds(run1_files_rm[subnum][0])
  cf2 = high_variance_confounds(run2_files_rm[subnum][0])
  sch_ts1 = schaefer_masker.fit_transform(run1_files_rm[subnum][0],confounds=cf1)
  sch_ts2 = schaefer_masker.fit_transform(run2_files_rm[subnum][0],confounds=cf2)
  pl_ts1 = pauli_masker.fit_transform(run1_files_rm[subnum][0],confounds=cf1)
  pl_ts2 = pauli_masker.fit_transform(run2_files_rm[subnum][0],confounds=cf2)
  ts_set1 = np.hstack((sch_ts1, pl_ts1))
  ts_set2 = np.hstack((sch_ts2, pl_ts2))
  cat_files.append(np.vstack((ts_set1, ts_set2)))

# Testing data size
len(cat_files) == len(run1_files_rm) # True
len(cat_files) #578 
  # all as expected 

# examining shape of each individuals timeseries 
for i in range(len(cat_files)):
  print(i, cat_files[i].shape)
  # some to consider
    # 42 only has 160 timepoints in the ts
    # 491 has 53 timepoints in the ts
  # most have ~760
  # there are a few that are ~450


  # list of those < 30% of what most are (760)
l30 = []
for i in range(len(cat_files)):
  if cat_files[i].shape[0] < (700 - (760 * 0.3)):
    l30.append(i)

l30
len(l30) # 18

low30ID = [i for j, i in enumerate(id_rm) if j in l30]

  # I need to decide what to do with these later


# Writing timeseries  
  # here we are making 
    ##> the columns the names of regions
    ##> saving by participant number
    ##> we start for sp_ for schaefer_pauli_ then subject number
# for i in range(len(cat_files)):
#   pd.DataFrame(cat_files[i], columns = list(np.hstack((schaefer_labels, pauli_labels)))).to_csv(r"C:\Users\wintersd\OneDrive - The University of Colorado Denver\1 Publications\ABCD DATA\schaefer_pauli_timeseries\sp_{}.csv".format(id_rm[i]))


# next to create matricies from timeseries
# cor_mat = []
# for i in cat_files:
#   cor_mat.append(correlation_measure_part.fit_transform([i])[0])















