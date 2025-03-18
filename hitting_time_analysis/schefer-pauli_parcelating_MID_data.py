
# MID Task Matricies 


import pandas as pd
import numpy as np
import os, glob, pathlib
import re # to manupulate varaibles 


Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\wintersd\\AppData\\Local\\Programs\\Python\\Python312")



aa = pd.read_csv(r"C:\Users\wintersd\OneDrive - The University of Colorado Denver\1 Publications\ABCD DATA\DL_3165_Code\3165DL_sub-IDs.csv", header= None)
len(aa) # 605


ab =  []
for i in aa.values.tolist():
  ab.append(str(i[0]))




run1_files = []
for subnum in range(len(ab)):
    #print(glob.glob("E:\\RKLND_data\\%s\\*%s\\REST_1400_00*\\swurest_1400*.nii" % (subjlist[subnum,0], subjlist[subnum,1])))
    # print(ab[subnum])
    if glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
      run1_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
    elif glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
      run1_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-1_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
    else:
      run1_files.append('None')

len(run1_files)
np.where(pd.Series(run1_files) == 'None')
len(np.where(pd.Series(run1_files) == 'None')[0]) # 22 participants without functional files
r1_miss = np.array(aa.iloc[np.where(pd.Series(run1_files) == 'None')].T)[0]



run2_files = []
for subnum in range(len(ab)):
    #print(glob.glob("E:\\RKLND_data\\%s\\*%s\\REST_1400_00*\\swurest_1400*.nii" % (subjlist[subnum,0], subjlist[subnum,1])))
    # print(ab[subnum])
    if glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
      run2_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n517\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
    elif glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])):
      run2_files.append(glob.glob(r"F:\Collection_3165\dat_3165_n527\derivatives\abcd-hcp-pipeline\%s\ses-baselineYear1Arm1\func\%s_ses-baselineYear1Arm1_task-MID_run-2_space-MNI_bold.nii.gz" % (ab[subnum],ab[subnum])))
    else:
      run2_files.append('None')


len(run2_files)
np.where(pd.Series(run2_files) == 'None')
len(np.where(pd.Series(run2_files) == 'None')[0]) # 27 participants without functional files
r2_miss = np.array(aa.iloc[np.where(pd.Series(run2_files) == 'None')].T)[0]


# seeing which IDs are not in both (should = 7)
  ## same
len(set(r2_miss).intersection(set(r1_miss)))
  ## different
len(set(r2_miss).difference(set(r1_miss)))





# making list of IDs to remove from both runs
id_rm = list(set(r2_miss).intersection(set(r1_miss))) + list(set(r2_miss).difference(set(r1_miss)))

ind_rm = ([i for i, e in enumerate(ab) if e in id_rm])



run1_files_rm = [i for j, i in enumerate(run1_files) if j not in ind_rm]
run2_files_rm = [i for j, i in enumerate(run2_files) if j not in ind_rm]

  ## Test of length
len(run1_files_rm) == len(run2_files_rm) ## TRUE

len(run1_files_rm) # 489 -- this is less than the rest n of 578







## ATLAS
atlas = "Schaefer200Yeo17Pauli" # msdl or haox or mmp # this is the pauli which is subcortical
atlas_c= "Schaefer200"# thsi is the schaefer atlas which is cortical
schaefer_atlas = datasets.fetch_atlas_schaefer_2018(n_rois=200, yeo_networks=17, resolution_mm=1, data_dir=None, base_url=None, resume=True, verbose=1) #atlas_filename = "MMP1_rois.nii" #Glasser et al., 2016
schaefer_filename = schaefer_atlas.maps
schaefer_labels = schaefer_atlas.labels
schaefer_masker = NiftiLabelsMasker(labels_img=schaefer_filename, standardize=True,
                           memory='nilearn_cache', verbose=5) 
pauli_atlas = datasets.fetch_atlas_pauli_2017()
pauli_filename = pauli_atlas.maps
pauli_labels = pauli_atlas.labels
pauli_masker = NiftiMapsMasker(maps_img=pauli_filename, standardize=True, verbose=5) 

all_labels = np.hstack([schaefer_labels, pauli_labels])
#print(all_labels)



# Turns out 4 is an odd file so I need to remove it
    ## making list of participant numbers  
# a=list(range(len(run1_files_rm_2)))
# a.pop(4)
# a.pop(42) # we do 42 bc we popped 4 already (so it'll be 43)

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

# making sure we have the correct sizes
len(cat_files) == len(run1_files_rm) # True
len(cat_files) #578


  # writing
# for i in range(len(cat_files)):
#  pd.DataFrame(cat_files[i], columns = list(np.hstack((schaefer_labels, pauli_labels)))).to_csv(r"C:\Users\wintersd\OneDrive - The University of Colorado Denver\1 Publications\ABCD DATA\schaefer_pauli_timeseries_MID\sp_{}.csv".format(ids[i]))
