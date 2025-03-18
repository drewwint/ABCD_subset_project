###############################################################################
### Hitting time analysis in the ABCD dataset #################################
###############################################################################
### Drew E. Winters, PhD. 
# 
# Dynamic diffefrential covaraince (DDC) calculation code
## Origional paper: https://doi.org/10.1073/pnas.2117234119
## Translated from matlab code: https://github.com/yschen13/DDC
## DDC recovery verification code collab notebook: 
##   https://colab.research.google.com/drive/1_e4nKuPXJteBdBFcY7ExtKEFtCHAEEX6
#
###############################################################################


# specifying python to use in reticulate
Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\wintersd\\AppData\\Local\\Programs\\Python\\Python312")


# Packages
import pandas as pd
import numpy as np
from scipy import stats
import os, glob, pathlib
import re # to manupulate varaibles 
from nilearn.connectome import ConnectivityMeasure
from sklearn.linear_model import Ridge
import time


# DDC calculation functions
  # set up functions
    #> derivative_123
    #> dCov_numerical
    #> estimators
  # functions for calculations 
    #> dCov_linear_Reg() for linear ddc
    #> dcov_relu() for nonlinear ddc
  # Note: each of these functions build on eachother
    # after specifying all functions you can just run 
      #> dCov_linear_Reg() for linear ddc
      #> dcov_relu() for nonlinear ddc



def derivative_123(f, dm, dt):
    """
    Compute the first, second, and third derivatives of a function.
    Parameters:
    - f (numpy.ndarray): Array of function values.
    - dm (int): Parameter for the derivative calculation.
    - dt (float): Time step size.
    Returns:
    - D1 (numpy.ndarray): First derivative.
    - D2 (numpy.ndarray): Second derivative.
    - D3 (numpy.ndarray): Third derivative.
    """
    t = np.arange(1 + dm, len(f) - dm)
    # for this function to work f has to be a numpy array
    f = np.array(f)
    # D1 (First Derivative)
    D1 = 0
    d = 0
    for n1 in range(1, dm + 1):
        for n2 in range(n1 + 1, dm + 1):
            d += 1
            numerator = -((f[t - n2] * n1**3 - f[t + n2] * n1**3 -
                f[t - n1] * n2**3 + f[t + n1] * n2**3))
            denominator = (2 * dt * n1**3 * n2 - 2 * dt * n1 * n2**3)
            D1 += numerator / denominator
    D1 /= d

    # D2 (Second Derivative)
    D2 = 0
    d = 0
    for n1 in range(1, dm + 1):
        for n2 in range(n1 + 1, dm + 1):
            d += 1
            numerator = (f[t - n2] * n1**4 + f[t + n2] * n1**4 -
                         f[t - n1] * n2**4 - f[t + n1] * n2**4 -
                         2 * f[t] * (n1**4 - n2**4))
            denominator = (dt**2 * n2**2 * (n1**4 - n1**2 * n2**2))
            D2 += numerator / denominator
    D2 /= d

    # D3 (Third Derivative)
    D3 = 0
    d = 0
    for n1 in range(1, dm + 1):
        for n2 in range(n1 + 1, dm + 1):
            for n3 in range(n2 + 1, dm + 1):
                d += 1
                numerator = (3 * (f[t - n3] * n1 * n2 * (n1**4 - n2**4) +
                                  f[t + n3] * (-n1**5 * n2 + n1 * n2**5) +
                                  n3 * ((f[t - n1] - f[t + n1]) * n2 * (n2**4 - n3**4) +
                                        f[t + n2] * (n1**5 - n1 * n3**4) +
                                        f[t - n2] * (-n1**5 + n1 * n3**4))))
                denominator = (dt**3 * n1 * (n1**2 - n2**2) * n3 * (n1**2 - n3**2) * (n2**3 - n2 * n3**2))
                D3 += numerator / denominator

    D3 /= d
    return D1, D2, D3



def dCov_numerical(cx, h, dm=4):
  """
  Compute various numerical differential covariance matrices using different finite difference methods (derivative approximations).

  Parameters:
  - cx (numpy.ndarray): Time series data matrix (T x N).
  - h (float): Step size for numerical differentiation.
  - dm (int, optional): Maximum window size for custom derivative computation. Default is 4.

  Returns:
  - dCov1 (numpy.ndarray): Differential covariance using first difference method (N x N) [approximation of first derivatie].
  - dCov2 (numpy.ndarray): Differential covariance using central difference method (N x N) [approximation of second derivative].
  - dCov5 (numpy.ndarray): Differential covariance using a higher-order difference method (N x N) [approximation of five-point stencil approximation].
  - dCov_center (numpy.ndarray): Differential covariance using a custom centered method (N x N) [the centered derivative from Taylor expansion].

  NOTE:
  	dCov = <dv/dt,v>
  	covariance is computed through np.cov()
  """
  T, N = cx.shape
  # Compute first-order derivative
  cx = np.array(cx)
  diff_cx = np.vstack(((cx[1:, :] - cx[:-1, :]) / h,
                       np.mean(cx[1:, :] - cx[:-1, :], axis=0)))
  Csample = np.cov(np.hstack((diff_cx, cx)).T, rowvar=False)
  dCov1 = zscore(Csample[:N, N:N + N])

  # Compute second-order derivative
  diff_cx = np.vstack(((1/2 * cx[2:, :] - 1/2 * cx[:-2, :]) / h,
                       np.mean(1/2 * cx[2:, :] - 1/2 * cx[:-2, :], axis=0),
                       diff_cx[-1, :]))
  Csample = np.cov(np.hstack((diff_cx, cx)), rowvar=False)
  dCov2 = zscore(Csample[:N, N:N + N])

  # Compute five-point stencil derivative
  diff_cx = np.zeros((T, N))
  for i in range(2, T - 2):
      for j in range(2, N - 2):
        diff_cx[i, j] = (
            -cx[i + 2, j] + 8 * cx[i + 1, j]
            - 8 * cx[i - 1, j] + cx[i - 2, j]
        ) / 12.0
  Csample = np.cov(np.hstack((diff_cx, cx)), rowvar=False)
  dCov5 = Csample[:N, N:N + N]

  # Compute centered derivative using derivative_123
  diff_cx = np.zeros(((T - 2 * dm)-1, N))
  for i in range(N):
      dx, _, _ = derivative_123(cx[:, i], dm, h)
      diff_cx[:, i] = dx
  cx_trunc = cx[1+dm:T - dm, :]
  Csample = np.cov(np.hstack((diff_cx, cx_trunc)), rowvar=False)
  dCov_center = Csample[:N, N:N + N]

  return dCov1, dCov2, dCov5, dCov_center


def estimators(V, thres, TR):
  """
	INPUT:
  		V_obs: time points x variables
  		thres: Relu offset
  		TR: sampling interval (seconds)
  OUTPUT:
  		B: ReLu(x),x
  		dCov: dx/dt,x
  """
  T, N = V.shape
  V_obs = zscore(V, axis=0)
  Cov = np.cov(V_obs, rowvar=False)
  precision = np.linalg.inv(Cov)
  # ReLU transformation
  Fx = np.maximum(V_obs - thres, 0)
  B = np.cov(Fx,rowvar=False)
  dV = np.diff(V_obs) / TR
  dV = np.vstack([np.mean(dV, axis=0), dV[1:-1], np.mean(dV, axis=0)])
  dCov = np.cov(np.hstack([dV, V_obs]), rowvar=False)

  return Cov, precision, B, dCov


def dCov_linear_Reg(V, TR, lambda_=0.01):
    """
    L2 Regularized version of deltaL.

    Parameters:
    - V (numpy.ndarray): Time series data (T x N).
    - TR (float): Sampling interval.
    - lambda_ (float): Regularization strength.

    Returns:
    - Delta_L (numpy.ndarray): Linear DDC with ridge regularization.
    """
    T, N = V.shape

    # Standardize V
    V_obs = zscore(V, axis=0)

    # Calculate numerical derivatives and estimators
    _, dCov2, _, _ = dCov_numerical(V_obs, TR)
    Cov, _, _, _ = estimators(V_obs, 0, TR)

    # Use Cov for B as per the MATLAB code
    B = Cov
    C = dCov2

    # Initialize the output matrix
    Delta_L = np.zeros_like(C)

    # Perform ridge regression for each row of C
    BtB_reg_inv = np.linalg.inv(B.T @ B + lambda_ * np.eye(B.shape[1]))
    Delta_L = (zscore(np.flipud((BtB_reg_inv @ B.T @ C.T).T))/(N-1))

    return Delta_L


def dcov_relu(V, TR, lambda_=0.01):
    """
    Regularized version of Delta_ReLU
    INPUT:
        V: time series
        TR: sampling interval
        lambda_: regularization strength
    OUTPUT:
        Delta_ReLU: nonlinear DDC with ReLU nonlinearity
    """
    T, N = V.shape

    # Standardize V
    V_obs = zscore(V, axis=0)

    # Calculate numerical derivatives and estimators
    _, dCov2, _, _ = dCov_numerical(V_obs, TR)

    # Apply ReLU non-linearity to the standardized observations
    V_relu = np.maximum(V_obs, 0)  # Applying ReLU

    # Re-compute the estimators using ReLU-transformed data
    Cov, _, B, _ = estimators(V_relu, 1, TR)

    # Ridge regression using the ReLU-transformed covariance matrix
    BtB_reg_inv = np.linalg.inv(B.T @ B + lambda_ * np.eye(B.shape[1]))
    Delta_ReLU = (zscore(np.flipud((BtB_reg_inv @ B.T @ dCov2.T).T))/(N-1))

    return Delta_ReLU



# Extracting tiemseries data
  # to get subject numbers
ids = pd.read_csv(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\DL_3165_Code\3165DL_sub-IDs_KEPT.csv", header= None)
len(ids) == 578 # 578
ids = ids.iloc[:,0].tolist()

ids_f_df = []
for i in ids:
  ids_f_df.append("NDAR_" + i.lstrip('sub-NDAR'))


  ## resting state timeseries data

ts_files = os.listdir(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\schaefer_pauli_timeseries")
len(ts_files) ==578 #578

ts = []
ts_path = r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\schaefer_pauli_timeseries"
for i in ts_files:
  ncols = len(open(os.path.join(ts_path,  i)).readline().split(','))
  ts.append(pd.read_csv(os.path.join(ts_path,  i), usecols=range(1,ncols),skiprows=1, header= None))


  ## MID timeseries data

ts_files_mid = os.listdir(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\schaefer_pauli_timeseries_MID")
len(ts_files_mid) #489

ts_mid = []
ts_path_mid = r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\schaefer_pauli_timeseries_MID"
for i in ts_files_mid:
  ncols = len(open(os.path.join(ts_path_mid,  i)).readline().split(','))
  ts_mid.append(pd.read_csv(os.path.join(ts_path_mid,  i), usecols=range(1,ncols),skiprows=1, header= None))




  # Calculating DDC
  ## Resting
cor_mat_ddc = []
from joblib import parallel_backend
with parallel_backend('threading', n_jobs=-1): 
  for i in ts:
    cor_mat_ddc.append(dcov_relu(i,0.1)) #dCov_linear_Reg

  ## MID
mid_mat_ddc = []
from joblib import parallel_backend
with parallel_backend('threading', n_jobs=-1): 
  for i in ts_mid:
    mid_mat_ddc.append(dcov_relu(i,0.1)) #dCov_linear_Reg



  ## writing DDC matricies
# for i in range(len(cor_mat_ddc)):
#  pd.DataFrame(cor_mat_ddc[i]).to_csv(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\correlation_matrix_rest_DDC\sp_{}.csv".format(ids[i]))

# for i in range(len(mid_mat_ddc)):
#  pd.DataFrame(mid_mat_ddc[i]).to_csv(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\correlation_matrix_MID_DDC\sp_{}.csv".format(ids[i]))



                        ### loading data ###
# After calculatin we can just load the data to reduce time
  # loading ddc matricies
cor_mat_ddc = []
cor_mat_ddc_path = r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\correlation_matrix_rest_DDC"
cor_mat_ddc_files = os.listdir(cor_mat_ddc_path)
for i in cor_mat_ddc_files:
  ncols = len(open(os.path.join(cor_mat_ddc_path,  i)).readline().split(','))
  cor_mat_ddc.append(pd.read_csv(os.path.join(cor_mat_ddc_path,  i), usecols=range(1,ncols),skiprows=1, header= None))

mid_mat_ddc = []
mid_mat_ddc_path = r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\correlation_matrix_MID_DDC"
mid_mat_ddc_files = os.listdir(mid_mat_ddc_path)
for i in mid_mat_ddc_files:
  ncols = len(open(os.path.join(mid_mat_ddc_path,  i)).readline().split(','))
  mid_mat_ddc.append(pd.read_csv(os.path.join(mid_mat_ddc_path,  i), usecols=range(1,ncols),skiprows=1, header= None))
 



                  ### Calculating Hitting matricies ###

# hitting matrix function
def hitting_matrix(correlation_matrix):
    from joblib import Parallel, delayed
    # start_time = time.perf_counter()
    correlation_matrix = np.array(abs(correlation_matrix))  # Ensure absolute values
    np.fill_diagonal(correlation_matrix, 0)  # Set diagonal to 0

    L = correlation_matrix.shape[0]
    A_matrix = correlation_matrix.copy()

    # Degree matrix
    row_sums = A_matrix.sum(axis=1)
    d_max = row_sums.max()

    # Ensure graph connectivity
    for j in range(L):
      if np.max(A_matrix[j,:]) < .05:
          A_matrix[j,j] = d_max - row_sums[j]
    
    row_sums = A_matrix.sum(axis=1)  # Recalculate after adjustment
    D_inv = np.diag(1.0 / row_sums)
    D_sqrt = np.diag(np.sqrt(row_sums))
    D_sqrt_inv = np.diag(1.0 / np.sqrt(row_sums))

    # Transition probability matrix and Graph Laplacian
    p_matrix = D_inv @ A_matrix
    eye_P = np.eye(L) - p_matrix
    G_Lap_n = D_sqrt @ eye_P @ D_sqrt_inv

    # Eigen decomposition
    eig_val, eig_vec = np.linalg.eigh(G_Lap_n)

    # Precompute reusable quantities
    eig_val_nonzero = eig_val[eig_val > eig_val.min()]
    eig_vec_squared = eig_vec ** 2
    d_total = row_sums.sum()

    def compute_H_row(i):
        H_row = np.zeros(L)
        deg_i = row_sums[i]
        for j in range(L):
            deg_j = row_sums[j]
            t_ij = (
                eig_vec_squared[i, eig_val > eig_val.min()] / deg_i
                - eig_vec[i, eig_val > eig_val.min()]
                * eig_vec[j, eig_val > eig_val.min()]
                / np.sqrt(deg_i * deg_j)
            )
            H_row[j] = np.sum(d_total * t_ij / eig_val_nonzero)
        return H_row

    # Parallelize computation of rows
    with Parallel(n_jobs=-1, backend="loky") as parallel:
        H = np.array(parallel(delayed(compute_H_row)(i) for i in range(L)))
    # end_time = time.perf_counter()
    # print(f"individual time: {end_time - start_time:.2f} seconds")
    return H



# Calculating hitting matrix with the DDC covariance
hit_rest_mat_ddc = []
from joblib import parallel_backend
with parallel_backend('threading', n_jobs=-1):
  start_time = time.perf_counter()
  for i in cor_mat_ddc:
    hit_rest_mat_ddc.append(hitting_matrix(i))
  end_time = time.perf_counter()
  print(f"rest total time: {end_time - start_time:.2f} seconds")


hit_mid_mat_ddc = []
from joblib import parallel_backend
with parallel_backend('threading', n_jobs=-1):
  start_time = time.perf_counter()
  for i in mid_mat_ddc:
    hit_mid_mat_ddc.append(hitting_matrix(i))
  end_time = time.perf_counter()
  print(f"MID total time: {end_time - start_time:.2f} seconds")




## need to figure out what happens with mid_mat_ddc[1] 

# writing

## hitting matricies 
# for i in range(len(hit_rest_mat_ddc)):
#  pd.DataFrame(hit_rest_mat_ddc[i]).to_csv(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\hitting_matrix_rest_DDC\sp_{}.csv".format(ids[i]), header= False, index = False)

# for i in range(len(hit_mid_mat_ddc)):
#  pd.DataFrame(hit_mid_mat_ddc[i]).to_csv(r"D:\CU OneDrive Sync\1 Publications\ABCD DATA\hitting_matrix_MID_DDC\sp_{}.csv".format(ids[i]), header= False, index = False)







