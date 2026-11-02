import numpy as np
import matplotlib.pyplot as plt
import mabs_funcs as mf
from scipy.stats import binned_statistic_2d



Nnc  = 2*(mf.nc10+1 - mf.n_nancut)-1 
Nnc_hf = Nnc//2+1
L   = Nnc * mf.dnu_mwa * mf.rc_prime
rc4 = mf.rc**4
Lsq = L**2
dkll = 2*np.pi/L



def blackman_nuttall(M):
    """
    Return the length‑M Blackman–Nuttall window:
      w[n] = a0 
             - a1 cos(2π n/(M−1))
             + a2 cos(4π n/(M−1))
             - a3 cos(6π n/(M−1))
    with the standard coefficients.
    """
    n = np.arange(M)
    a0, a1, a2, a3 = 0.3635819, 0.4891775, 0.1365995, 0.0106411
    return (a0
            - a1 * np.cos(2*np.pi*n/(M-1))
            + a2 * np.cos(4*np.pi*n/(M-1))
            - a3 * np.cos(6*np.pi*n/(M-1)))


## mabs cylbs without windows
def mabs(r1,r2,r3,file_loc_v2, d_type,Rlzn):  
    #ell1 = 2*np.pi*Uavg_phy[r1]
    #ell2 = 2*np.pi*Uavg_phy[r2]
    #ell3 = 2*np.pi*Uavg_phy[r3]
    # print(r'l1=%f, kp1=%f'%(ell1,ell1/mf.rc))
    # print(r'l2=%f, kp2=%f'%(ell2,ell2/mf.rc))
    # print(r'l3=%f, kp3=%f'%(ell3,ell3/mf.rc))
    # filename = "%s/bsout_final_ra%d_rb%d_rc%d_R%d" % (file_loc_v1, r1, r2, r3, 1)
    # Ntri1, mabs1, numc1 = mf.read_bsout_float(filename)
    filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc_v2, r1, r2, r3, d_type,Rlzn)
    Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
    mabs_F1 = mf.fold_quadrants(mabs1)
    numc_mean_1 = mf.fold_quadrants(numc1)
    mabs_T_1 = mabs_F1/mf.mabs_norm/numc_mean_1
    mabs_final_1=  mf.cut_nan(mabs_T_1,mf.n_nancut)

    return mabs_final_1


def cylbs(r1,r2,r3,file_loc_v2, d_type,Rlzn,bin_flag,n_bins=30):
    mabs_final_F = mabs(r1,r2,r3,file_loc_v2, d_type,Rlzn)
    res_F = mf.Bdnu2Bk_fft(mabs_final_F)
    if (bin_flag==True):
        final_res = bin_2D_cylbs(res_F,n_bins)
    elif (bin_flag==False):
        final_res = res_F
    else: print('Give bin flag\n')

    return final_res
    
    

def mabs_RSP(r1,r2,r3,file_loc_v1, d_type,Rlzn):  ## rings from single plane
    filename = "%s/bsout_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc_v1, r1, r2, r3,  d_type,Rlzn)
    Ntri1, mabs1, numc1 = mf.read_bsout_float(filename)
    #filename = "%s/bsout_final2_ra%d_rb%d_rc%d_R%d" % (file_loc_v2, r1, r2, r3, i)
    #Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
    mabs_F1 = mf.fold_quadrants(mabs1)
    numc_mean_1 = mf.fold_quadrants(numc1)
    mabs_T_1 = mabs_F1/mf.mabs_norm/numc_mean_1
    mabs_final_1=  mf.cut_nan(mabs_T_1,mf.n_nancut)

    return mabs_final_1
    
    
def cylbs_RSP(r1,r2,r3,file_loc_v2, d_type,Rlzn,bin_flag,n_bins=30):
    mabs_final_F = mabs_RSP(r1,r2,r3,file_loc_v2, d_type,Rlzn)
    res_F = mf.Bdnu2Bk_fft(mabs_final_F)
    if (bin_flag==True):
        final_res = bin_2D_cylbs(res_F,n_bins)
    elif (bin_flag==False):
        final_res = res_F
    else: print('Give bin flag\n')

    return final_res   
    
# Windowed Mabs and Cylbs (with folding)

def mabs_win(r1,r2,r3,file_loc_v2, d_type,Rlzn):
    #filename = "%s/bsout_final2_ra%d_rb%d_rc%d_R%d" % (file_loc_v2, r1, r2, r3, 1)
    filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc_v2, r1, r2, r3, d_type,Rlzn)
    Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
    mabs_F1 = mf.fold_quadrants(mabs1)
    numc_mean_1 = mf.fold_quadrants(numc1)
    mabs_T_1 = mabs_F1/mf.mabs_norm/numc_mean_1
    
    M, N = mabs1.shape
    wy   = blackman_nuttall(M - 2*mf.n_nancut)
    wx   = blackman_nuttall(N - 2*mf.n_nancut)
    W2d  = np.outer(wy, wx)
    wind = W2d[wx.shape[0]//2:,:] #mf.fold_quadrants(W2d)/2
    
    mabs_final_1=  mf.cut_nan(mabs_T_1,mf.n_nancut)*wind

    return mabs_final_1
    
    
    
def mabs_RSP_win(r1,r2,r3,file_loc_v1, d_type,Rlzn):  ## rings from single plane
    filename = "%s/bsout_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc_v1, r1, r2, r3,  d_type,Rlzn)
    Ntri1, mabs1, numc1 = mf.read_bsout_float(filename)
    #filename = "%s/bsout_final2_ra%d_rb%d_rc%d_R%d" % (file_loc_v2, r1, r2, r3, i)
    #Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
    mabs_F1 = mf.fold_quadrants(mabs1)
    numc_mean_1 = mf.fold_quadrants(numc1)
    mabs_T_1 = mabs_F1/mf.mabs_norm/numc_mean_1
    
    M, N = mabs1.shape
    wy   = blackman_nuttall(M - 2*mf.n_nancut)
    wx   = blackman_nuttall(N - 2*mf.n_nancut)
    W2d  = np.outer(wy, wx)
    wind = W2d[wx.shape[0]//2:,:] #mf.fold_quadrants(W2d)/2
    
    mabs_final_1=  mf.cut_nan(mabs_T_1,mf.n_nancut)*wind
    

    return mabs_final_1

# def mabs_win2(r1,r2,r3): ### it clips data including window leaving non-zero at edges
#     filename = "%s/bsout_final2_ra%d_rb%d_rc%d_R%d" % (file_loc_v2, r1, r2, r3, 1)
#     Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
#     M, N   = mabs1.shape
#     wy     = blackman_nuttall(M)
#     wx     = blackman_nuttall(N)
#     W2d    = np.outer(wy, wx)
#     mabs_F1= mf.fold_quadrants(mabs1*W2d)
#     # mabs_F1 = mf.fold_quadrants(W2d)
#     numc_mean_1 = mf.fold_quadrants(numc1)
#     mabs_T_1 = mabs_F1/mf.mabs_norm/numc_mean_1
    
#     mabs_final_1=  mf.cut_nan(mabs_T_1,mf.n_nancut)

#     return mabs_final_1


def cylbs_win(r1,r2,r3,file_loc_v2, d_type,Rlzn):
    data_win = mabs_win(r1,r2,r3,file_loc_v2, d_type,Rlzn)
    M, N = data_win.shape
    res = np.fft.fftshift(np.fft.irfft2(np.fft.ifftshift(data_win, axes=1),axes=(1, 0),s=(N, N),))
    
    rc4   = mf.rc**4
    L     = N * mf.dnu_mwa * mf.rc_prime 
    Lsq   = L**2
    res  *= (rc4 * Lsq)

    return res
    
def cylbs_RSP_win(r1,r2,r3,file_loc_v2, d_type,Rlzn):
    data_win = mabs_RSP_win(r1,r2,r3,file_loc_v2, d_type,Rlzn)
    M, N = data_win.shape
    res = np.fft.fftshift(np.fft.irfft2(np.fft.ifftshift(data_win, axes=1),axes=(1, 0),s=(N, N),))
    
    rc4   = mf.rc**4
    L     = N * mf.dnu_mwa * mf.rc_prime 
    Lsq   = L**2
    res  *= (rc4 * Lsq)

    return res   
    
    
    #Windowed mabs and cylbs in Full plane (without fold)
def cut_nan3(arr,n):
    return arr[n:-n,n:-n]

def mabs_win3(r1,r2,r3,file_loc_v2, d_type,Rlzn):
    #filename = "%s/bsout_final2_ra%d_rb%d_rc%d_R%d" % (file_loc_v2, r1, r2, r3, 1)
    filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc_v2, r1, r2, r3, d_type,Rlzn)
    Ntri2, mabs1, numc1 = mf.read_bsout_float_v2(filename)
    M, N = mabs1.shape
    wy = blackman_nuttall(M - 2*mf.n_nancut)
    wx = blackman_nuttall(N - 2*mf.n_nancut)
    W2d      = np.outer(wy, wx)
    mabs_T_1 = mabs1/mf.mabs_norm/numc1
 
    mabs_final_1=  cut_nan3(mabs_T_1,mf.n_nancut)*W2d

    return mabs_final_1
    
def cylbs_win3(r1,r2,r3,file_loc_v2, d_type,Rlzn):
    data_win = mabs_win3(r1,r2,r3,file_loc_v2, d_type,Rlzn)
    M, N = data_win.shape
    res = np.fft.fftshift(
             np.fft.ifft2(
               np.fft.ifftshift(data_win, axes=(0,1)),
               axes=(1,0),
               s=(N, N),
             )
           ).real
    
    rc4   = mf.rc**4
    L     = N * mf.dnu_mwa * mf.rc_prime 
    Lsq   = L**2
    res  *= (rc4 * Lsq)

    return res
    
    
    
    
    
def bin_2D_cylbs(cylbs_mean_NF,n_bins):
    x = np.linspace(-dkll*(Nnc_hf-1), dkll*(Nnc_hf-1), cylbs_mean_NF.shape[0])  #np.linspace(-dkll*(Nnc_hf-1),dkll*(Nnc_hf-1),Nnc )
    y = np.linspace(-dkll*(Nnc_hf-1), dkll*(Nnc_hf-1), cylbs_mean_NF.shape[1])
    X, Y = np.meshgrid(x, y)
    
    # Flatten the arrays to feed into the binned statistic function
    x_flat = X.ravel()
    y_flat = Y.ravel()
    data_flat_NF = cylbs_mean_NF.ravel()
    
    # Bin the data into 30x30 bins using the mean as the statistic
    binned_clbs_NF, x_edges, y_edges, binnumber = binned_statistic_2d( x_flat, y_flat, data_flat_NF, statistic='mean', bins=n_bins)

    return binned_clbs_NF

    
    
    
  
    
from scipy.signal.windows import (
    nuttall,
    blackmanharris,
    flattop,
    parzen,
    bohman,
    kaiser,
    chebwin,
)
def winds(M, win='nutt', **kwargs):
    """
    data_dnu : 2D array of shape (M, N)
    win : one of {'bn','nutt','hann','hamming','blackman','bartlett',
                  'blackharris','flattop','parzen','bohman','kaiser','chebwin'}
    kwargs : for 'kaiser' you can pass beta=…, for 'chebwin' at=…
    """
    N = M

    # pick the 1D windows
    if win == 'bn':
        wy = blackman_nuttall(M)
        wx = blackman_nuttall(N)

    elif win == 'nutt':
        wy = nuttall(M)
        wx = nuttall(N)

    elif win == 'hann':
        wy = np.hanning(M)
        wx = np.hanning(N)

    elif win == 'hamming':
        wy = np.hamming(M)
        wx = np.hamming(N)

    elif win == 'blackman':
        wy = np.blackman(M)
        wx = np.blackman(N)

    elif win == 'bartlett':
        wy = np.bartlett(M)
        wx = np.bartlett(N)

    elif win == 'blackharris':
        wy = blackmanharris(M)
        wx = blackmanharris(N)

    elif win == 'flattop':
        wy = flattop(M)
        wx = flattop(N)

    elif win == 'parzen':
        wy = parzen(M)
        wx = parzen(N)

    elif win == 'bohman':
        wy = bohman(M)
        wx = bohman(N)

    elif win == 'kaiser':
        beta = kwargs.get('beta', 14.0)      # you can tweak this
        wy   = kaiser(M, beta)
        wx   = kaiser(N, beta)

    elif win == 'chebwin':
        at = kwargs.get('at', 100.0)         # sidelobe attenuation in dB
        wy  = chebwin(M, at)
        wx  = chebwin(N, at)

    else:
        raise ValueError(f"Unknown window `{win}`")

    # form 2D separable window and apply
    W2d     = np.outer(wy, wx)

    return W2d
    
    
    
  
  
def get_binned_3DBS(bin_edges, start_shape, n_shapes, rdata, Uavg_phy, k2_para, k1_para, file_loc, mask, d_type, Rlzn, rsp_key):
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    for rp in range(start_shape,n_shapes):
        r1 = int(rdata[rp][0])
        r2 = int(rdata[rp][1])
        r3 = int(rdata[rp][2])
        Ntri = np.round(  rdata[rp][3]/(1176*1176)  )
        
        if (r1==0 or r2==0 or r3 ==0): continue

        #data_k, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)
        #filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc, r1, r2, r3, d_type,Rlzn)
        if rsp_key==0:
        	data_k = cylbs_win(r1,r2,r3,file_loc, d_type,Rlzn) 
        elif rsp_key==1: 
        	data_k = cylbs_RSP_win(r1,r2,r3,file_loc, d_type,Rlzn) 
        else: print('wrong RSP key')

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / mf.rc
        k2_perp = ell2 / mf.rc
        k3_perp = ell3 / mf.rc
        k3_para = - k1_para - k2_para

        k1 = np.sqrt(k1_perp**2 + k1_para**2)
        k2 = np.sqrt(k2_perp**2 + k2_para**2)
        k3 = np.sqrt(k3_perp**2 + k3_para**2)

        K_stack = np.stack((k1, k2, k3), axis=0)

        K_sorted  = np.sort(K_stack, axis=0)
        k1_sorted = K_sorted[2]  # largest value
        k2_sorted = K_sorted[1]  # middle value
        k3_sorted = K_sorted[0]  # smallest value

        t  = k2_sorted / k1_sorted
        mu = 0.5 * (t + 1. / t - (k3_sorted / k1_sorted) * (k3_sorted / k2_sorted))

        k1_flat   = k1_sorted.flatten()
        mu_flat   = mu.flatten()
        t_flat    = t.flatten()
        data_k_flat = data_k.flatten()
        mask_flatten= mask.flatten()

        samples = np.column_stack([k1_flat, t_flat, mu_flat])

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat*mask_flatten)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges, weights=mask_flatten)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
      
  
  
 
  
def get_binned_3DBS_wow(bin_edges, start_shape, n_shapes, rdata, Uavg_phy, k2_para, k1_para, file_loc, mask, d_type, Rlzn, rsp_key):  ## without window
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    for rp in range(start_shape,n_shapes):
        r1 = int(rdata[rp][0])
        r2 = int(rdata[rp][1])
        r3 = int(rdata[rp][2])
        Ntri = np.round(  rdata[rp][3]/(1176*1176)  )
        
        if (r1==0 or r2==0 or r3 ==0): continue

        #data_k, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)
        #filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc, r1, r2, r3, d_type,Rlzn)
        if rsp_key==0:
        	data_k = cylbs(r1,r2,r3,file_loc, d_type,Rlzn,False,n_bins=30) 
        elif rsp_key==1: 
        	data_k = cylbs_RSP(r1,r2,r3,file_loc, d_type,Rlzn,False,n_bins=30) 
        else: print('wrong RSP key')

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / mf.rc
        k2_perp = ell2 / mf.rc
        k3_perp = ell3 / mf.rc
        k3_para = - k1_para - k2_para

        k1 = np.sqrt(k1_perp**2 + k1_para**2)
        k2 = np.sqrt(k2_perp**2 + k2_para**2)
        k3 = np.sqrt(k3_perp**2 + k3_para**2)

        K_stack = np.stack((k1, k2, k3), axis=0)

        K_sorted  = np.sort(K_stack, axis=0)
        k1_sorted = K_sorted[2]  # largest value
        k2_sorted = K_sorted[1]  # middle value
        k3_sorted = K_sorted[0]  # smallest value

        t  = k2_sorted / k1_sorted
        mu = 0.5 * (t + 1. / t - (k3_sorted / k1_sorted) * (k3_sorted / k2_sorted))

        k1_flat   = k1_sorted.flatten()
        mu_flat   = mu.flatten()
        t_flat    = t.flatten()
        data_k_flat = data_k.flatten()
        mask_flatten= mask.flatten()

        samples = np.column_stack([k1_flat, t_flat, mu_flat])

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat*mask_flatten)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges, weights=mask_flatten)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
      
  
  
  
  
 
def get_binned_3DBS_WdgAvd(bin_edges, start_shape, n_shapes, rdata, Uavg_phy, k2_para, k1_para, file_loc, mask_key, d_type, Rlzn, rsp_key):
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    
    WF = 3.5121801445642684
    for rp in range(start_shape,n_shapes):
        r1 = int(rdata[rp][0])
        r2 = int(rdata[rp][1])
        r3 = int(rdata[rp][2])
        Ntri = np.round(  rdata[rp][3]/(1176*1176)  )
        
        if (r1==0 or r2==0 or r3 ==0): continue

        #data_k, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)
        #filename = "%s/bsout2_ra%d_rb%d_rc%d_GV_%s_%d_LR" % (file_loc, r1, r2, r3, d_type,Rlzn)
        if rsp_key==0:
        	data_k = cylbs_win(r1,r2,r3,file_loc, d_type,Rlzn) 
        elif rsp_key==1: 
        	data_k = cylbs_RSP_win(r1,r2,r3,file_loc, d_type,Rlzn) 
        else: print('wrong RSP key')

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / mf.rc
        k2_perp = ell2 / mf.rc
        k3_perp = ell3 / mf.rc
        k3_para = - k1_para - k2_para
        
        ############### Preparing Forground avidance mask #################
        klv = np.linspace( -dkll*(Nnc_hf-1), dkll*(Nnc_hf-1), Nnc )
        kl1m, kl2m = np.meshgrid(klv, klv, indexing='ij')
        kl3m  = -kl1m - kl2m
        kl1m_cut = np.where(np.abs(kl1m) < WF*k1_perp, 0, 1)
        kl2m_cut = np.where(np.abs(kl2m) < WF*k2_perp, 0, 1)
        kl3m_cut = np.where(np.abs(kl3m) < WF*k3_perp, 0, 1)

        total_info = kl1m_cut+kl2m_cut+kl3m_cut
	# A1: all 1,2,3 → 1; 0 → 0
        A1 = (total_info > 0).astype(int)
	# A2: 2,3 → 1; 1 → 0; 0 → 0
        A2 = (total_info >= 2).astype(int)
	# A3: 3 → 1; everything else → 0
        A3 = (total_info == 3).astype(int)

        if mask_key == 'A1':
                mask = A1
        elif mask_key == 'A2':
                mask = A2
        elif mask_key == 'A3':
                mask = A3
        else:
                print("Invalid mask_key\n")
       ############## Done #######################################################   

        k1 = np.sqrt(k1_perp**2 + k1_para**2)
        k2 = np.sqrt(k2_perp**2 + k2_para**2)
        k3 = np.sqrt(k3_perp**2 + k3_para**2)

        K_stack = np.stack((k1, k2, k3), axis=0)

        K_sorted  = np.sort(K_stack, axis=0)
        k1_sorted = K_sorted[2]  # largest value
        k2_sorted = K_sorted[1]  # middle value
        k3_sorted = K_sorted[0]  # smallest value

        t  = k2_sorted / k1_sorted
        mu = 0.5 * (t + 1. / t - (k3_sorted / k1_sorted) * (k3_sorted / k2_sorted))

        k1_flat   = k1_sorted.flatten()
        mu_flat   = mu.flatten()
        t_flat    = t.flatten()
        data_k_flat = data_k.flatten()
        mask_flatten= mask.flatten()

        samples = np.column_stack([k1_flat, t_flat, mu_flat])

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat*mask_flatten)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges, weights=mask_flatten)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
      
      
