import numpy as np
import matplotlib.pyplot as plt

LL_sim = 4. # Mpc
LL_perp = 4.
LL_para = 0.68
theta_0=0.6*24.68
theta_0*=(np.pi/180)
fact = np.pi*theta_0*theta_0/3

Nchan = 768#192
nc10  = 383#95
nclen = 2*nc10+1
ncsize = nclen*nclen
uf_nc  =  552#144
n_nancut = 8

nu0=154.255e6
k_B=1.38e3
c=3e8
rjfac = 2.*k_B*nu0*nu0/(c*c)/1e3; ## mK/Jy
mabs_norm = fact*(rjfac**3)

rc = 9210.11291939594
rc_prime = 17. ## Mpc/MHz

dnu_mwa = 0.04  ## in MHz
dnu_sim = 4.*LL_para/rc_prime

du_phy = 1.0254384902642049
########## Rings ##########################

def define_rings():
    intervals = [
        [0, 50, 5],
        [50, 100, 10],
        [100, 200, 20],
        [200, 270, 35],
    ]
    num_intervals = len(intervals)

    # Parameters
    Nbin = 22  # Define the number of bins as per your requirement
    Umax = 270  # Maximum U value (adjust as needed)
    du = 1.0254384902642049
    # Initialize arrays
    ULow = np.zeros(Nbin, dtype=float)
    UHigh = np.zeros(Nbin, dtype=float)
    Uavg = np.zeros(Nbin, dtype=float)

    # Fill ULow and UHigh
    index = 0
    for i in range(num_intervals):
        start, end, step = intervals[i]
        for j in range(start, end, step):
            if index < Nbin:
                ULow[index] = j
                if index > 0:
                    UHigh[index - 1] = ULow[index]
                index += 1

    # Set boundary conditions
    ULow[0] = 1.0
    UHigh[Nbin - 1] = 1.0 * Umax

    # Compute Uavg
    for i in range(Nbin):
        if UHigh[i] > ULow[i]:  # Ensure valid range for computation
            Uavg[i] = (2.0 / 3.0) * (UHigh[i]**3 - ULow[i]**3) / (UHigh[i]**2 - ULow[i]**2)
    #     print(f"\t{ULow[i]:.6f}\t{UHigh[i]:.6f}\t{Uavg[i]:.6f}")
    Uavg_phy = Uavg
    return Uavg_phy
    
    
  
def CalCR(Uavg, Nr=22):
    valid_combinations = []
    rcount = 0
    for r1 in range(Nr):
        for r2 in range(r1 + 1):  # r2 goes from 0 to r1 (inclusive)
            tt = Uavg[r2] / Uavg[r1]
            if 0.5 <= tt <= 1.0:
                for r3 in range(r2 + 1):  # r3 goes from 0 to r2 (inclusive)
                    mu = 0.5 * (tt + 1.0/tt - (Uavg[r3] / Uavg[r1]) * (Uavg[r3] / Uavg[r2]))
                    if 0.5/tt <= mu <= 1.0:
                        valid_combinations.append((r1, r2, r3))
                        rcount += 1
    return valid_combinations, rcount  
    



import os
def read_rings(res_path):
    """
    Reads ring‐configuration data from two files in res_path:
      - 'rcount'        : a text file containing a single integer count
      - 'closed_rings'  : a binary file containing 4 * count doubles
    
    Returns
    -------
    rcount : int
        Number of ring configurations.
    rdata : numpy.ndarray, shape (rcount, 4)
        The ring data, one row per configuration.
    """
    # -- read the count --
    filename_count = os.path.join(res_path, "rcount")
    try:
        with open(filename_count, 'r') as f:
            rcount = int(f.read().strip())
    except Exception as e:
        raise IOError(f"Error reading ring count from '{filename_count}': {e}")

    # -- read the binary ring data --
    filename_rings = os.path.join(res_path, "closed_rings")
    try:
        # read exactly 4*rcount doubles
        raw = np.fromfile(filename_rings, dtype=np.double, count=rcount * 4)
    except Exception as e:
        raise IOError(f"Error reading ring data from '{filename_rings}': {e}")

    if raw.size != rcount * 4:
        raise ValueError(
            f"Expected {rcount*4} doubles in '{filename_rings}', "
            f"but got {raw.size}"
        )

    # reshape into (rcount, 4) for easy indexing
    rdata = raw.reshape(rcount, 4)

    print(f"Successfully read {rcount} ring configurations.")
    return rcount, rdata









def kp(Uavg_phy,ring):
    ell = 2*np.pi*Uavg_phy[ring]
    return ell/rc 
    
    
def kllH(ring):
    res = rc*kp(Uavg_phy,ring)/(rc_prime*nu0*1e-6)
    return res   
  ############  Read data  ###############
  
def read_bsout(filename):
    with open(filename, 'rb') as f:
        # Read the first header value: Ntri (double)
        Ntri = np.fromfile(f, dtype=np.float64, count=1)[0]
        # Read the normalization value: Num (double)
        Num  = np.fromfile(f, dtype=np.float64, count=1)[0]
        # Read the bldnu array (ncsize doubles)
        bldnu = np.fromfile(f, dtype=np.float64, count=ncsize)
        # Read the numdnu array (ncsize ints)
        numdnu = np.fromfile(f, dtype=np.int32, count=ncsize)

    # Reshape the arrays into 2D grids if desired
    bldnu = bldnu.reshape((nclen, nclen))
    numdnu = numdnu.reshape((nclen, nclen))
    
    return Ntri, bldnu/Num, numdnu
    
    
   
def read_bsout_float(filename):
    with open(filename, 'rb') as f:
        # Read the first header value: Ntri (double)
        Ntri = np.fromfile(f, dtype=np.float64, count=1)[0]
        # Read the normalization value: Num (double)
        Num  = np.fromfile(f, dtype=np.float64, count=1)[0]
        # Read the bldnu array (ncsize doubles)
        bldnu = np.fromfile(f, dtype=np.float32, count=ncsize)
        # Read the numdnu array (ncsize ints)
        numdnu = np.fromfile(f, dtype=np.int32, count=ncsize)

    # Reshape the arrays into 2D grids if desired
    bldnu = bldnu.reshape((nclen, nclen))
    numdnu = numdnu.reshape((nclen, nclen))
    
    return Ntri, bldnu/Num, numdnu
    #return Ntri, bldnu, numdnu*Num ## delete it, just for checking
    
    
def read_bsout_float_v2(filename):
    with open(filename, 'rb') as f:
        # Read the normalization value: Num (double)
        Num  = np.fromfile(f, dtype=np.float64, count=1)[0]
        # Read the bldnu array (ncsize doubles)
        bldnu = np.fromfile(f, dtype=np.float32, count=ncsize)
        # Read the numdnu array (ncsize ints)
        numdnu = np.fromfile(f, dtype=np.float32, count=ncsize)

    # Reshape the arrays into 2D grids if desired
    bldnu = bldnu.reshape((nclen, nclen))
    numdnu = numdnu.reshape((nclen, nclen))
    Ntri=1
    
    return Ntri, bldnu, numdnu


############  Fold data  ###############

def fold_quadrants(arr):
    """Folds the (-x,-y) quadrant into (+x,+y) and (+x,-y) into (-x,+y),
    returning an array of shape (400, 200)."""
    
    center_uv = nc10
    
    upper_half = arr[:center_uv, :]
    
    lower_half = arr[center_uv:, :]
    
    lower_half[1:, :] += upper_half[::-1,::-1]
    
    return lower_half

def fold_quadrants_old(arr):
    """Folds the (-x,-y) quadrant into (+x,+y) and (+x,-y) into (-x,+y),
    returning an array of shape (400, 200)."""
    
    center_uv = nc10
    
    upper_half = arr[:center_uv+1, :]
    
    lower_half = arr[center_uv+1:, :]
    
    upper_half[:center_uv, :] += lower_half[::-1,::-1]
    
    return upper_half



def dnu_grid():
    i_indices = np.arange(nclen)
    j_indices = np.arange(nclen)

    i_minus_nc10 = i_indices - nc10
    j_minus_nc10 = j_indices - nc10

    I, J = np.meshgrid(i_minus_nc10, j_minus_nc10)
    
    return J,I


    
 
def cut_half(mabs_T):
    Nnc_hf = nc10//2+1
    return mabs_T[:Nnc_hf,Nnc_hf-1:Nnc_hf+nc10]   
    
   
def cut_nan(arr,n):
	return arr[:arr.shape[0]-n,n:-n]
    
 #######  BS theory #############
 
def input_pk(kll,kperp,pkn):
    ksq = (kll*kll + kperp*kperp)
    return pow(ksq,pkn/2.)

def Bk(k2ll, k3ll, ell1, ell2, ell3, fng=1., pkn=-2.):
    k1ll   = -k2ll-k3ll
    k1perp = ell1/rc
    k2perp = ell2/rc
    k3perp = ell3/rc
    
    pk1 = input_pk(k1ll,k1perp,pkn)
    pk2 = input_pk(k2ll,k2perp,pkn)
    pk3 = input_pk(k3ll,k3perp,pkn)
    
    pfact = pk1*pk2 + pk2*pk3 + pk3*pk1
    
    return 2.*fng*pfact
    
    

pi = np.pi
def pk_func(ell, i, j, L):
    kperp = ell / rc
    tpl_sq = (2 * pi / L)**2
    kpara = tpl_sq * (i*i + j*j + 2 * i * j)
    ktot = kperp**2 + kpara
    return 1.0 / ktot
    
    
def Bk_func(ell1, ell2, ell3, i, j, L, fng):
    pk1 = pk_func(ell1, i, j,  L)
    pk2 = pk_func(ell2, i,  0,  L)
    pk3 = pk_func(ell3, 0,  j,  L)
    res = 2.*fng*(pk1 * pk2 + pk2 * pk3 + pk3 * pk1)
    return res 

def Bk_func_old(ell1, ell2, ell3, i, j, L, fng):
    pk1 = pk_func(ell1, i, 0,  L)
    pk2 = pk_func(ell2, 0,  j,  L)
    pk3 = pk_func(ell3, i,  j,  L)
    res = 2.*fng*(pk1 * pk2 + pk2 * pk3 + pk3 * pk1)
    return res 
    
 
    
    
########  DFTs #################


def Bk2Bdnu_fft(data_k):
    N = data_k.shape[1]
    data_dnu_r =np.fft.fftshift(np.fft.rfft2(np.fft.ifftshift(data_k), axes=(1, 0)), axes=(1))

    rc4  = rc**4
    L    = N * dnu_mwa*rc_prime 
    Lsq  = L**2
    data_dnu_r *= 1.0 / (rc4 * Lsq)
    
    return data_dnu_r.real


def Bdnu2Bk_fft(data_dnu):
    N = data_dnu.shape[1]
    res = np.fft.fftshift(np.fft.irfft2(np.fft.ifftshift(data_dnu, axes=(1)) ,axes=(1,0), s=(N,N)))
 
    rc4 = rc**4
    L    = N * dnu_mwa*rc_prime 
    Lsq = L**2
    
    res *= (rc4 * Lsq)
    return res


def Bk2Bdnu_fft_old(data_k,N):
    data_dnu_r =np.fft.fftshift(np.fft.rfft2(np.fft.ifftshift(data_k), axes=(1, 0)), axes=(1))
    #data_dnu_r =np.fft.fftshift(np.fft.rfft2(np.fft.ifftshift(data_k)), axes=(0))  ## if you use this then transpose the output data
    #data_dnu_c = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(data_k)))  ## it do full 2D fft instead of half plane
    rc4  = rc**4
    #L   = Nchan* BW*rc_prime 
    L    = N * dnu_sim*rc_prime 
    Lsq  = L**2
    data_dnu_r *= 1.0 / (rc4 * Lsq)
    
    return data_dnu_r.real


def Bdnu2Bk_fft_old(data_dnu,N):
    #ss =  N
    res = np.fft.fftshift(np.fft.irfft2(np.fft.ifftshift(data_dnu, axes=(1)) ,axes=(1,0), s=(N,N)))
    #res = np.fft.fftshift(np.fft.irfft2(np.fft.ifftshift(data_dnu.T, axes=(0)) , s=(ss,ss)))  ## if you use this then transpose the input data
    
    rc4 = rc**4
    #L   = Nchan* BW*rc_prime 
    L    = N * dnu_sim*rc_prime 
    Lsq = L**2
    
    res *= (rc4 * Lsq)
    return res
    
  
pi=np.pi
def pk2cdnu_fft(data_k,N):
#     data_dnu_r =np.fft.fftshift(np.fft.fft(np.fft.ifftshift(data_k)))
    data_dnu_r =np.fft.rfft(np.fft.ifftshift(data_k))
    rc2 = rc**2
    L   = N * dnu_sim *rc_prime 
    data_dnu_r *= 1.0 / (rc2 * L)
    
    return data_dnu_r.real  
  

######## B(kll) from specified ring and dnu-space data file ###########
def get_cyl_BS(r1,r2,r3,dnu_file_loc,Rlzn):
    
    filename = "%s/bsout_final_ra%d_rb%d_rc%d_R%d" % (dnu_file_loc, r1, r2, r3, Rlzn)
    Ntri, mabs, numc = read_bsout_float(filename)
    mabs_F    = fold_quadrants(mabs)
    numc_mean = fold_quadrants(numc)

    mabs_T    = mabs_F/mabs_norm/numc_mean

    mabs_final=  cut_nan(mabs_T,n_nancut)  #cut_half(mabs_T)

    res       = Bdnu2Bk_fft(mabs_final)
    
    return res, np.round(Ntri)   
    
    
  ######## 2D Bispectrum from cyl B(kll) ###########  
def get_binned_3DBS(bin_edges, n_shapes, rdata, Uavg_phy, k2_para, k3_para, file_loc, Rlzn):
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    for rp in range(n_shapes):
        r1 = rdata[rp][0]
        r2 = rdata[rp][1]
        r3 = rdata[rp][2]

        data_k, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / rc
        k2_perp = ell2 / rc
        k3_perp = ell3 / rc
        k1_para = - k3_para - k2_para

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

        samples = np.column_stack([k1_flat, t_flat, mu_flat])

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
    
    
 
def get_Ana_3DBS(bin_edges, n_shapes, rdata, Uavg_phy, k2_para, k3_para, file_loc, Rlzn):
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    for rp in range(n_shapes):
        r1 = rdata[rp][0]
        r2 = rdata[rp][1]
        r3 = rdata[rp][2]

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / rc
        k2_perp = ell2 / rc
        k3_perp = ell3 / rc
        k1_para = - k3_para - k2_para
        
        _, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)
        
        #data_k = mf.Bk_func(ell1, ell2, ell3, i_vals, j_vals, L, 1.0)
        data_k = Bk(k2_para, k3_para, ell1, ell2, ell3, 1., -2.)

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

        samples = np.column_stack([k1_flat, t_flat, mu_flat])

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
    
    
    
    

def get_samples(n_shapes, rdata, Uavg_phy, k2_para, k3_para):
    all_samples=[]
    for rp in range(n_shapes):
        r1 = rdata[rp][0]
        r2 = rdata[rp][1]
        r3 = rdata[rp][2]

        ell1 = 2. * np.pi * Uavg_phy[r1]
        ell2 = 2. * np.pi * Uavg_phy[r2]
        ell3 = 2. * np.pi * Uavg_phy[r3]

        k1_perp = ell1 / rc
        k2_perp = ell2 / rc
        k3_perp = ell3 / rc
        k1_para = - k3_para - k2_para

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

        samples = np.column_stack([k1_flat, t_flat, mu_flat])
        
        all_samples.append(samples)


    return np.array(all_samples)
    
    

def get_binned_3DBS_fast(bin_edges, rdata, n_shapes, all_samples, file_loc, Rlzn):
    
    k1_bins, t_bins, mu_bins = bin_edges[0].size - 1, bin_edges[1].size - 1, bin_edges[2].size - 1
    bispec = np.zeros((k1_bins,t_bins,mu_bins))
    num    = np.zeros((k1_bins,t_bins,mu_bins))
    for rp in range(n_shapes):
        r1 = rdata[rp][0]
        r2 = rdata[rp][1]
        r3 = rdata[rp][2]

        data_k, Ntri = get_cyl_BS(r1,r2,r3,file_loc,Rlzn)

        data_k_flat = data_k.flatten()

        samples = all_samples[rp]

        binned_data, edges = np.histogramdd(samples, bins=bin_edges, weights=data_k_flat)
        binned_num, edges = np.histogramdd(samples, bins=bin_edges)

        bispec += (binned_data * Ntri)
        num    += (binned_num * Ntri)

    return bispec, num
    
    
    
######## B(dnu) from B(k) using numerical integration ###########

def Bk2Bdnu_Int(DeltaNu2, DeltaNu3,rc, rcprime,kpar1_array, kpar3_array, BHI_array):

    # Create a 2D mesh for kpar1 and kpar3
    K1, K3 = np.meshgrid(kpar1_array, kpar3_array, indexing='ij')
    
    # Compute the cosine factor: cos(k1*rc'*DeltaNu2 + k3*rc'*DeltaNu3)
    cos_factor = np.cos(rcprime * (K1 * DeltaNu2 + K3 * DeltaNu3))
    
    # The integrand of the double integral
    integrand = BHI_array * cos_factor
    
    # Perform the 2D integration:
    # 1. Integrate along kpar3 with trapezoidal rule
    int_over_kpar3 = np.trapz(integrand, x=kpar3_array, axis=1)
    
    # 2. Integrate the result over kpar1
    integral_value = np.trapz(int_over_kpar3, x=kpar1_array)
    
    # Multiply by the overall prefactor 1/(4*pi^2 * rc^2 * rc'^2)
    prefactor = 1.0 / (4.0 * np.pi**2 * rc**4)
    
    return prefactor * integral_value
    
    
    
def Pk2cl_Int(DeltaNu,rc, rcprime,kpar1_array, PHI_array):

    # Create a 2D mesh for kpar1 and kpar3
    K1 = kpar1_array
    
    # Compute the cosine factor: cos(k1*rc'*DeltaNu2 + k3*rc'*DeltaNu3)
    cos_factor = np.cos(rcprime * K1 * DeltaNu)
    
    # The integrand of the double integral
    integrand = PHI_array * cos_factor
    
    # Perform the 2D integration:
    # 1. Integrate along kpar3 with trapezoidal rule
    integral_value = np.trapz(integrand, x=K1)
    
    
    # Multiply by the overall prefactor 1/(4*pi^2 * rc^2 * rc'^2)
    prefactor = 1.0 / (2.0 * np.pi * rc**2)
    
    return prefactor * integral_value
    
    
    
    
    
    
  
