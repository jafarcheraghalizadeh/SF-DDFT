import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fftn, fftshift, fftfreq
import matplotlib

# Plotting settings
matplotlib.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# Parameters
Nx, Ny, Nz = 16, 16, 16  # 3D grid size
sample_size = Nx * Ny * Nz
grid_shape = (Nx, Ny, Nz)
sigma = 0.1
sizex = 6.6
q_int = 1
samples_to_try = [0]  # You can set this to any list of indices you want

files = [  'last_sample3.dat', 'last_sample4.dat', 'last_sample5.dat',
         'last_sample6.dat', 'last_sample7.dat', 'last_sample9.dat', 'last_sample10.dat',
         'last_sample11.dat', 'last_sample12.dat', 'last_sample13.dat', 'last_sample14.dat', 'last_sample15.dat',
         'last_sample16.dat', 'last_sample18.dat', 'last_sample20.dat']

def calculate_S0_3D(phi_A, phi_B):
    phi_A_squared = np.mean(phi_A ** 2)
    phi_B_squared = np.mean(phi_B ** 2)
    cross_correlation = np.mean(phi_A * phi_B)
    numerator = phi_A_squared + phi_B_squared + 2 * cross_correlation
    denominator = phi_A_squared + phi_B_squared
    S0 = numerator / denominator
    return S0

def compute_structure_factor(phi_a_3d, phi_b_3d):
    psi = phi_a_3d - phi_b_3d
    psi_q = fftn(psi)
    S_q = np.abs(psi_q)**2
    return S_q, psi

# Frequency grid
kx = fftfreq(Nx, sizex) * 2 * np.pi * Nx
ky = fftfreq(Ny, sizex) * 2 * np.pi * Ny
kz = fftfreq(Nz, sizex) * 2 * np.pi * Nz
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
k = np.sqrt(KX**2 + KY**2 + KZ**2)

# Radial profile calculation
def radial_profile_3d(data, r):
    r_max = np.max(r) * 1.4
    num_bins = int(np.ceil(r_max)) + 1
    bins = np.linspace(0, r_max, num_bins)
    radial_sum = np.zeros(num_bins - 1)
    count = np.zeros(num_bins - 1)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                radius = int(np.round(r[i, j, k]))
                if radius < num_bins:
                    radial_sum[radius] += data[i, j, k]
                    count[radius] += 1

    radial_mean = np.divide(radial_sum, count, out=np.zeros_like(radial_sum), where=(count != 0))
    valid_bins = count > 2
    return bins[:-1][valid_bins], radial_mean[valid_bins]

# Loop over different sample indices
for sample_index in samples_to_try:
    print(f"\nProcessing sample index: {sample_index}")
    all_structure_factors = []
    all_variances = []

    for file in files:
    
        data = np.loadtxt(file)
   
        start = sample_index * sample_size
        end = (sample_index + 1) * sample_size

        if end > data.shape[0]:
            print(f"Skipping file {file} (not enough samples)")
            continue

        phi_a = data[start:end, 0].reshape(grid_shape)
        phi_b = data[start:end, 1].reshape(grid_shape)

        S_q, psi = compute_structure_factor(phi_a, phi_b)
        S0 = calculate_S0_3D(phi_a, phi_b)
        all_structure_factors.append(S_q)

        psi_mean = np.mean(psi)
        psi_variance = np.mean(psi**2) - psi_mean**2
        all_variances.append(S0)

    if not all_structure_factors:
        print(f"No valid samples found for sample index {sample_index}")
        continue

    S_q_avg = np.mean(all_structure_factors, axis=0)
    psi_variance_avg = np.mean(all_variances)
    print(f"Average Variance of psi (S0): {psi_variance_avg:.4f}")

    # Compute radial profile
    k_profile, P_k_radial = radial_profile_3d(S_q_avg, k * q_int)

    # Save the result
    out_filename = f"radial_profile_sample{sample_index}.dat"
    np.savetxt(out_filename, np.column_stack((k_profile/q_int, P_k_radial/(2 * 0.8)**2 / (Nx * Ny * Nz))))
    print(f"Saved radial profile to {out_filename}")
