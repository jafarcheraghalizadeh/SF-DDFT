import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
colors =colors = [
    'blue', 'orange', 'forestgreen', 'red', 'purple', 'olive', 'deeppink',
    'gray', 'cyan', 'brown', 'black', 'pink', 'gold', 'lightblue', 'darkviolet',
    'darkorange', 'yellow', 'teal', 'navy','lime', 'turquoise'
]
markers = ['s', 'v','1','o','>','*','<','*','D','H','P','X','|','_','1','2','3','4','8','s']
Nn=1024
Lx, Ly, Lz = 6.6, 6.6, 6.6  # Physical size of the domain along x and y
Nx, Ny, Nz = 18, 18, 18  # Grid size
def compute_gradient(phi, dx, dy, dz):
    """Compute the gradient of phi using finite differences in 3D."""
    dphi_dx = np.gradient(phi, dx, axis=0)
    dphi_dy = np.gradient(phi, dy, axis=1)
    dphi_dz = np.gradient(phi, dz, axis=2)
    return dphi_dx, dphi_dy, dphi_dz

# Compute the expectation value
def compute_tensor(phi, dx, dy, dz):
    """Compute sqrt(<u_i u_j - u^2 delta_ij>) for a 3D field."""
    dphi_dx, dphi_dy, dphi_dz = compute_gradient(phi, dx, dy, dz)
    
    u_x = dphi_dx
    u_y = dphi_dy
    u_z = dphi_dz
    u2 = u_x**2 + u_y**2 + u_z**2

    u_xx = u_x * u_x
    u_yy = u_y * u_y
    u_zz = u_z * u_z
    u_xy = u_x * u_y
    u_xz = u_x * u_z
    u_yz = u_y * u_z

    avg_u_xx = np.mean(u_xx)
    avg_u_yy = np.mean(u_yy)
    avg_u_zz = np.mean(u_zz)
    avg_u_xy = np.mean(u_xy)
    avg_u_xz = np.mean(u_xz)
    avg_u_yz = np.mean(u_yz)
    avg_u2 = np.mean(u2)

    delta_ij = np.eye(3)
    tensor = 0.5 * np.array([
        [avg_u_xx  , avg_u_xy, avg_u_xz],
        [avg_u_xy, avg_u_yy  , avg_u_yz],
        [avg_u_xz, avg_u_yz, avg_u_zz  ]
    ])
    
    sqrt_tensor = np.sqrt(np.abs(tensor))
    return sqrt_tensor

shape = (Nx, Ny, Nz)
dx = Lx / Nx
dy = Ly / Ny
dz = Lz / Nz
files = [f"sample{i}.dat" for i in range(1, 21)]
if __name__ == "__main__":
    j = 0
    jj=0
    order11= np.empty((0, (Nn))) 
    order22= np.empty((0, (Nn))) 
    order33= np.empty((0, (Nn))) 
    for name in files:

        ang = []
        print(name)

        with open(name, 'r') as file:
            order1 = []
            order2 = []
            order3 = []
            for i in range(0, Nn):
                data = np.loadtxt(islice(file, Nx * Ny * Nz), usecols=[0])
                if len(data) < Nx * Ny * Nz:
                    print('Not enough data in', name)
                    break
                
                Z = np.reshape(data, shape)
                sqrt_tensor = compute_tensor(Z, dx, dy, dz)
                eigenvalues, eigenvectors = np.linalg.eigh(sqrt_tensor)
                ang.append((0.5 * np.arctan2(2 * sqrt_tensor[0, 1], (sqrt_tensor[1, 1] - sqrt_tensor[0, 0]))) * 180.0 / np.pi)
                order1.append((eigenvalues[1]-eigenvalues[0]) / (eigenvalues[2]))
                order2.append((eigenvalues[2]-eigenvalues[1]) / (eigenvalues[2]))
                order3.append((eigenvalues[0]) / (eigenvalues[2]))
        order11=np.vstack((order11,order1))
        order22=np.vstack((order22,order2))
        order33=np.vstack((order33,order3))
        xx=np.linspace(0.0025,3,Nn)
        if (j+1)%2==0:
             
             axes[0].plot(xx,order1, marker=markers[jj], color=colors[jj],markevery=50, lw=1.0, ms=0, alpha=1)
             axes[1].plot(xx,order2, marker=markers[jj],color=colors[jj], markevery=50, lw=1.0, ms=0, alpha=1)
             axes[2].plot(xx,order3, marker=markers[jj],color=colors[jj], markevery=50, lw=1.0, ms=0, alpha=1)
             axes[0].set_ylabel(r'$\frac{(\lambda_2-\lambda_1)}{\lambda_3}$', fontsize=20)
             axes[1].set_ylabel(r'$\frac{(\lambda_3-\lambda_2)}{\lambda_3}$', fontsize=20)
             axes[2].set_ylabel(r'$\frac{(\lambda_1)}{\lambda_3}$', fontsize=20)
             jj+=1
        j+=1
    
    axes[0].plot(xx,np.mean(order11,axis=0),ls='--',lw=2.0,label='Average',color='k')
    axes[1].plot(xx,np.mean(order22,axis=0),ls='--',lw=2.0,color='k')
    axes[2].plot(xx,np.mean(order33,axis=0),ls='--',lw=2.0,color='k')
    axes[0].legend()
    axes[0].set_title('linearity')
    axes[1].set_title('planarity')
    axes[2].set_title('sphericity')
    for ax in axes[0:3]:
        ax.set_ylim(0, 1.1)
        ax.set_xlabel(r'$t$', fontsize=20)
        
        ax.tick_params(which='major', direction='in', length=8, width=1, grid_color='r', labelsize=18, right=True, top=True)
        ax.tick_params(which='minor', direction='in', length=3, width=1, grid_color='r', labelsize=15, right=True)
        ax.grid(color='k', alpha=0.3) 
    plt.tight_layout()
    plt.savefig('lambda_alpha_3D.jpg', dpi=800)
    plt.show()
    
