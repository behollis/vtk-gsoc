from numpy import genfromtxt

PREFIX_U = 'u_r_'
PREFIX_V = 'v_r_'
PREFIX_RHO = 'rho_r_'

PATH = '/home/behollis/DATA/in/ts00050/'

'''
Lock-exchange data parameters.

> Nx: 128
> Ny: 128
> T: 15  (non-dimensional)
> dt: 4.0000e-04
> nu: 2.5000e-04
> kappa: 2.5000e-04
> S: 15 or 20  (total number of modes)
> PlotIntrvl: 50  (number of time-step/interval between plot times)
> var: 0.0500
> MC: 1000
> PrRa: 1
> The realizations are sample path evolution of the flow: I believe that
> there are 1000 of them in what we give you.
'''

def readMember(mem):
    ''' Read numpy matrix from disk and return tuple of realization data. '''
    filename_u = PREFIX_U + str(mem) + '.txt'
    filename_v = PREFIX_V + str(mem) + '.txt'
    #filename_rho = PREFIX_RHO + str(mem) + '.txt'
    
    u = genfromtxt(PATH + filename_u)
    u = -1.0*u[0:127,0:127]
    
    v = genfromtxt(PATH + filename_v)
    #rho = genfromtxt(PATH + filename_rho)
    
    return (u, v)# , rho)
    
if __name__ == '__main__':
    ''' Reads one ensemble member, (u,v) and density. '''
    mem01 = readMember(10)
    