#    constants                                                                 !
#------------------------------------------------------------------------------!
from input import *
import numpy as np
from Cyress_fuctions import gaussian
q_e=4.8e-10_8               # statcoulomb
m_e=9.1e-28_8               # g
c=3.0e10_8                  # cm/s
amp_to_statamp=3.0e9_8      # amp_to_statamp
watt_to_erg_per_s=1.0e7_8   # watt_to_erg/s

ne=nep*neu
pi=np.pi
rs=2.*pi/h*zsep
seed_slids=sigma_r/(lambdar*zsep)
overlap_point=overlap_point/(lambdar*zsep)
k=lambda_w/lambdar
bunchslide=bunch_l/(lambdar*zsep)
# from here we change the unit on lambdar and lambda_w
lambda_w=lambda_w*100
lambdar=lambdar*100
delta_uz= delta_gamma
num_step=np.int(nwp*2.*pi/h)
k_w=2.0_8*pi/lambda_w
r_b=np.sqrt(r_bx*r_by)
r_a=r_b
u0=np.sqrt(gamma**2-1)
omega=k
nslipage =num_step/rs
omegaw01=Kvalue

#----------------
# intializing
#-------------------
slides=list(range(1,nslice+1))
p_in=list(range(1,nslice+1))
dela0=np.zeros((nslice))
omegab=list(range(1,nslice+1))
actual_ne1=list(range(1,nslice+1))
dela00=np.zeros((2,max_h, nslice))
dela=np.zeros((2,max_h))
U=np.zeros((4,ne))
powerco=(c/(4)*(r_a**2)*(k_w*m_e*c**2/q_e)**2)/2.
for i in range(nslice):
    p_in[i]=p00in*gaussian(i+1, overlap_point, seed_slids)
    dela0[i]=q_e/(m_e*c**2)*np.sqrt( 2*pi*c*p_in[i]*watt_to_erg_per_s/
                (omega**2*(r_a**2)*pi) )*np.sqrt(2.)
    omegab[i]=(1/(c*k_w*r_b))*np.sqrt(4*gamma*q_e*i_b1*
                (gaussian(i+1,4*bunchslide , bunchslide))*
                    amp_to_statamp/(c*m_e*np.sqrt(gamma**2-1)))

    actual_ne1[i]=i_b1*(gaussian(i+1,4*bunchslide , bunchslide))*\
                    amp_to_statamp*lambdar/( q_e*c*(1-1/gamma**2) )
