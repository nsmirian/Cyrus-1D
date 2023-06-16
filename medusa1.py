###########
''' Nun-average free electron laser code ... we called Cyruss 1D
... programmer:  Najmeh Mirian'''
''' The code is similar to MEDUSA
for more information contact nsmirian64@gmail.com '''

''' 16 june 2023 '''
''' the code translated from fortran to python, still require some work  ''''


###########
from input import *
from parameters import *
import numpy as np
from Cyress_fuctions import *
from scipy import integrate
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')
s=-1
step=20
Nstep=int(num_step/step)
power=np.zeros((Nstep+1,max_h, nslice))
#
#--------------------------------------------#
#                 Body of the Program        #
#--------------------------------------------#

#______________________________________________
#				Pierce Parameter /roh
ro1= 1/gamma*omegaw01 * 0.25 * np.max(omegab  )** 2.0 /3.0

#_________________________________________________
#		Diffrection effect
diffracton_effect=diffractioneffect(ro1)
for i in range(nslice): omegab[i]=omegab[i]*diffracton_effect


#call write_data
#print*, diffracton_effect, omegab(nslice/2)
#omegab=omegab*(diffracton_effect**(3./2.))
#print*, diffracton_effect, omegab
#___________________________________________________
# testing resonance condition

lambdarp=(lambda_w/(2*gamma**2))*(1+(omegaw01**2)/2)
lambdarh=(lambda_w/(2*gamma**2))*(1+omegaw01**2)
#___________________________________________________________
#-----------------------------------------------------------
#-----------------------------------------------------------

      #----------------------------------------#
      #              Runge-Kutta               #
      # Yn(i+1)=Yn(i)+h/6*(Kn1+2Kn2+2Kn3+Kn4)  #
      #----------------------------------------#
#------assigning initial values to the radiation and electron beam parameters    #
U00=filling(1)   # matris 4*macroparticles*nslice
#-----------# adding shot noise #------------------------------#
if (shotnoise_choice==1): sai0_noise=add_shot_noise(actual_ne1, U00[0,:,:])
U00[0,:,:]=sai0_noise
#------------# adding dispersion section ----------------------#
if (Disp_sec==1):
    for j in range(neu*nep):
        relco=np.sqrt(1+U00[1,j,:]**2+U00[2,j,:]**2+U00[3,j,:]**2)
        U00[0,j,:]=(2*np.pi/lambdar)* R56*(relco[:]-gamma)/gamma+U00[0,j,:]
#-----------------------------------------------------------
U00[0,:,:]=harmconvert* U00[0,:,:]

for i in range(max_h):
    dela00[0,i, : ]=dela0[:]*(i+1)*0.1  # real part
    dela00[1,i, : ]=dela0[:]*(i+1)*0.1  #imaginery parameters
for ns in range(nslice):
    z=0
    s=-1
    for i in range(num_step):
        dela=dela00[:,:, ns]
        U=U00[:,:,ns]

        #1
        m, mu=rungkutta(z,U, dela , omegab[ns],ns)
        U=U00[:,:,ns]+mu*0.5
        dela=dela00[:,:, ns]+m*0.5
        #print(m)
        #2
        p, pu=rungkutta(z+0.5*h,U, dela,omegab[ns],ns )
        U=U00[:,:,ns]+pu*0.5
        dela=dela00[:,:, ns]+p*0.5
        #3
        q, qu=rungkutta(z+0.5*h,U, dela ,omegab[ns],ns)
        U=U00[:,:,ns]+qu
        dela=dela00[:,:, ns]+q
        #4
        r, ru=rungkutta(z+h,U, dela,omegab[ns] ,ns)
        U=U00[:,:,ns]+( mu+ 2*pu + 2*qu + ru )/6.0
        dela=dela00[:,:, ns]+( m+ 2*p+ 2*q + r )/6.0

        U00[:,:,ns]=U
        #slippage added later
        dela00[:,:,ns]=dela

        z=z+h
        if (np.mod(i,step)==0):
            s=s+1
            print(s)
            for hn in range(max_h):
                power[s, hn, ns]=powerco * (hn* omega )**2.0 /1.0e7*(dela[0, hn]**2+dela[1,hn]**2)

print(power)
