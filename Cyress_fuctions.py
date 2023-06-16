#fuctions

def gaussian(x, mu, sig):
    import numpy as np
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
#===============================================================
#		Diffraction Effect in Dattoli booklet
#  Prof. Datoli function
#====================================
#def diffractioneffect(F, ro1):
    ##mux=lambda_w *lambda/((4.*pi)**2.*betax*emitx*ro1)
    ##muy=lambda_w *lambda/((4.*pi)**2.*betay*emity*ro1)
    #F=( (1.+mux)*(1.+muy) )**(-1./6.)

#return F
#======================================
####  M. Xie function
##=====================================
def diffractioneffect(ro1):
    import numpy as np
    from input import r_by, r_bx, emitx, emity, gamma, delta_gamma, lambda_w, lambdar
    a=list(range(20))
    betay=(r_by)**2/(emity/gamma)
    betax=(r_bx)**2/(emitx/gamma)
    a[1] =0.45
    a[2] =0.56
    a[3] =0.55
    a[4] =1.6
    a[5] =3
    a[6] =2
    a[7] =0.35
    a[8] =2.9
    a[9] =2.4
    a[10]=51
    a[11]=0.95
    a[12]=3
    a[13]=5.4
    a[14]=0.7
    a[15]=1.9
    a[16]=1140
    a[17]=2.2
    a[18]=2.9
    a[19]=3.2
    beta=np.sqrt(betax*betay)
    emitt=np.sqrt(emitx*emity)/gamma
    Ld=lambda_w/(4*np.pi*np.sqrt(3.)*ro1)
    Lr=4*np.pi*r_bx**2/lambdar
    md=Ld/Lr
    me=Ld/(beta) *(4*np.pi*emitt)/lambdar
    mg=4*np.pi*(Ld/lambda_w)*( delta_gamma/gamma)
    mu=a[1]*md**a[2] + a[3]*me**a[4]  + a[5]*mg**a[6]+\
        a[7]*me**a[8]*mg**a[9] + a[10]*md**a[11] * mg**a[12]+\
        a[13]*md**a[14]*me**a[15]+\
        a[16]*md**a[17]* mg**a[19]
    return 1./( 1+mu)




      #--------------------------------------------------------------------------------#
      #            filled the phase space by using Gauss Quadrature for slow beam      #
      #--------------------------------------------------------------------------------#

def filling( self):
    import numpy as np
    import input as inpt
    macroparticles=inpt.nep*inpt.neu
    u00=np.zeros((4,macroparticles, inpt.nslice))
    # Gauss-Legendre (default interval is [-1, 1])
    xp, wp = np.polynomial.legendre.leggauss(inpt.nep)
    xu, wu = np.polynomial.legendre.leggauss(inpt.neu)
    sai_width=inpt.bunchingco*np.pi
    delta_uz= inpt.delta_gamma* (inpt.gamma /np.sqrt(inpt.gamma**2-1))
    u0=np.sqrt(inpt.gamma**2-1)  # without chirp , all slides has same energy and energy spread
    #a=u0-delta_uz
    #b=u0+delta_uz
    a=0.0
    b=sai_width
    for j in range(inpt.nep):
    	for i in range(inpt.neu):
            u00[0, inpt.neu * (j) + i, :]=(xp[j]+1)*sai_width/2 #point_sai[j]
            #u00[1,neu * (j-1) + i, inpt.nslices]=0.0
            #u00[2,neu * (j-1) + i, inpt.nslices]=0.0
            u00[3,inpt.neu * (j) + i,:]=(xu[i])* delta_uz/2+u0#point_u(i)
    return u00




#if(pf==1)then
#do i=1, ne
#READ(UNIT=101, FMT="(4(F10.4,1x))") sai00(i), ux00(i) , uy00(i), uz00(i)
#ux00(i)=0.0

#sai00(i)=harmconvert* sai00(i)
#write(2000, FMT="(2(F10.4,1x))") sai00(i),  uz00(i)
#enddo

def add_shot_noise(actual_ne1,sai0):
    import numpy as np
    import input as inpt
    from random import random
    ne=inpt.nep*inpt.neu
    sai0_noise=np.zeros((ne, inpt.nslice))

 #call random_number(harvest)
    for j in range(inpt.nslice):
        for jj in range(ne):
            ii=1  # harmonic number
            shotnoise=0.0+\
	                1/( (2*ii-1) * np.sqrt(actual_ne1[j]) ) *\
					np.sin( (2*ii-1) * (sai0[jj,j]-2*np.pi*random() ))
            sai0_noise[jj,j]=sai0[jj, j]+shotnoise
    return sai0_noise

#_____________________________________________________________
#_____________________________________________________________
#______________________________________________________________

def rungkutta(z,U, dela, omegb, ns ):
    from input import max_h, h
    from parameters import ne, k, omega
    import numpy as np
    t=np.zeros((2,max_h))
    U2=np.zeros((4,ne))

#----------------------------------------------------------#
# This part produces the current matrixes i.e., s1 and s2  #
#----------------------------------------------------------#
    for hn in range(max_h):
        ave1, ave2=average( U, hn, ns)
        s1=(omegb**2/((hn+1)*k))*( ave1)
        s2=(omegb**2/((hn+1)*k))*(-ave2)
        t[0,hn]=h*( s1 )
        t[1,hn]=h*( s2 )

    for j in range(ne):
	##Dynamics Eq
        relco=np.sqrt(1+U[1,j]**2+U[2,j]**2+U[3,j]**2)
        radittion_effect1=0.0
        radittion_effect2=0.0
        for hn in range(max_h):
            radittion_effect1=((hn+1)*k)*(1-relco/U[3,j])*\
                (dela[0,hn]*np.cos(U[0,j])-dela[1,hn]*np.sin(U[0,j]) )
            radittion_effect2=-((hn+1)*k)*(U[1,j]/U[3,j])*\
                ( dela[0,hn]*np.cos(U[0,j])-dela[1,hn]*np.sin(U[0,j]) )
        U2[1,j]=h*( omegaw(z)*np.sin(z)+radittion_effect1)
        U2[3,j]=h*( -omegaw(z)*np.sin(z)*(U[1,j]/U[3,j])+radittion_effect2 )
        U2[0,j]=h*(k-omega*relco/U[3,j])

    return  t, U2
#--------------------------------------------------------------------------------#
#                    Calculation averages in current for beam               #
#--------------------------------------------------------------------------------#
#                     uz00,sai00, sai,ux,uz,ave1,ave2
                    #V[4,::], V[0,::] U[0,::]
#average(  U, hn, ns  )
def average(U1,hn,ns):
    from input import nep, neu, gamma,delta_gamma, bunchingco
    import numpy as np
    V=filling( 1)
    #uz00=V[3,:,:]
    #t10=V[0,:,:]*hn
    # t1=sai, t2=ux, t3=uz
    delta_uz= 2*delta_gamma* (gamma /np.sqrt(gamma**2-1))
    sai_width=bunchingco*np.pi
    xp, wp = np.polynomial.legendre.leggauss(nep)
    xu, wu = np.polynomial.legendre.leggauss(neu)
    sai_width=bunchingco*np.pi
    delta_uz= delta_gamma* (gamma /np.sqrt(gamma**2-1))

    s1=0
    s2=0

    for i in range(nep):
        for j in range(neu):
            weight1=wp[i]*wu[j]

            s1=s1+weight1*((sigma(V[0,neu*i+j, ns]*(hn+1))*\
                in_dis(V[3,neu*i+j, ns])*U1[1,neu*i+j]*\
                np.cos(U1[0,neu*i+j]*(hn+1)))/np.abs(U1[3,neu*i+j]))

            s2=s2+weight1*((sigma(V[0,neu*i+j,ns]*(hn+1)) *\
                in_dis(V[3,neu*i+j, ns])*U1[1,(neu*i+j)] *\
                np.sin(U1[0,(neu*i+j)]*(hn+1)))/np.abs(U1[3,neu*i+j]))


    ave1=np.sqrt(gamma**2-1)/gamma *s1/(sai_width*2.*delta_uz)
    ave2=np.sqrt(gamma**2-1)/gamma *s2/(sai_width*2.*delta_uz)
    #print(ave1, ave2)
    return ave1,ave2


#****************************************************************************************
#***************************************************************************************
      #-------------------------------------------------------------------------#
      #      Definition of Gousian distribution functions  in entry times	    #
      # Note that for thermal beam it is used				      	            #
      #-------------------------------------------------------------------------#
def in_dis(v0):
#.........NOTE: ****  FOR NONTHERMAL BEAM, THIS COD CAN'T BE UESED & IT SHOULD BE MODIFIED  ****............#
#....................for example, betaz0 is omitted in equations and is added to gousian dis. func.........#
#...........................so it should be correct for nonthermal state...................................#
    #u0=np.sqrt(gamma**2-1)
    in_disf=1.0#(dsqrt(2/pi))*(1/delta_uz) * (v0/dsqrt(1 + v0**2))*exp(-2*((v0-u0)/delta_uz)**2)
#in_disf=(dsqrt(2/pi))*(1/delta_uz) * exp(-2*((v0-u0)/delta_uz)**2)  # distribution for test_integral
    #in_disf=(v0/np.sqrt(1.0+v0**2.))/(2*delta_uz)
#(dsqrt(2/pi))*(1/(delta_uz)) * (v0/dsqrt(1 + v0**2))*exp(-2*((v0-u0)/(delta_uz))**2)
#in_disf=(dsqrt(2/pi))*(1/delta_uz) * exp(-2*((v0-u0)/delta_uz)**2)  #normalized distribution for test_integral

    return in_disf
#***********************************************************************************************************
#***********************************************************************************************************
#--------------------------------------------------------------------------------#
#      Definition of functions that describes distribution in entry times        #
#     Note that for bunched beam it is not uniform                               #
#--------------------------------------------------------------------------------#
def sigma (tt):
    from input import bunchingco
    sai_width=bunchingco
    sigmaf=(2./sai_width)
    return sigmaf
#------------------------------------------#
#    taperred magnetic field amplitute     #
#------------------------------------------#
#                z
def omegaw (t1):
    from input import Kvalue, n_w
    import numpy as np
    if (t1<=n_w*2*np.pi): omegawf=Kvalue*(np.sin(t1/(4.*n_w)))**2
    else: omegawf=Kvalue
    return omegawf
