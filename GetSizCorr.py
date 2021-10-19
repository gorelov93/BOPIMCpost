import numpy
#import routines
import matplotlib
#matplotlib.use('svg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

eV=27.2114

"vmad from Markus' code"
vmad=2.64209138414202/10.267181959

#vmadC2C=[vmadC2CP300,vmadC2CP350,vmadC2CP400,vmadC2CP250]

"function that fits S(k) and calculates size corrections using v_mad"

def Quad(x,a,b):
    return a+b*x**2

def line(x,a,b):
    return a+x*b

"this is how the size corrections are done"
#I compute corrections with sigma and without and output the ratio to check the
#accuracy

def getSizeCorr(SofK,Kgrid,cutoff,vmad,Sigma=[None]):
    # SofK=Ne*See(k,Ne)-Np*See(k,Np)
    # Sigma=Var(SofK)
    ind=numpy.where(Kgrid**2<cutoff)[0]
    sind=0
    if(Sigma[0]!=None): #Sigma is optional
        sind=ind
    sigmafit=curve_fit(line,Kgrid[ind]**2,SofK[ind],sigma=Sigma[sind])
    nosigmafit=curve_fit(line,Kgrid[ind]**2,SofK[ind])

    x=numpy.arange(0,cutoff,0.1)
    plt.figure()
    plt.errorbar(Kgrid[ind]**2,SofK[ind],yerr=Sigma[sind],fmt='o')
    plt.plot(x**2,line(x**2,*sigmafit[0]),'-')
    plt.plot(x**2,line(x**2,*nosigmafit[0]),'-')
    plt.xlim(0,cutoff+0.5)
    plt.legend(['sigma','nosigma'])
    plt.show(block=True)
    return numpy.asarray((sigmafit[0][0]*vmad/2., nosigmafit[0][0]*vmad/2.))


"BragPicks corrected S(K) for C2C"
M=32
name='rmc'
nconf=1
nelec=4
out_name='rmc_corr'
Scordata=numpy.stack((numpy.stack(((numpy.loadtxt(('c%d/%s%d.rhok2twist' % (c,name,k)),skiprows=2,
                                                               usecols=(range(10,M*4+8,4))) -
                                    numpy.loadtxt(('c%d/%s%d.rhok2twist' % (c,name,k)),skiprows=2,
                                                               usecols=(range(11,M*4+8,4))))[:150]
                              for k in range(-nelec,nelec+1,2))
                         )for c in range(1,nconf+1)))
#KgridP300 = numpy.loadtxt('C2CP300/RMC/rmc0.rhok2twist',skiprows=2,usecols=(0))
"K grid "
Kgrid=numpy.stack((numpy.stack(((numpy.loadtxt(('c%d/%s%d.rhok2twist' % (c,name,k)),skiprows=2,
                                                               usecols=(0)))[:150]
                              for k in range(-nelec,nelec+1,2))
                         )for c in range(1,nconf+1)))

"get correction "
CorrAv=numpy.zeros((nelec+1,2))
for n in range(nelec+1):              
    CorrAv[n]=getSizeCorr((numpy.average(Scordata[:,n],axis=(0,2))-numpy.average(Scordata[:,int(nelec/2)],axis=(0,2))),
                           Kgrid[0,n],cutoff=1.5,vmad=vmad,
                           Sigma=numpy.sqrt((numpy.var(Scordata[:,n],axis=(0,2))+numpy.var(Scordata[:,int(nelec/2)],
                               axis=(0,2)))/(M-1.)))

numpy.savetxt('%s.out' % out_name,CorrAv[:,0])
