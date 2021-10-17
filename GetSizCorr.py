import numpy
#import routines
import matplotlib
#matplotlib.use('svg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

eV=27.2114

"vmad"
vmad=3.0787488096677436/11.120467409715490

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
M=108

Scordata=numpy.stack((numpy.stack(((numpy.loadtxt(('c%d/C2cp400t02np32%d.rhok2twist' % (c,k)),skiprows=2,
                                                               usecols=(range(10,M*4+8,4))) -
                                    numpy.loadtxt(('c%d/C2cp400t02np32%d.rhok2twist' % (c,k)),skiprows=2,
                                                               usecols=(range(11,M*4+8,4))))[:150]
                              for k in range(-6,7,2))
                         )for c in range(1,25)))
#KgridP300 = numpy.loadtxt('C2CP300/RMC/rmc0.rhok2twist',skiprows=2,usecols=(0))
"K grid for C2C"
Kgrid=numpy.stack((numpy.stack(((numpy.loadtxt(('c%d/C2cp400t02np32%d.rhok2twist' % (c,k)),skiprows=2,
                                                               usecols=(0)))[:150]
                              for k in range(-6,7,2))
                         )for c in range(1,25)))

"get correction for C2C"
CorrAv=numpy.zeros((7,2))
for n in range(7):
    CorrAv[n]=getSizeCorr((numpy.average(Scordata[:,n],axis=(0,2))-numpy.average(Scordata[:,3],axis=(0,2))),
                           Kgrid[0,n],cutoff=1.5,vmad=vmad,
                           Sigma=numpy.sqrt((numpy.var(Scordata[:,n],axis=(0,2))+numpy.var(Scordata[:,3],axis=(0,2)))/107.))

numpy.savetxt('C2CT200P400corr_av.out',CorrAvC2C[:,0])
