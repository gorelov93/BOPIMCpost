import numpy
import routines
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

eV=27.2114
"function that does everything"
def GC_QMC(data,Vdata,mbeg,mend,mstep,Theta):
    mu = numpy.arange(mbeg,mend,mstep)
    nn = numpy.arange(-4,5,2)
    N = numpy.zeros((mu.shape[0],Theta))
    E = numpy.zeros((mu.shape[0],Theta))
    VE = numpy.zeros((mu.shape[0],Theta))

    for i,m in enumerate(mu):
        for theta in range(Theta):
            a=True
            for n,ni in enumerate(nn):
                if (data[n+1,theta]-data[n+1-1,theta])/2.<=m and (data[n+1+1,theta]-data[n+1,theta])/2.>m and a:
#                    a = False
                    N[i,theta]=96+ni
                    E[i,theta]=data[n+1,theta]
                    VE[i,theta]=Vdata[n+1,theta]
    N[N == 0.] = numpy.nan
    E[E == 0.] = numpy.nan
    VN = numpy.nanmean(N**2, axis=1)
    N = numpy.nanmean(N, axis=1)
    E = numpy.nanmean(E, axis=1)
    VE = numpy.sqrt(numpy.sum(VE**2, axis=1))/Theta
    VN = numpy.sqrt((VN-N*N)/Theta)
    a=0.
    for i,m in enumerate(mu):
        if i<mu.shape[0]-2:
            if N[i]==96.: #E[i] == E[i+1] == E[i+2] == E[i+3]:
                #print m*27.2114
                a+=mstep
    return N,E,VE,VN, mu,(a)*27.2114

"load data for different P and T"
M=108
nconf=24
dataT200P400=[]
VdataT200P400=[]
for i in range(nconf):
    dataT200P400.append(numpy.stack((routines.AverageAndVar(numpy.loadtxt(('c%d/C2cp400t02np32%d.np.1.excite' % (i+1, k)),
                                                               usecols=(range(0,M*2,2)), skiprows=100))[0]
                              for k in range(-6,7,2))
                            ))
    VdataT200P400.append(numpy.stack((routines.AverageAndVar(numpy.loadtxt(('c%d/C2cp400t02np32%d.np.1.excite' % (i+1, k)),
                                                          usecols=(range(0,M*2,2)), skiprows=100))[1]
                               for k in range(-6,7,2))
                             ))
corrC2CT200P400=numpy.loadtxt('C2CT200P400corr_av.out')
numpy.save('C2CT200P400_corr',numpy.stack((numpy.asarray(dataT200P400)+numpy.transpose(numpy.tile(corrC2CT200P400*2,(24,108,1)),axes=(0,2,1)))))
