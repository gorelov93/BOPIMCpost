#!/usr/bin/python
import sys
import os
import numpy

def openAver(Fname):
    with open(Fname) as fn:
        a=numpy.genfromtxt(fn,skip_header=3,usecols=range(0,19*19*12),invalid_raise=False)
    ndel=numpy.where(a[:,1]==1)[0]
    print(ndel)
    a_del=numpy.delete(a,ndel,axis=0)
    return numpy.asarray([numpy.average(a_del,axis=0),numpy.std(a_del,axis=0)/(a_del.shape[0]-1)**0.5])

if __name__ == '__main__':
    b=numpy.zeros((8,2,19*19*12))
    #b[0]=openAver('rmc_exc0.nt.1.optex')
    for i in range(int(sys.argv[2]),int(sys.argv[2])+8):
        print(i)
        b[i-int(sys.argv[2])]=openAver('rmc_exc7.nt.%d.optex' % i)
    numpy.save('optex7_%s_aver' % sys.argv[1],b)
