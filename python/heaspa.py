import heasp,os,pyfits
from numpy import *

# always clobber!

class ARF:
    def __init__(self,e1,e2,effa):
        self.e1=e1
        self.e2=e2
        self.effa=effa
    
    def write(self,fn):
        arf=heasp.arf()
        nchan=len(self.e1)
        arf.initChannels(nchan)
        print arf.disp()
        for i in range(nchan):
            arf.LowEnergy[i]=self.e1[i]
            arf.HighEnergy[i]=self.e2[i]
            arf.EffArea[i]=self.effa[i]

        if os.path.exists(fn):
            os.remove(fn)
        arf.write(fn)

class PHA:
    def __init__(self,pha,staterr,exposure,syserr=None,datatype="COUNTS",response=None):
        self.pha=pha
        self.staterr=staterr
        self.exposure=exposure
        self.datatype=datatype
        self.response=response
    
    def write(self,fn):
        sp=heasp.pha()
        sp.initChannels(len(self.pha))
        sp.Exposure=self.exposure
        sp.Datatype=self.datatype

        print dir(sp)    

        print sp.disp()
        for i in range(len(self.pha)):
            sp.Pha[i]=self.pha[i]
            sp.StatError[i]=self.staterr[i]
        if os.path.exists(fn):
            os.remove(fn)
        sp.write(fn)
        
        if self.response is not None:
            f=pyfits.open(fn)
            f[1].header['RESPFILE']=self.response
            f.writeto(fn,clobber=True)

class RMF:
    def __init__(self,ce1,ce2,e1,e2,matrix):
        self.ce1=ce1
        self.ce2=ce2
        self.e1=e1
        self.e2=e2
        self.matrix=matrix
    
    def write(self,fn):
        rmf=heasp.rmf()
        rmf.FirstChannel=0
        rmf.AreaScaling=1
        rmf.ResponseThreshold=0
        rmf.RMFVersion="1.3.0"

        print "writing rmf with channels",len(self.ce1)
        for i in range(len(self.ce1)):
            print i,self.ce1[i],self.ce2[i]
            rmf.ChannelLowEnergy.push_back(float(self.ce1[i]))
            rmf.ChannelHighEnergy.push_back(float(self.ce2[i]))

        print "writing rmf with rows",len(self.e1)
        for i in range(len(self.e1)):
            y=map(float,self.matrix[i,:])

            print "rmf row",self.e1[i],self.e2[i],len(y)
            rmf.addRow(y,float(self.e1[i]),float(self.e2[i]))
            

        if os.path.exists(fn):
            os.remove(fn)
        rmf.write(fn)

