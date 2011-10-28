###
## Representative volume element maker
###
import matplotlib.pyplot as plt
import numpy as np
import math
rand = np.random.rand
randi = np.random.random_integers


def steglich_format(filename=None):
    """
    2011-Sept-21
    Convert Dr. Steglich's OCD format into LABOTEX's

    Reference file: 'ODF numerisch.txt'
    """
    f = open(filename, 'r')
    contents = f.read()

    ## Assumes that descrete COD is provided by slice of sections
    ## that are perpendicular to phi axis
    blocks = contents.split('Phi1=')
    header = blocks[0]
    planes = blocks[1:]

    axis_p1 = []
    axis_P = []
    axis_p2 = []
    cod = []

    for i in range(len(planes)): #each block of phi=constant plane
        clines = planes[i].split('\n')
        block = clines[1:][:-1:] #tail off
        block = np.array(block)
        dum = []
        for i in range(len(block)): #phi2
            if i!=0 and len(block[i]) > 3: #PHI
                dum.append(
                    map(float,
                        block[i].split()[1:]
                        )
                    ) #remove the first row
                pass
            pass
        dum = np.array(dum) # dum: (phi2, PHI)
        dum = dum.T         # dum: (PHI, phi2)
        # dum = dum[0:]
        dum = dum.tolist()  # make numpy array into list type
        cod.append(dum) # cod: (phi1, PHI, phi2)
        pass
    
    rst = np.zeros((len(cod), len(cod[0]), len(cod[0][0])))
    for i in range(len(cod)): #phi1
        for j in range(len(cod[i])): #PHI
            for k in range(len(cod[i][j])): #phi2
                rst[i][j][k] = cod[i][j][k]
                pass
            pass
        pass
    print 'rst shape:', rst.shape

    ## write this into LABOTEX descrete COD format file
    
    ##  phi1 phi phi2 COD
    ##   0   0    0    0.002
    ##   5   0    0    0.012
    ##   ..
    ## 360   0    0    0.023
    ##   0   5    0    0.100
    ##   5   5    0    0.123
    ##   ..
    ##   0   0    5    0.603

    # permute the rst(phi1, phi, phi2) -> temp(phi, phi2, phi1)
    temp = np.transpose(rst, (1,2,0))
    print 'temp shape:', temp.shape
    fout = open('%s_labo.txt'%filename.split('.')[0], 'w')
    fout.writelines('%s %s %s %s \n'%('PHI1','PHI2','PHI', 'COD'))
    for i in range(len(temp)): #phi
        for j in range(len(temp[i])): #phi2
            for k in range(len(temp[i][j])): #phi1
                fout.writelines(
                    '  %6.2f  %6.2f  %6.2f  %12.7e\n'%(
                        k*5., j*5., i*5., temp[i][j][k]
                        )
                    )
                pass
            pass
        pass
    return rst

def random(phi1=90, phi2=90, phi=90, ngrain=1000, iplot=False, filename=None):
    """
    Random iostropic texture maker
    Bunge convention

    gr consists of phi1, phi, phi2, and the volume fraction in order.
    """
    gr = []
    #coef = math.atan(1.0)/45.
    for i in range(ngrain):
        cp1 = rand() * phi1  #phi1
        cp2 = rand() * phi2  #phi2
        cp = rand()
        cp = math.acos(cp) * 180./ math.pi
        if phi==180:
            if randi(0,1)==0: cp = cp
            else: cp = 180 - cp
        elif phi==90:
            pass
        else:
            #raise IOError, "Unexpected phi range is given"
            pass
        gr.append([cp1, cp, cp2, 1./ngrain])
        pass

    gr = np.array(gr)
    p1 = gr.transpose()[0]
    p = gr.transpose()[1]
    p2 = gr.transpose()[2]    

    if filename!=None and type(filename).__name__=='str':
        FILE = open(filename, 'w')
        FILE.writelines(
            'dummy\ndummy\ndummy\n%s %i\n'%
            ('B', ngrain))
        for i in range(ngrain):
            FILE.writelines(
                "%10.4e %10.4e %10.4e %10.4e\n"%
                (p1[i], p[i], p2[i], 1./ngrain)
                )
            pass
        pass
    
    if iplot==True:
        fig = plt.figure()
        ax1 = fig.add_subplot(311); ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax1.plot(p1, p, '.')
        ax2.plot(p1, p2, '.')
        ax3.plot(p2, p, '.')
        pass
    return gr


def rve_ortho(cod, rve):
    """
    Apply the orthonormal sample symmetry to the given discrete
    crystallographic orientation distribution (COD)

    Orthonormal (phi1 = 90)
    Monoclinic (phi1=180)
    None (Triclinic) (phi1=360)
    """
    from euler import euler

    codt = cod.transpose()
    ## information ------------------
    p1max = max(codt[0])  #phi1
    print 'p1max: %4.1f'%p1max
    # phi1 = codt[0]
    # phi2 = codt[1]
    # phi = cot[2]
    ## ------------------------------

    if p1max==90: ssym="Orth"
    elif p1max==180: ssym="Mono"
    elif p1max==360: ssym="Tric"
    else: raise IOError, "Unexpected maximum phi1 anlge"
    print 'symmetry: %s'%ssym

    new_rve = [ ]
    for igr in range(len(rve)):
        ## Phi1, Phi, Phi2 angles and volume fraction
        p1 = rve[igr][0]; p = rve[igr][1]
        p2 = rve[igr][2]; vf = rve[igr][3]

        ## rotation matrix of the current grain
        amat = euler(p1, p, p2, echo=False)        
        amat_t = amat.transpose()
        amat_new = []
        if ssym=="Orth":
            ## multiplication of the matrix according to the symmetry

            # x-mirror
            oldt = amat_t.copy()
            oldt[1] = oldt[1]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())

            # y-mirror
            oldt = amat_t.copy()
            oldt[0] = oldt[0]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())

            # x and y-mirror
            oldt = amat_t.copy()
            oldt[0] = oldt[0]*-1
            oldt[1] = oldt[1]*-1
            amat_new.append(oldt.transpose())

            nvol = 4
            pass
        
        elif ssym=="Mono":
            # x-mirror (along TD)
            oldt = amat_t.copy()
            oldt[1] = oldt[1]*-1
            oldt[2] = oldt[2]*-1
            amat_new.append(oldt.transpose())
            nvol = 2

            pass
        
        elif ssym=="Tric":
            nvol=1
            #no mirror axis
            pass

        ## assigns the newly multiplied A-matrix to the new_rve
        temp = rve[igr].copy(); temp[3] = vf/nvol
        new_rve.append(temp)
        for i in range(len(amat_new)):
            ph1, ph, ph2 = euler(a=amat_new[i],echo=False)
            new_rve.append([ph1,ph,ph2,vf/nvol])
            pass
        pass
    return np.array(new_rve)


class RVE:
    """
    Generates the random iostropic aggregate
    and combines with the discrete COD.

    It assumes that the give COD is the reduced space
    in that the boundary of the Euler space of the COD
    indicates the sample symmetry

    Arguments:
       ngrain : number of grains in the RVE
       odf : Orientation Distribution File (Labotex format)
       cmbfile : output file of the created RVE
    """
    def __init__(self, ngrain=100, odf=None, cmbfile='temp.cmb'):
        ## globally shared info ----------------------------- ##
        self.cod = np.loadtxt(odf,skiprows=1)
        self.codt = self.cod.transpose()
        self.resolution = self.codt[0][1] - self.codt[0][0]
        self.inc = self.resolution
        self.p1max = max(self.codt[0])  #phi1
        self.p2max = max(self.codt[1])  #phi2
        self.p3max = max(self.codt[2])  #phi
        self.ngrain = ngrain
        ## -------------------------------------------------- ##

        ## Whether the sample symmetry is applied ----------- ##
        nrot = int(360./self.p1max)        
        print "The max phi1 angle in the COD is %f"%self.p1max
        ngrain = ngrain / nrot
        ## the isotropic random grain number is set to be 1 of nrot
        ## since the number of grain will be multiplied during
        ## sample symmetry application (see def rve_orth above)
        ## -------------------------------------------------- ##
        
        print 'max values %5.1f %5.1f %5.1f'%(
            self.p1max, self.p3max, self.p2max)
        self.cmbfile = cmbfile

        ## isotropic sampling, to which the COD is mapped.
        print 'maximum angles along each axis'
        print 'phi1: %6.2f'%max(self.codt[0])
        print 'phi2: %6.2f'%max(self.codt[1])
        print 'phi : %6.2f'%max(self.codt[2])
        self.gr = random(
            phi1=max(self.codt[0]), phi2=max(self.codt[1]),
            phi=max(self.codt[2]), ngrain=ngrain, iplot=False)

        ## Calculates the RVE
        self.rve = self.cmb()
        pass

    
    def index_odf(self, i, j, k):
        """
        Returns the value of the cod
        corresponding index i,j,k for
        phi1, phi, and phi2, respectively.
        """
        isize= int(self.p1max/self.resolution) + 1
        jsize= int(self.p2max/self.resolution) + 1
        ksize= int(self.p3max/self.resolution) + 1
        
        a = k * isize * jsize
        b = j * isize
        c = i
        ind = int(a) + int(b) + int(c)
        return self.cod[ind][3]        
    
    def search(self,  phi1, phi2, phi):
        """
        Provided the phi1, phi2, and phi,
        find the Euler cube which contains the point.
        Then, estimiates the volume fraction... 
        """
        if phi1>self.p1max or phi2>self.p2max or phi>self.p3max:
            print 'phi1, phi2, phi: %6.3f  %6.3f  %6.3f'%(
                phi1, phi2, phi)
            raise IOError, "The given angle is not available"
    
        i = int(phi1/self.resolution)
        j = int(phi2/self.resolution)
        k = int(phi/self.resolution)
        return self.index_odf(i,j,k)

    def interpolate(self, phi1, phi2, phi):
        """
        Interpolates the volume fraction
        """
        r = [ [], [], [] ]
        r[0].append(phi1 - phi1%self.inc)             #phi
        r[0].append(phi1 - phi1%self.inc  +self.inc)
        r[1].append(phi2 - phi2%self.inc)             #phi2
        r[1].append(phi2 - phi2%self.inc  +self.inc)
        r[2].append(phi  - phi %self.inc)             #phi
        r[2].append(phi  - phi %self.inc  +self.inc)
        value = 0
        x, y, z = 0, 0, 0
        for j in range(2):              #phi1
            for k in range(2):          #phi2
                for l in range(2):      #phi
                    x = self.inc - abs(phi1 - r[0][j])
                    y = self.inc - abs(phi2 - r[1][k])
                    z = self.inc - abs(phi  - r[2][l])
                    value = value + self.search(
                        phi1=r[0][j], phi2=r[1][k],
                        phi=r[2][l]) * x * y * z / (self.inc**3)
                    pass
                pass
            pass
        return value    

    def cmb(self):
        """
        Interpolates the individual grains' volume fraction
        from the discrete crystallographic orientation distribution file
        and returns the representative volume element.
        """
        RVE = [ ]
        tiny = 10**-6
        for i in range(len(self.gr)):
            ## Trick to bypass the case in which the random grain is at the edge
            ## of the given Euler space of the COD.
            if self.gr[i][0]==self.p1max: self.gr[i][0]=self.gr[i][0]-tiny#phi1
            if self.gr[i][1]==self.p3max: self.gr[i][1]=self.gr[i][1]-tiny#phi
            if self.gr[i][2]==self.p2max: self.gr[i][2]=self.gr[i][2]-tiny#phi2
            ## -------------------------------------------------------------- ##

            ## interpolation of the COD value for each of the grain
            value = self.interpolate(phi1=self.gr[i][0],
                                     phi2=self.gr[i][2],
                                     phi=self.gr[i][1])
            RVE.append([self.gr[i][0], self.gr[i][1],
                        self.gr[i][2], value])
            ## ----------------------------------------------------
            pass

        RVE = np.array(RVE)
        ## rve orthogonal symmetry application
        RVE = rve_ortho(cod=self.cod, rve=RVE)

        FILE = open(self.cmbfile, 'w')
        FILE.writelines(
            'dummy\ndummy\ndummy\n%s %i \n'%
            ('B', self.ngrain))
        for i in range(len(RVE)):
            FILE.writelines(
                '%13.7f %13.7f %13.7f   %13.7e \n'%
                (RVE[i][0], RVE[i][1], RVE[i][2], RVE[i][3]))
            pass
        FILE.close()
        return RVE
    pass # end of the class RVE

def main(odf, ngrain, outputfile, iplot=True, irandom=True):
    """
    Arguments

      odf: od filename
      ngrain : number of grains in the RVE
      output : combined file
      iplot =True : flag for plotting
    """
    import upf
    import matplotlib.pyplot as plt
    temp = RVE(ngrain=ngrain, odf=odf, cmbfile=outputfile)

    if iplot==True:
        mypf = upf.polefigure(grains =  temp.rve, csym='cubic')
        mypf.pf(pole=[[1,0,0],[1,1,0],[1,1,1]], mode='contourf', ifig=2)
        plt.show()
        pass

    if irandom==True:
        filename='iso.cmb'
        FILE = open(filename, 'w')
        FILE.writelines('dummy\ndummy\ndummy\n %s %i\n'%('B', int(ngrain)))
        for i in range(len(temp.gr)):
            FILE.writelines('%11.7f %11.7f %11.7f  %10.4e \n'%
                            (temp.gr[i][0], temp.gr[i][1],
                             temp.gr[i][2], temp.gr[i][3]))
            pass
        FILE.close()
        #np.savetxt(filename, temp.gr)
        print 'The random isotropic aggregate is written to %s'%filename
        pass
    pass

if __name__=='__main__':
    import getopt, sys
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                    'i:n:o:sr')
        pass
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        pass

    ##  default options  ##
    iplot = False
    ngrain = 2000
    outputfile ='temp.cmb'
    irandom = False    
    ## ----------------- ##
    
    for o, a in opts:
        if o in ('-i'): inputfile = a
        elif o in ('-n'): ngrain = int(a)
        elif o in ('-o'): outputfile = a
        elif o in ('-s'): iplot=True
        elif o in ('-r'): irandom=True
        else: assert False, 'Unhandled option'
        pass

    main(odf=inputfile, ngrain=ngrain,
         outputfile=outputfile, iplot=iplot,
         irandom=irandom
         )
    
    pass

