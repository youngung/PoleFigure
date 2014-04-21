"""
Convert Sectioned OD to LABOTEXT format

- example
>>> import main(fn='A_steel.cod ', odfn='labo.txt')
"""

print __doc__

def reader(fn=None):
    """
    Reads a file consisting of sections of 'phi2=constant'
    and return total section blocks and maximum phi2

    Argument
    ========
    fn = None
    """
    datlines = open(fn, 'r').readlines()
    blocks   = []
    temp     = []
    for i in range(len(datlines)):
        temp.append(datlines[i])
        if len(datlines[i]) < 3:
            blocks.append(temp)
            temp = []

    maxphi2 = float(blocks[-2][1].split('phi2=')[1])
    dum = float(blocks[-3][1].split('phi2=')[1])
    dphi2 = maxphi2 - dum
    print 'maximum phi2', maxphi2
    print 'phi2 increment:', dphi2
    return blocks, maxphi2

def readblock(block):
    """
    Read a block and returns
    1) section's intensity mesh
    2) Phi  max
    3) phi1 max

    Argument
    ========
    block
    """
    import numpy as np
    header = block[0].split()
    dphi    = float(block[1][5:10])
    phimx   = float(block[1][10:15])
    dphi1   = float(block[1][15:20])
    phi1mx  = float(block[1][20:25])
    print 'dphi, phimx, dphi1, phi1mx',
    print dphi, phimx, dphi1, phi1mx

    if phi1mx==180: pass
    elif phi1mx==90: pass
    else: raise IOError, 'Not expected maximum phi value...'

    phi  = np.zeros(np.arange(0., phimx  + 0.001, dphi ).shape)
    phi1 = np.zeros(np.arange(0., phi1mx + 0.001, dphi1).shape)

    section = np.zeros((len(phi), len(phi1)))

    block = block[2:]

    if phi1mx==180:
        for i in range((len(block)-1)/2):
            arim = block[i*2][:18*4+1][1:] + block[i*2+1][:19*4+1][1:]
            for j in range(len(arim[::4])):
                section[i,j] = float(arim[4*j:4*j+4])
    elif phi1mx==90:
        for i in range(len(block)-1):
            arim = block[i][:19*4+1][1:]
            for j in range(len(arim[::4])):
                section[i,j] = float(arim[4*j:4*j+4])

    # # block = block[::-1][0:]
    # # block = block[::-1]

    # for i in range(len(block)-1):
    #     dum = block[i].split()
    #     section[i] = map(float, dum)

    section = section.T # (phi, phi1) -> (phi1, phi)
    return section, phimx, phi1mx

def __cod2labo__(fn=None):
    """
    Argument
    =========
    fn = None
    """
    import numpy as np
    blocks, phi2mx = reader(fn=fn)

    sct, phimx, phi1mx = readblock(blocks[0])

    nphi2 = len(blocks)-1
    nphi1 = len(sct)
    nphi  = len(sct[0])

    cod = np.zeros((nphi1, nphi, nphi2))

    # len(blocks)
    for ib in range(len(blocks)-1): # nphi2
        section, pmx, pmx1 = readblock(blocks[ib])
        for i in range(len(section)): # nphi1
            for j in range(len(section[i])): #nphi
                cod[i,j,ib] = section[i,j]

    return cod, phi1mx, phimx, phi2mx

def main(fn=None, odfn ='labo.txt', ssym=None):
    """
    Writes COD to Labotex convention file

    arguments
    =========
    fn   = None
    odfn = 'labo.txt'
    ssym = None, 0, 1, 2
    """
    import numpy as np
    cod, phi1mx, phimx, phi2mx = __cod2labo__(fn=fn)
    print phi1mx, phimx, phi2mx

    f = open(odfn, 'w')
    phi1 = np.linspace(0, phi1mx+0.001, len(cod)) # phi2
    phi  = np.linspace(0, phimx+0.001, len(cod[0])) #
    phi2  = np.linspace(0, phi2mx+0.001, len(cod[0][0]))

    cod = np.array(cod).swapaxes(2,0) # phi2, phi, phi1
    cod = np.array(cod).swapaxes(0,1) # phi, phi2, phi1

    print 'cod.shape', cod.shape

    f.write('phi1   phi2  phi  COD\n')
    for i in range(len(cod)): # phi
        p = phi[i]
        for j in range(len(cod[i])): #phi2
            p2 = phi2[j]
            for k in range(len(cod[i][j])): #phi1
                p1 = phi1[k]
                f.write('%4.1f  %4.1f  %4.1f %6i \n'%(
                        p1, p2, p, cod[i][j][k]))

    f.close()
    import os
    os.sys.path.append('/Users/yj/Repo/EVPSC_TRIP/pyscripts/pf')
    # import cmb
    # cmb.RVE(odf=odfn, ngrain=2000,
    #         outputfile='dum.tex',ssym=3)

