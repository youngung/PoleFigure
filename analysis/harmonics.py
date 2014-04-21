#----------------------------------------------------------------------#
# Fitting grid pole figure with the Harmonics series...                #
def harm_QP(PF=None, dalpha=None, dbeta=None, alpha=None,
           beta=None, l=23, m=3):
    """
    Given the order l, m for the harmonics with theta and phi grid
    returns the coefficients Q_lm

    m shouldn't be larger than l
    l shouldn't be a negative value

    This def returns Q_lm, the harminics coefficients.

    Arguments
    ========
    PF     = Pole figure intensity grid[m, n]
    dalpha = delta alpha
    dbeta  = delta beta
    alpha  = array of alpha   (The tilting angle from the centeral pole)
    beta   = array of thetas  (Rotating angle)
    l      = 23
    m      = 3
    """
    import numpy as np
    P_lm = harm_P(l=l, m=m, beta=beta, alpha=alpha)
    Q_lm = np.zeros((len(beta),len(alpha)), dtype='complex')
    im = complex(0.,1.)

    # for l0=0, l
    # for m0=0, m
    #    for l0 in range(l):        # degree of harmonics
    #        for m0 in range(m):    # order of the harmonic
    for i in range(len(beta)):
        b = beta[i]
        for j in range(len(alpha)):
            a = alpha[j]
            Q_lm[i,j] = Q_lm[i,j] + \
                PF[i,j] *\
                P_lm[i,j] *\
                np.cos(a) *\
                np.exp(-im * m * b) *\
                np.sin(a) * \
                dbeta * \
                dalpha

    return Q_lm, P_lm

def harm_P(l, m, beta, alpha):
    """
    Associated Legendre polynomial

    Arguments
    =========
    l     =
    m     =
    beta  =
    alpha =
    """
    import numpy as np
    P_lm = np.zeros((len(beta), len(alpha)), dtype='complex')
    from scipy.special import sph_harm

    mt, lt, bt, at = [], [], [], []
    for i in range(len(beta)):
        for j in range(len(alpha)):
            mt.append(m); lt.append(l)
            bt.append(beta[i]); at.append(alpha[j])
    dat = sph_harm(mt, lt, bt, at)

    k = 0
    for i in range(len(beta)):
        for j in range(len(alpha)):
            P_lm[i,j] = dat[k]
            k = k + 1
    return P_lm

def harm_PF(PF, l=23):
    """
    PF is assumed to be a 'node' as shown below:

    --------------------------------------------------
    ** nodes = np.zeros((mgrid+1,ngrid+1)) , m: rot, n: tilt
     The mesh nodes is on the central point, 'o', of
    each grid as depicted below:

    90o---o---o
      |   |   |
      .
      .
      .

      o---o---o
      |   |   |
    10o---o---o
      |   |   |   . . .
     5o---o---o
      |   |   |
     0o---o---o
      0   5  10 ....  360

    beta (rotating)
    alpha(tilting)

    Arguments
    =========
    PF : Pole figure nodes from polefigure class
    l  : Order of harmonics
    """
    import numpy as np
    mn, nn  = PF.shape
    dbeta  = 360. / (mn-1)
    dalpha =  90. / (nn-1)
    beta = np.linspace(0.,360., mn) # 0.~360.
    alpha = np.linspace(0.,90., nn) # 0.~90.

    beta = beta     * np.pi / 180.
    alpha = alpha   * np.pi / 180.
    dbeta = dbeta   * np.pi / 180.
    dalpha = dalpha * np.pi / 180.

    im = complex(0, 1)
    PF_ab = np.zeros(np.array(PF).shape, dtype='complex')

    for l0 in range(l):
        print 'l=',l0
        print 'm=',
        for m0 in np.arange(-l0, l0+1):
            print m0,
            Q_lm, P_lm = harm_QP(PF=PF, dalpha=dalpha, dbeta=dbeta,
                                 alpha=alpha, beta=beta, l=l0, m=m0)
            for b0 in range(len(beta)):
                bet = beta[b0]
                for a0 in range(len(alpha)):
                    alph = alpha[a0]
                    PF_ab[b0,a0] = PF_ab[b0,a0] +  \
                        Q_lm[b0,a0] * P_lm[b0,a0] * \
                        np.cos(alph) * np.exp(im * m0 * bet)
        print '\n'
    return PF_ab

def harm_pf(grains=None, filename=None, l=10, dm=7.5, dn=7.5,
            pole=[[1,0,0]], csym='cubic'):
    """
    Something is wrong and I couldn't figure out.
    Further development is defferred.

    Arguments
    =========
    grains   = None
    filename = None
    l        = 10
    dm       = 7.5
    dn       = 7.5
    pole     = [[1, 0, 0]]
    csym     = 'cubic'
    """
    import numpy as np
    import cmb
    from upf import polefigure, circle
    pi = np.pi

    if grains==None and filename==None:
        gr = cmb.random(filename='dum.tex',ngrain=1000,
                        phi1=360., phi=90., phi2=180.)
        mypf = polefigure(grains=grains, csym=csym)
    else: mypf = polefigure(grains=grains, filename=filename, csym=csym)

    mypf.pf(pole=pole, dm=dm, dn=dn)
    hpf = []

    import matplotlib.pyplot as plt
    fact = 2. #size factor
    figsize = (len(pole)*2.*fact, 1.*2.*fact)
    fig = plt.figure(33, figsize=figsize)
    ax = plt.gca()

    for ip in range(len(pole)):
        hpf.append(harm_PF(PF=mypf.pfnodes[ip], l=l))
        print 'max:', np.real(max(hpf[ip].flatten()))
        print 'min:', np.real(min(hpf[ip].flatten()))

    theta = np.linspace(pi, pi/2., (90.)/dn + 1)
    phi = np.linspace(0.,2.*pi,  (360.)/dm+1)
    r = np.sin(theta)/(1-np.cos(theta))
    R, PHI = np.meshgrid(r,phi)
    PHI = PHI + pi/2.
    x = R*np.cos(PHI); y = R*np.sin(PHI)

    #    return x, y, mypf.pfnode[0]
    cnt = ax.contourf(x,y, hpf[0])
    ax.set_frame_on(False)
    ax.set_axis_off()
    ax.set_aspect('equal')
    rx, ry = circle()
    ax.plot(rx, ry, 'k')
    ax.set_xlim(-1.2,1.5)
    ax.set_ylim(-1.2,1.5)

    tcolors = cnt.tcolors
    clev = cnt._levels
    for i in range(len(tcolors)):
        cc = tcolors[i][0][0:3]
                    #if levels==None:
        if ip==len(pole)-1:
                        ## level line
            ax.plot([1.28, 1.35],
                    [1. - i * 0.2, 1. - i * 0.2],
                    color=cc)
                        ## level text
            ax.text(x=1.40, y= 1. - i*0.2 - 0.05,
                    s='%3.2f'%(clev[i]),
                    fontsize=4*fact)
