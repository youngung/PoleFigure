## kurjumov sachs

""" {111}<110> // {110}<-111> """

import numpy as np
def rot_nb(n1,b1,n2,b2):
    r_ato1 = np.zeros((3,3))
    r_mto2 = np.zeros((3,3))
    r      = np.zeros((3,3))

    t1 = np.cross(n1,b1)
    t2 = np.cross(n2,b2)

    for l in range(3):
        r_ato1[0,l]=b1[l]
        r_ato1[1,l]=t1[l]
        r_ato1[2,l]=n1[l]

        r_mto2[0,l]=b2[l]
        r_mto2[1,l]=t2[l]
        r_mto2[2,l]=n2[l]

    #r = np.dot(r_mto2.T,r_ato1) # martensite <- austenite

    # for i in range(3):
    #     for j in range(3):
    #         for k in range(3):
    #             r[i,j] = r[i,j] + r_mto2[k,i] * \
    #                      r_ato1[k,j]


    r_2tom = r_mto2.T
    for i in range(3):
        for j in range(3):
            for k in range(3):
                r[i,j] = r[i,j] + r_2tom[i,k] * r_ato1[k,j]

    return r

def ksr():
    """
    return a ks \delta g
    """
    n1 = np.array([ 1., 1., 1.])/np.sqrt(3.)
    b1 = np.array([ 1.,-1., 0.])/np.sqrt(2.)
    n2 = np.array([ 0.,1., 1.])/np.sqrt(2.)
    b2 = np.array([ 1., -1., 1.])/np.sqrt(3.)

    n1 = n1 / np.sqrt(sum(n1**2.))
    b1 = b1 / np.sqrt(sum(b1**2.))
    n2 = n2 / np.sqrt(sum(n2**2.))
    b2 = b2 / np.sqrt(sum(b2**2.))

    t1=np.cross(n1,b1)
    t2=np.cross(n2,b2)

    print n1,b1,t1
    print n2,b2,t2

    print np.dot(n1,b1), np.dot(n1,t1), np.dot(b1,t1)
    print np.dot(n2,b2), np.dot(n2,t2), np.dot(b2,t2)

    r = rot_nb(n1,b1,n2,b2)


    print '\n\n\n'
    print np.dot(r,n1)
    print np.dot(r,b1)
    print '\n'
    print np.dot(r.T,n2)
    print np.dot(r.T,b2)


    return r # dr_m<-a

def main():
    """
    Return all 24 variant orientations
    """
    r = ksr() # dr_m<-a
    from symmetry import sym
    hs = sym.cubic()

    mats=[]
    for i in range(len(hs)):
        mat=np.dot(r,hs)
        mats.append(mat)
    
    return np.array(mats)

def aust2ksalpha(a_xtal):
    """
    xtal in the form of [[phi1,phi,phi2,int],[...],...]
    """
    from euler import euler
    from symmetry import sym
    hs = sym.cubic()

    eul = euler.euler
    ks1 = ksr() # m<-a
    #ksr = main()
    m_xtal = []
    for ig in range(len(a_xtal)):
        phi1,phi,phi2,it = a_xtal[ig] # ca<-sa
        a = eul(ph=phi1,th=phi,tm=phi2,echo=False) # ca<-sa
        m1 = np.dot(ks1,a)  # m1 <- ca_aust <- sa

        for i in range(len(hs)):
            mi = np.dot(hs[i],m1)
            phi1,phi,phi2=eul(a=mi,echo=False)
            m_xtal.append([phi1,phi,phi2,it])

        # for iv in range(len(ksr)):
        #     b = np.dot(ksr[iv],a)
        #     b = eul(a=b,echo=False)
        #     m_xtal.append([b[0],b[1],b[2],it/len(ksr)])

    return np.array(m_xtal)
        
def pf_ks(aust_hkl=[],aust_uvw=[],w0=10.,ngr=100):
    """
    Given for the ideal component of austenite,
    plot pole figures of martensite in the KS relationship
    """

    from euler import euler
    eul = euler.euler
    
    from analysis import ideal_comp
    from pf import upf

    if w0!=0:
        a_xtal = ideal_comp.main(
            aust_hkl,aust_uvw,w0,ngr=ngr,ifig=1)
    elif w0==0:
        # mat = ideal_comp.miller2mat(hkl=aust_hkl,uvw=aust_uvw)
        # phi1,phi,phi2 = eul(a=mat,echo=False)
        phi1,phi,phi2 = ideal_comp.miller2euler(aust_hkl,aust_uvw)
        a_xtal = [[phi1,phi,phi2,1]]

    #a_xtal [[phi1,phi,phi2,intensity],[..],...]
    
    m_xtal = aust2ksalpha(a_xtal)
    #return m_xtal
    mypf=upf.polefigure(grains=m_xtal,csym='cubic')
    mypf.pf(cmode='gray_r',ifig=2,mode='dot')
