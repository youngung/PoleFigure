import numpy as np
import matplotlib.pyplot as plt
import os
os.sys.path.append('/Users/yj/repo/evpsc-dev/analysis/lib/')
import mpl_lib
import COD_conv
def main(fn='../dat/exp/MV8/FCC.COD',nbin=10):

    odfn='labo_%s.txt'%(fn.split(os.sep)[-1].split('.COD')[0])
    COD_conv.main(fn=fn,odfn=odfn)

    cod, phi1mx, phimx,phi2mx = COD_conv.__cod2labo__(
        fn=fn)

    od_plot(cod,nbin,fn)

    cod = reader(odfn)
    od_plot(cod,nbin,fn)

def od_plot(cod,nbin,fn):
    mx=cod.max();mn=cod.min()

    lev = np.linspace(mn,mx,nbin)

    z, x, y = cod.shape
    xv = np.arange(0,x)*5; yv = np.arange(0,y)*5
    #xv = xv[::-1]; yv=yv[::-1]
    X, Y = np.meshgrid(xv,yv)
    fig = mpl_lib.wide_fig(nw=5,nh=4,h0=0.20,h1=0,w0=0,
                           w1=0,uw=1.2,uh=1.4,down=0.2,
                           right=0.2,
                           iarange=True)
    ax = fig.axes
    for i in range(z):
        l = ax[i].contourf(X,Y,cod[i],
                           lev,cmap=plt.cm.cmap_d['gray_r'])


    ax[-1].set_axis_off()
    ax[-2]
    cbar=plt.colorbar(l,orientation='vertical',shrink=0.8,
                      ticks=map(int,lev[::2]))
    cbar.ax.set_ylabel('Intensity', fontsize=9)
    cl=plt.getp(cbar.ax,'ymajorticklabels')
    plt.setp(cl,fontsize=8)

    plt.draw()
    mpl_lib.rm_all_lab(ax)

    for i in range(z):
        ax[i].set_xlabel(r'$\Phi=%3.1f$'%(xv[i]),
                         dict(fontsize=9))


    fig.text(0.02,0.10,'filename: %s'%fn.split(os.sep)[-1],
             dict(fontsize=8))
    return l

def reader(fn='labo_BCC.txt'):
    cod=np.zeros((19,19,19))
    lines=open(fn,'r').readlines()[1:]
    for i in range(len(lines)):
        phi1,phi2,phi,w=map(float,lines[i].split())
        j = phi1/5.
        i = phi2/5.
        k = phi/5.
        cod[i,j,k] = w # phi2, phi1, phi

    return cod



def ex1():
    fn=['../dat/exp/MV8/FCC.COD','../dat/exp/MV8/BCC.COD']
    l = main(fn[0])
    l = main(fn[1])
    return l
