import numpy as np, sys
sys.path.insert(0,'.')
import wfcheck as W
from scipy.integrate import simpson
h=0.1; v,R,s=0.5,4.0,1.5
xs,ys,X,Y,Z=W.grid(10.0,10.0,h)
def tot(F,sl):
    xi=xs[sl[0]]; yi=ys[sl[1]]; F2=F[sl]; mm=yi>=0
    return simpson(simpson(F2[:,mm]*2*np.pi*yi[mm][None,:],x=yi[mm],axis=1),x=xi)
for A in (2.0,0.5,4.0):
    g3=W.modified_time_comoving(X,Y,Z,v,R,s,A,as_coded=False)
    out=W.analyse(g3,h,h,h,n_dir=50,eigen=False)
    sl=out['slice']; Tf=out['Tf']; N=out['alpha']
    rho=Tf[0,0]; S=Tf[1,1]+Tf[2,2]+Tf[3,3]
    I_null=tot(N*(rho+S/3),sl); I_rho=tot(N*rho,sl); I_rho2=tot(N**2*rho,sl); I_rho0=tot(rho,sl)
    nec=np.nanmin(out['sampled']['NEC'])
    print('A=%.1f: int N(rho+S/3) = %.4e, 2 int N rho = %.4e, ratio %.4f | int N^2 rho = %.3e | int rho = %.3e | NEC min %.2e | lapse range %.3f..%.3f'%(A,I_null,2*I_rho,I_null/(2*I_rho),I_rho2,I_rho0,nec,np.nanmin(N[sl]),np.nanmax(N[sl])))
