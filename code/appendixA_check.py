import numpy as np, sys
sys.path.insert(0,'.')
import wfcheck as W
from scipy.integrate import simpson
h=0.1; v,R,s=0.5,4.0,1.5
xs,ys,X,Y,Z=W.grid(10.0,10.0,h)
g3=W.alcubierre_comoving(X,Y,Z,v,R,s)
out=W.analyse(g3,h,h,h,n_dir=50,eigen=False)
sl=out['slice']; Tf=out['Tf']
rho=Tf[0,0]; S=Tf[1,1]+Tf[2,2]+Tf[3,3]
# K_ij = (d_i b_j + d_j b_i)/2 with b = (v(1-f),0,0), 4th-order derivatives from the same slab
r=np.sqrt(X**2+Y**2+Z**2); f=W.shape_alcubierre(r,R,s); bx=v*(1-f)
def d(F,ax):
    out_=np.full_like(F,np.nan)
    if ax==0: out_[2:-2]=(-(F[4:]-F[:-4])+8*(F[3:-1]-F[1:-3]))/(12*h)
    if ax==1: out_[:,2:-2]=(-(F[:,4:]-F[:,:-4])+8*(F[:,3:-1]-F[:,1:-3]))/(12*h)
    if ax==2: out_[:,:,2:-2]=(-(F[:,:,4:]-F[:,:,:-4])+8*(F[:,:,3:-1]-F[:,:,1:-3]))/(12*h)
    return out_
bxx=d(bx,0); bxy=d(bx,1); bxz=d(bx,2)
K=bxx                      # trace
KK=bxx**2+0.5*(bxy**2+bxz**2)
dK=d(bxx,0)                # d_x K
bdK=bx*dK
k=2
rho_f=(K**2-KK)[:,:,k]/(16*np.pi)*W.C**4/W.GN
S_f=(-rho_f*W.GN/W.C**4 - (KK+bdK)[:,:,k]/(4*np.pi))*W.C**4/W.GN
m=(slice(4,-4),slice(4,-4))
print('rho_E: pipeline vs Hamiltonian constraint, rel err %.1e'%(np.nanmax(abs(rho[m]-rho_f[m]))/np.nanmax(abs(rho_f[m]))))
print('S_E  : pipeline vs trace-evolution identity, rel err %.1e'%(np.nanmax(abs(S[m]-S_f[m]))/np.nanmax(abs(S_f[m]))))
def tot(F):
    xi=xs[m[0]]; yi=ys[m[1]]; F2=F[m]; mm=yi>=0
    return simpson(simpson(F2[:,mm]*2*np.pi*yi[mm][None,:],x=yi[mm],axis=1),x=xi)
I1=tot(rho); I2=tot(rho+S/3)
print('int rho_E dV = %.4e J ; int (rho_E + S/3) dV = %.4e J ; ratio = %.4f (identity predicts 2)'%(I1,I2,I2/I1))
print('max over plane of (rho_E+S/3): %.3e ; min: %.3e'%(np.nanmax((rho+S/3)[m]),np.nanmin((rho+S/3)[m])))
