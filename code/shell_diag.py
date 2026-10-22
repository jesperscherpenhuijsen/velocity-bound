import numpy as np, sys
sys.path.insert(0,'.')
import wfcheck as W
R1,R2=10.0,20.0; m=R2/(2*W.GN)*W.C**2/3
def run(v,h,xmax=26.0):
    x=np.arange(-xmax,xmax+h/2,h); y=np.arange(-1.0,xmax+h/2,h); z=h*np.arange(-2,3)
    X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
    worldSize=np.sqrt(2*30.0**2+0.36)
    prof=W.warpshell_profiles(m,R1,R2,0.0,0.0,4000,worldSize*1.2)
    g3=W.warpshell_comoving(X,Y,Z,prof,v,doWarp=True)
    out=W.analyse(g3,h,h,h,n_dir=2000)
    sl=out['slice']; Tf=out['Tf'][(slice(None),slice(None))+sl]
    mag=np.max(np.abs(Tf),axis=(0,1)); nec=out['sampled']['NEC']
    i,j=np.unravel_index(np.nanargmin(nec),nec.shape)
    xs=x[sl[0]]; ys=y[sl[1]]
    r=np.hypot(xs[i],ys[j])
    print('v=%.2f h=%.2f: NEC min sampled %.3e at (x,y)=(%.1f,%.1f) r=%.2f, max|T| %.3e, ratio %.3f'%(v,h,nec[i,j],xs[i],ys[j],r,mag.max(),nec[i,j]/mag.max()))
    T=Tf[:,:,i,j]
    np.set_printoptions(precision=3,linewidth=150)
    print(' frame tensor / 1e40:\n',T/1e40)
    eta=np.diag([-1.,1,1,1]); w,V=np.linalg.eig(eta@T)
    print(' eigenvalues of T^a_b /1e40:',w/1e40)
    # also: fraction of masked points with sampled NEC < -1e-3 max|T|
    mask=mag>1e-8*mag.max()
    print(' fraction of masked points with NEC < -1e-3 max|T|: %.4f'%np.mean((nec< -1e-3*mag.max())[mask]))
    # radial location histogram of violating points
    XX,YY=np.meshgrid(xs,ys,indexing='ij'); RR=np.hypot(XX,YY)
    viol=(nec< -1e-3*mag.max())&mask
    if viol.any(): print(' violating points r range: %.2f .. %.2f m'%(RR[viol].min(),RR[viol].max()))
for v,h in ((0.2,0.2),(0.2,0.1),(0.5,0.2),(0.02,0.2)):
    run(v,h)
