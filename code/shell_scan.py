import numpy as np, sys
sys.path.insert(0,'.')
import wfcheck as W
R1,R2=10.0,20.0
h=0.2
x=np.arange(-26,26+h/2,h); y=np.arange(-1.0,26+h/2,h); z=h*np.arange(-2,3)
X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
worldSize=np.sqrt(2*30.0**2+0.36)
def run(v,factor,Rbuff=0.0):
    m=R2/(2*W.GN)*W.C**2*factor
    prof=W.warpshell_profiles(m,R1,R2,Rbuff,0.0,4000,worldSize*1.2)
    g3=W.warpshell_comoving(X,Y,Z,prof,v,doWarp=True)
    out=W.analyse(g3,h,h,h,n_dir=800,eigen=False)
    sl=out['slice']; Tf=out['Tf'][(slice(None),slice(None))+sl]
    mag=np.max(np.abs(Tf),axis=(0,1)); s=out['sampled']
    mask=mag>1e-8*mag.max()
    def mn(k): return np.nanmin(np.where(mask,s[k],np.nan))
    print('v=%.2f compactness=%.2f Rbuff=%.1f: NEC %.2e WEC %.2e SEC %.2e (ratios to max|T| %.3f %.3f %.3f) max|T| %.2e; frac NEC<-1e-3max: %.3f'%(
        v,factor,Rbuff,mn('NEC'),mn('WEC'),mn('SEC'),mn('NEC')/mag.max(),mn('WEC')/mag.max(),mn('SEC')/mag.max(),mag.max(),np.mean((s['NEC']<-1e-3*mag.max())[mask])))
    sys.stdout.flush()
for v,f,rb in ((0.05,1/3,0),(0.1,1/3,0),(0.2,1/3,0),(0.2,0.5,0),(0.2,0.7,0),(0.2,1/3,3.0),(0.5,0.7,0)):
    run(v,f,rb)
