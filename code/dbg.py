import numpy as np, sys
sys.path.insert(0,'.')
import wfcheck as W
exec(open('validate.py').read().split('# ---------- 3.')[0].replace("print(","(lambda *a,**k: None)("))
# pick a point in the wall: x=0, y=R
i=np.argmin(abs(xs-0.0)); j=np.argmin(abs(ys-5.0))
Gn=Gnum[:,:,i,j]; Ge=Gex[:,:,i,j]
beta=vv*(1-np.array(sp.lambdify((x,y,z,v,R,s),f,'numpy')(0.0,5.0,0.0,vv,RR,ss)))
n=np.array([1,-beta,0,0])
print('beta',beta)
print('G_num n n =', n@Gn@n, ' G_ex n n =', n@Ge@n, ' rho_exact/KAPPA =', rho_exact[i,j]/W.KAPPA)
print('pipeline rho_E/KAPPA =', out['rho_E'][i,j]/W.KAPPA)
e,alpha,bu=W.eulerian_frame(out['g'])
print('alpha',alpha[i,j],'beta_up',bu[i,j],'e0',e[0,:,i,j])
print('Gex matrix\n',Ge)
