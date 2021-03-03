      common/fitdat/dz(4),tux,tuy,tus,cx,cy,qx,qy,hh,vv,tt,hv,ht,vt,
     *ax,bx,gx,ay,by,gy,at,bt,gt,ex,ey,et,wex,wey,wet,fx,xb,xa,xu,xd,
c     *fy,yb,ya,yu,yd,ft,tb,ta,tu,td 
c AJD 6/5/93
     *fy,yb,ya,yu,yd,fte,teb,tea,teu,ted
      dimension fitval(47)
      equivalence (fitval,dz)
c Note: The quantities ft,ta,tb,tu,td have been changed to fte,tea,teb,teu,ted
c in order to avoid name conflicts with with various "ta" arrays in various 
c analysis routines.  However, these quantities are still referenced by the keys
c 'ft','tb','ta','tu','td' defined in the include file keyset.inc - AJD
c 6/5/93.
