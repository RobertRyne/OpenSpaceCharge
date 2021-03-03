*******************************************************************************
* DUMMY                                                                       *
* These are dummy routines to please AFRO and OPTI                            *
*******************************************************************************
c
*****************************************************************************
c Routines fo please AFRO
*****************************************************************************
c
      subroutine bell
      write(6,*)' BEEP'
      return
      end
c
****************************************************************************
c
      subroutine wmrt(p)
      dimension p(*)
      write(6,*)' in wmrt'
      return
      end
c
***********************************************************************
c
      subroutine subctr(p)
      dimension p(*)
      write(6,*)' in subctr'
      return
      end
c
*****************************************************************************
c
      subroutine fwa(p)
      write(6,*) 'in dummy routine fwa' 
      return
      end
c
*****************************************************************************
c
      subroutine dism(p,h,mh)
      write(6,*) 'in dummy routine dism' 
      return
      end
c
*****************************************************************************
c
      subroutine rmap(p,h,mh)
      write(6,*) 'in dummy routine rmap' 
      return
      end
c
*****************************************************************************
c
      subroutine flag(p)
      write(6,*) 'in dummy routine flag'
      return
      end      
c
************************************************************************
c Routines to please OPTI
**************************************************************************
c
cryne 8/17/02      subroutine dcopy(nx,fu,maxf,y,i)
cryne 8/17/02      write(6,*) 'in dummy routine dcopy'
cryne 8/17/02      return
cryne 8/17/02      end     
c
***************************************************************************
c
      subroutine dposl( fjac,maxf,nx,step)
      write(6,*) 'in dummy routine dposl'
      return
      end     
c
***************************************************************************
c
      subroutine dqrdc(fjac,maxf,nf,nx,qraux,jdum,dum,i)
      write(6,*) 'in dummy routine dqrdc'
      return
      end     
c
***************************************************************************
c
      subroutine dqrsl(fj,ma,nf,nx,qr,fc,du,r,sp,rs,d,i,in)
      write(6,*) 'in dummy routine dqrsl'
      return
      end     
