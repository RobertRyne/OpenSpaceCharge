      subroutine rotv(a,theta,vin,vout)
c
c This subroutine rotates a vector by an angle theta about the axis a.
c
      include 'impli.inc'
c
      dimension aj1(3,3), aj2(3,3), aj3(3,3)
      dimension adotj(3,3), adotj2(3,3), aiden(3,3)
      dimension rotm(3,3)
      dimension a(3), vin(3), vint(3), vout(3)
c
c set up j and identity matrices and store vin
c
      do i=1,3
      do j=1,3
      aj1(i,j)=0.d0
      aj2(i,j)=0.d0
      aj3(i,j)=0.d0
      aiden(i,j)=0.d0
      rotm(i,j)=0.d0
      end do
      end do
      aj1(2,3)=-1.d0
      aj1(3,2)=1.d0
      aj2(1,3)=1.d0
      aj2(3,1)=-1.d0
      aj3(1,2)=-1.d0
      aj3(2,1)=1.d0
      do i=1,3
      aiden(i,i)=1.d0
      vint(i)=vin(i)
      end do
c
c compute adotj
c
      do i=1,3
      do j=1,3
      adotj(i,j)=a(1)*aj1(i,j)+a(2)*aj2(i,j)+a(3)*aj3(i,j)
      end do
      end do
c
c compute adotj2
c
      do i=1,3
      do j=1,3
      adotj2(i,j)=0.d0
      do k=1,3
      adotj2(i,j)=adotj2(i,j)+adotj(i,k)*adotj(k,j)
      end do
      end do
      end do
c
c compute rotation matrix
c
      do i=1,3
      do j=1,3
      rotm(i,j)=aiden(i,j)+(dsin(theta))*adotj(i,j)
      rotm(i,j)=rotm(i,j)+(1.d0-dcos(theta))*adotj2(i,j)
      end do
      end do
c
c comput final result
c
      do i=1,3
      vout(i)=0.d0
      do j=1,3
      vout(i)=vout(i)+rotm(i,j)*vint(j)
      end do
      end do
c
      return
      end
c
***************************************************************************
c
      subroutine erotv(phi,theta,psi,vin,vout)
c
c This subroutine rotates a vector by a rotation described by euler angles.
c
      include 'impli.inc'
c
c
      dimension vin(3), vout(3), vtemp1(3), vtemp2(3)
      dimension ex(3), ey(3), ez(3)
c
c set up unit vectors
c
      do i=1,3
      ex(i)=0.d0
      ey(i)=0.d0
      ez(i)=0.d0
      end do
      ex(1)=1.d0
      ey(2)=1.d0
      ez(3)=1.d0
c
c carry out rotations
c
      call rotv(ez,psi,vin,vtemp1)       
      call rotv(ey,theta,vtemp1,vtemp2) 
      call rotv(ez,phi,vtemp2,vout) 
c
      return
      end
c
******************************************************
c
      subroutine wrtdo
c
c This subroutine writes out points on the design orbit
c
      return
      end
c
