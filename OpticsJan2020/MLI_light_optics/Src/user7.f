      subroutine user7(p,g,hmh)
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'ind.inc'
      include 'prodex.inc' 
      dimension p(*), g(*), hmh(*)
      write(6,*)' In USER7 test'
      lun = p(1)
      mbgn = p(2)
      if(mbgn.le.0) mbgn = 1
      num = p(3)
      if(num.le.0) num = monoms
      write(6,*)'PRODEX:'
      do i=0,923
      write(96,1234)i,prodex(1:6,i)
      enddo
 1234 format(i5,4x,6i10)
c      call sixi(lun,prodex,monoms)
      write(6,*)imaxi,' = IMAXI'
      write(6,*)' JV, INDEX1, INDEX2'
      call onei(lun,jv(mbgn),index1(mbgn),index2(mbgn),num)
      return
      end
c--------------------------------
      subroutine sixi(lun,ival,num)
      dimension ival(6,*), ibuf(6)
      do 100 k = 1,num+1
         kv = 0
         do 50 j = 1,6
           ibuf(j) = ival(j,k)
           if(ibuf(j).ne.0) kv = kv + 1
  50     continue
         if(kv.gt.0) write(lun,77) k,ibuf
  77     format(7i5)
 100  continue
      return
      end
c--------------------------------------
      subroutine onei(lun,ival,jval,kval,num)
      dimension ival(*), jval(*), kval(*)
      do 100 k = 1,num
      if(ival(k).ne.0) write(lun,77) k,ival(k),jval(k),kval(k)
  77  format(4i8)
 100  continue
      return
      end
c
c     subroutine user7(p,g,hmh)
c     dimension p(*), g(*), hmh(*)
c     return
c     end
