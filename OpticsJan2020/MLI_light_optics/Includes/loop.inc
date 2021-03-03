c mim(i) points to entities (lumps, elements, user routines, etc.) in a loop
c joy    is the number of things in mim
c nloop  is the index of the actual loop in /items/
      parameter (joymax=1000)
      common/loop/mim(joymax),joy,nloop
