c mstack = maximum depth of stack
      parameter (mstack=20)
c np    = stack pointer
c ntype = type of top stack element
c ith   = index of "   "       "    in its array
c mtype = type of actual slot of ith
c jth   = index of  "     "    "  "  in its array
c nslot(k) = actual slot of kth stack element
c loop(k)  = repetition factor of kth stack element
c lstac(k) = name of kth stack element
      common /stack/ np,ntype,ith,mtype,jth,nslot(mstack),loop(mstack)
      common /stac/  lstac(mstack)
      character*16 lstac
