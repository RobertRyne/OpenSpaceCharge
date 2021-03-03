       double precision zinitial(6),zfinal(6)
       common/zz/zinitial, zfinal
cryneneriwalstrom       data zinitial/0.,0.,0.,0.,0.,0./
cryne I commented this out because zinitial appears to be unused.
cryne To put it back, the data needs to be initialized in the
cryne block data subroutine in afro.f
