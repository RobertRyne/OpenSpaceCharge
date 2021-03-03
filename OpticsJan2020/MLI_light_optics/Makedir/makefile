VPATH=../Src

.SUFFIXES:
.SUFFIXES: .o .f .f90


# implicit compile rules
.f.o:
	$(Fortran) $<


.f90.o:
	$(FortranFree) $<

Fortran     = gfortran -c -I../Includes
FortranFree = gfortran -c -I../Includes
FortranMain = gfortran -o ../mli.x -I../Includes

OBJS = afro_mod.o parallel_mod.o constants_mod.o liea_mod.o spch3d_mod.o \
       e_gengrad_mod.o gengrad_mod.o timer_mod.o multitrack_mod.o \
       curve_fit.o greenfn_mod.o \
       afro.o anal.o base.o bessjm.o book.o boundp3d.o \
       cfbdang.o cfqd.o coil.o comm.o cons.o \
       diagnostics.o dist.o dummy.o dumpin.o \
       ebcomp.o elem.o env.o euclid.o fftpkgq.o fparser.o \
       genm.o gensol.o hamdrift.o inpu.o integ.o iron.o \
       liea.o linpak_old.o myblas.o xerbla.o \
       magnet.o math.o meri.o mygenrec5.o opti.o \
       parameters.o proc.o pure.o rfgap.o \
       setbound.o sif.o spch2d.o spch3d.o sss.o trac.o \
       user.o user7.o usubs.o wakefld.o xtra_notgnu.o \
       spch3d_dummy.o spch3d_chombo_dummy.o mpi.o optics.o 

# explicit rule to build executable
../mli.x: $(OBJS)
	$(FortranMain) $(OBJS)

# other targets
.PHONY: clean
clean:
	@\rm -f *.o *.mod *.d


# module dependencies
afro.o anal.o book.o boundp3d.o cfbdang.o cfqd.o coil.o comm.o diagnostics.o dist.o dumpin.o elem.o env.o gensol.o hamdrift.o inpu.o integ.o liea.o math.o mygenrec5.o proc.o rfgap.o setbound.o sif.o spch2d.o spch3d.o spch3d_essl.o trac.o user.o usubs.o wakefld.o : afro_mod.o

ebcomp.o e_gengrad_mod.o : constants_mod.o

afro.o ebcomp.o : e_gengrad_mod.o

sif.o : fparser.o

e_gengrad_mod.o : gengrad_mod.o

afro.o anal.o base.o book.o cfbdang.o comm.o cons.o diagnostics.o dist.o dumpin.o elem.o env.o genm.o gensol.o hamdrift.o inpu.o liea.o math.o mygenrec5.o proc.o pure.o rfgap.o spch3d.o sss.o trac.o user.o user7.o usubs.o : liea_mod.o

afro.o afro_mod.o anal.o book.o diagnostics.o dist.o dumpin.o env.o fftpkgq.o inpu.o pure.o rfgap.o sif.o spch3d.o spch3d_essl.o timer_mod.o : parallel_mod.o

fparser.o sif.o : parameters.o

afro.o spch3d.o : spch3d_mod.o

afro.o diagnostics.o fftpkgq.o spch3d.o spch3d_essl.o trac.o : timer_mod.o
 

# include file dependencies -- use the following target to regenerate the dependencies
make.depends:
	@grep -i '^[ \t]*include[ \t]' ../Src/*.f | grep -v 'mpif.h' | sort | uniq | awk '{sub("../Src/","");sub(".f:",".o:")gsub("'\''","");$$2=""; print $$1,"../Includes/" $$3}'
	@grep -i '^[ \t]*include[ \t]' ../Src/*.f90 | grep -v 'mpif.h' | sort | uniq | awk '{sub("../Src/","");sub(".f90:",".o:")gsub("'\''","");$$2=""; print $$1,"../Includes/" $$3}'

afro.o: ../Includes/codes.inc
afro.o: ../Includes/core.inc
afro.o: ../Includes/deriv.inc
afro.o: ../Includes/expon.inc
afro.o: ../Includes/files.inc
afro.o: ../Includes/fitbuf.inc
afro.o: ../Includes/frnt.inc
afro.o: ../Includes/impli.inc
afro.o: ../Includes/ind.inc
afro.o: ../Includes/infin.inc
afro.o: ../Includes/labpnt.inc
afro.o: ../Includes/lims.inc
afro.o: ../Includes/loop.inc
afro.o: ../Includes/map.inc
afro.o: ../Includes/maxcat.inc
afro.o: ../Includes/nturn.inc
afro.o: ../Includes/parset.inc
afro.o: ../Includes/pie.inc
afro.o: ../Includes/previous.inc
afro.o: ../Includes/setref.inc
afro.o: ../Includes/stack.inc
afro.o: ../Includes/talk.inc
afro.o: ../Includes/time.inc
afro.o: ../Includes/zeroes.inc
anal.o: ../Includes/buffer.inc
anal.o: ../Includes/codes.inc
anal.o: ../Includes/core.inc
anal.o: ../Includes/extalk.inc
anal.o: ../Includes/files.inc
anal.o: ../Includes/fitdat.inc
anal.o: ../Includes/hmflag.inc
anal.o: ../Includes/impli.inc
anal.o: ../Includes/loop.inc
anal.o: ../Includes/param.inc
anal.o: ../Includes/parset.inc
anal.o: ../Includes/pie.inc
anal.o: ../Includes/supres.inc
anal.o: ../Includes/usrdat.inc
book.o: ../Includes/dr.inc
book.o: ../Includes/drl.inc
book.o: ../Includes/expon.inc
book.o: ../Includes/files.inc
book.o: ../Includes/frnt.inc
book.o: ../Includes/id.inc
book.o: ../Includes/impli.inc
book.o: ../Includes/ind.inc
book.o: ../Includes/ind3.inc
book.o: ../Includes/iprod.inc
book.o: ../Includes/ja3.inc
book.o: ../Includes/len.inc
book.o: ../Includes/len3.inc
book.o: ../Includes/lims.inc
book.o: ../Includes/maxcat.inc
book.o: ../Includes/order.inc
book.o: ../Includes/pie.inc
book.o: ../Includes/pq.inc
book.o: ../Includes/prodex.inc
book.o: ../Includes/sr.inc
book.o: ../Includes/srl.inc
book.o: ../Includes/symp.inc
book.o: ../Includes/time.inc
book.o: ../Includes/vblist.inc
cfbdang.o: ../Includes/impli.inc
cfbdang.o: ../Includes/parset.inc
comm.o: ../Includes/expon.inc
comm.o: ../Includes/files.inc
comm.o: ../Includes/impli.inc
comm.o: ../Includes/lims.inc
comm.o: ../Includes/maxcat.inc
comm.o: ../Includes/param.inc
comm.o: ../Includes/stmap.inc
cons.o: ../Includes/impli.inc
cons.o: ../Includes/parset.inc
diagnostics.o: ../Includes/files.inc
diagnostics.o: ../Includes/impli.inc
diagnostics.o: ../Includes/map.inc
dist.o: ../Includes/buffer.inc
dist.o: ../Includes/files.inc
dist.o: ../Includes/impli.inc
dist.o: ../Includes/parset.inc
dumpin.o: ../Includes/codes.inc
dumpin.o: ../Includes/files.inc
dumpin.o: ../Includes/impli.inc
dumpin.o: ../Includes/incmif.inc
ebcomp.o: ../Includes/ebdata.inc
ebcomp.o: ../Includes/impli.inc
elem.o: ../Includes/files.inc
elem.o: ../Includes/impli.inc
elem.o: ../Includes/parset.inc
elem.o: ../Includes/pie.inc
elem.o: ../Includes/symp.inc
env.o: ../Includes/impli.inc
euclid.o: ../Includes/impli.inc
genm.o: ../Includes/hmflag.inc
genm.o: ../Includes/impli.inc
genm.o: ../Includes/lims.inc
gensol.o: ../Includes/combs.inc
gensol.o: ../Includes/files.inc
gensol.o: ../Includes/hmflag.inc
gensol.o: ../Includes/impli.inc
gensol.o: ../Includes/parset.inc
gensol.o: ../Includes/sol.inc
hamdrift.o: ../Includes/impli.inc
inpu.o: ../Includes/buffer.inc
inpu.o: ../Includes/codes.inc
inpu.o: ../Includes/core.inc
inpu.o: ../Includes/drl.inc
inpu.o: ../Includes/expon.inc
inpu.o: ../Includes/files.inc
inpu.o: ../Includes/impli.inc
inpu.o: ../Includes/incmif.inc
inpu.o: ../Includes/loop.inc
inpu.o: ../Includes/parset.inc
inpu.o: ../Includes/pbkh.inc
inpu.o: ../Includes/srl.inc
liea.o: ../Includes/expon.inc
liea.o: ../Includes/impli.inc
liea.o: ../Includes/ind.inc
liea.o: ../Includes/iprod.inc
liea.o: ../Includes/len.inc
liea.o: ../Includes/lims.inc
liea.o: ../Includes/pbkh.inc
liea.o: ../Includes/prodex.inc
liea.o: ../Includes/symp.inc
liea.o: ../Includes/vblist.inc
math.o: ../Includes/buffer.inc
math.o: ../Includes/expon.inc
math.o: ../Includes/id.inc
math.o: ../Includes/impli.inc
math.o: ../Includes/len.inc
math.o: ../Includes/vblist.inc
mygenrec5.o: ../Includes/combs.inc
mygenrec5.o: ../Includes/files.inc
mygenrec5.o: ../Includes/hmflag.inc
mygenrec5.o: ../Includes/impli.inc
mygenrec5.o: ../Includes/parset.inc
mygenrec5.o: ../Includes/quadpn.inc
mygenrec5.o: ../Includes/recmul.inc
myprot5.o: ../Includes/impli.inc
myprot5.o: ../Includes/param.inc
myprot5.o: ../Includes/parm.inc
myprot5.o: ../Includes/pie.inc
opti.o: ../Includes/impli.inc
opti.o: ../Includes/minvar.inc
opti.o: ../Includes/nlsvar.inc
proc.o: ../Includes/aimdef.inc
proc.o: ../Includes/amdiip.inc
proc.o: ../Includes/buffer.inc
proc.o: ../Includes/codes.inc
proc.o: ../Includes/files.inc
proc.o: ../Includes/fitbuf.inc
proc.o: ../Includes/fitdat.inc
proc.o: ../Includes/impli.inc
proc.o: ../Includes/keyset.inc
proc.o: ../Includes/labpnt.inc
proc.o: ../Includes/map.inc
proc.o: ../Includes/merit.inc
proc.o: ../Includes/parset.inc
proc.o: ../Includes/psflag.inc
proc.o: ../Includes/status.inc
proc.o: ../Includes/stmap.inc
proc.o: ../Includes/usrdat.inc
proc.o: ../Includes/xvary.inc
pure.o: ../Includes/dr.inc
pure.o: ../Includes/impli.inc
pure.o: ../Includes/sr.inc
rfgap.o: ../Includes/files.inc
rfgap.o: ../Includes/fitbuf.inc
rfgap.o: ../Includes/impli.inc
rfgap.o: ../Includes/map.inc
rfgap.o: ../Includes/usrdat.inc
setbound.o: ../Includes/impli.inc
sif.o: ../Includes/codes.inc
sif.o: ../Includes/files.inc
sif.o: ../Includes/impli.inc
sif.o: ../Includes/pie.inc
spch3d.o: ../Includes/files.inc
spch3d.o: ../Includes/impli.inc
spch3d.o: ../Includes/map.inc
sss.o: ../Includes/files.inc
trac.o: ../Includes/deriv.inc
trac.o: ../Includes/expon.inc
trac.o: ../Includes/files.inc
trac.o: ../Includes/impli.inc
trac.o: ../Includes/ind.inc
trac.o: ../Includes/ind3.inc
trac.o: ../Includes/ja3.inc
trac.o: ../Includes/len.inc
trac.o: ../Includes/len3.inc
trac.o: ../Includes/lims.inc
trac.o: ../Includes/pbkh.inc
trac.o: ../Includes/talk.inc
trac.o: ../Includes/vblist.inc
user.o: ../Includes/expon.inc
user.o: ../Includes/files.inc
user.o: ../Includes/fitdat.inc
user.o: ../Includes/impli.inc
user.o: ../Includes/linbuf.inc
user.o: ../Includes/map.inc
user.o: ../Includes/parset.inc
user.o: ../Includes/sigbuf.inc
user.o: ../Includes/sincos.inc
user.o: ../Includes/taylor.inc
user.o: ../Includes/usrdat.inc
user7.o: ../Includes/impli.inc
user7.o: ../Includes/ind.inc
user7.o: ../Includes/prodex.inc
usubs.o: ../Includes/actpar.inc
usubs.o: ../Includes/codes.inc
usubs.o: ../Includes/core.inc
usubs.o: ../Includes/expon.inc
usubs.o: ../Includes/files.inc
usubs.o: ../Includes/impli.inc
usubs.o: ../Includes/linbuf.inc
usubs.o: ../Includes/loop.inc
usubs.o: ../Includes/map.inc
usubs.o: ../Includes/parset.inc
usubs.o: ../Includes/pbkh.inc
usubs.o: ../Includes/sigbuf.inc
usubs.o: ../Includes/taylor.inc
usubs.o: ../Includes/usrdat.inc
usubs.o: ../Includes/vblist.inc
xtra_notgnu.o: ../Includes/files.inc
xtra_notgnu.o: ../Includes/impli.inc
xtra_notgnu.o: ../Includes/time.inc
