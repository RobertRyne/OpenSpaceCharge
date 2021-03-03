! 12/15/2004: Original version commented out by RDR
!     common/stmap/sf1(monoms),sf2(monoms),sf3(monoms),
!    #             sf4(monoms),sf5(monoms),
!    #             sm1(6,6),sm2(6,6),sm3(6,6),
!    #             sm4(6,6),sm5(6,6)
!     dimension sfa(monoms,5), sma(6,6,5)
!     equivalence (sfa,sf1), (sma,sm1)
      common/stmap/storedpoly(monoms,20),storedmat(6,6,20)
