      PROGRAM KINISOT 
C     Calculation of harmonic isotope effects based on the normal vibrational wavenumbers
C     of Reactant and an isotopomer, and  transition state and the equivalent isotopomer.
C
C     Original version: H S Rzepa 1975, Austin Texas.
C     Literature references: doi: 10.1021/ja00493a008  and doi: 10.1021/ja00486a013 (1978).
C     Formatting and I/O changes: H S Rzepa, 1980, Imperial College.
C     Comments:  H S Rzepa, July, 2015, Imperial College.
C
C     Compilation:   gfortran kinisot.f  -o kinisot.exe (Tested 3 July, 2015 on OS X 10.10.4). 
C     Compiler from http://hpc.sourceforge.net/   Version: gcc version 5.1.0 (GCC) 
C
C     Inputs:  isotope.dat.  Output: isotope.out
C
C     Format of Input file
C     1:  48 288 298 300 310 320 330 340 350 360 370     where N=48 is the number of atoms
C                                                        followed by  10 temperatures. If 
C                                                        number of atoms is zero, program stops.
C     2: Title for normal isotope reactant
C     3: 3N-6 values for the normal mode wavenumbers for reactant.
C     4: Title for isotopomer of reactant
C     5: 3N-6 values for the normal mode wavenumbers for isotopomeric reactant
C     6: Title for normal isotope transition state
C     7: 3N-6 values for the normal mode wavenumbers for transition state, with imaginary mode
C        listed first (as a -ve number)
C     8: 3N-6 values for the normal mode wavenumbers for isotopomeric transition state, with imaginary mode
C        listed first (as a -ve number)
C     9: Repeat card 1 for new isotopomers.
C
C     Program, input and output archived as  doi: 10.5281/zenodo.19272 3 July, 2015.
C
      DIMENSION WR1(500), WR2(500), WT1(500),WT2(500),TEMPS(10)
      character name(4)*70,filename*31
      LOGICAL INTER 
      DATA H/1.43863/ 
      open (7, file='isotope.dat',status='unknown')
      open (6, file='isotope.out',status='unknown')
 59   READ (7,*,END=194) N,TEMPS 
      IF (N.eq. 0) CALL EXIT         
      KT=0
      N6=3*N-6
      N7=3*N-7
      READ (7,110,END=197) name(1) 
      READ (7,*,END=198) (WR1(I),I=1,N6)
      READ (7,110,END=197) name(2)
      READ (7,*,END=198) (WR2(I),I=1,N6)
      READ (7,110,END=197) name(3)
      READ (7,*,END=198) (WT1(I),I=1,N6)
      READ (7,110,END=197) name(4)
      READ (7,*,END=198) (WT2(I),I=1,N6)
   50 CONTINUE
C
C     Hints for obtaining a list of normal mode wavenumbers from a  Gaussian  09 log file.
C     It is possible to write simple macros or programs to parse the output file for this 
C     information.  I find it almost as quick to adopt this simple procedure:
C     1: Open the log file using a good text editor. I use BBedit for  Mac.
C     2: Search for all occurrences of the string Frequencies --
C     3: Copy all occurrences to a blank document
C     4: Using column sensitive (area selection) mode (by holding the  ALT key down in BBedit),
C        drag out a selection box to include up to the text  Frequencies --  and delete the box.
C     5: Delete the top two lines. You should be left with  3N-6/3 lines of wavenumber values.
C     6: Select and copy these into the  isotope.dat file above at the appropriate location
      KT=KT+1 
      IF (TEMPS(KT).EQ.0.0.OR.KT.EQ.11) GO TO 95
      T=TEMPS(KT) 
      CONST=H/T 
      T1=1.0
      T2=1.0
      T3=0. 
      B1=1.0
      B2=1.0
      B3=0. 
      DO 60 I=1,N6
         U1=WR1(I)*CONST
         U2=WR2(I)*CONST
         T1=T1*U2/U1
         T2=T2*(1-EXP(-U1))/(1-EXP(-U2))
         T3=T3+U1-U2
         IF (I.EQ.1) GO TO 60 
         U1=WT1(I)*CONST
         U2=WT2(I)*CONST
         B1=B1*U2/U1
         B2=B2*(1-EXP(-U1))/(1-EXP(-U2))
         B3=B3+U1-U2
   60 CONTINUE
      U1STAR= WT1(1)*CONST/2. 
      U2STAR= WT2(1)*CONST/2. 
      IF(U1STAR.GT.4.0) TUN1 = 0.0
      IF(U1STAR.LT.4.0) TUN1 = (U1STAR)/(SIN(U1STAR)) 
      IF(U2STAR.GT.4.0) TUN2 = 0.0
      IF(U2STAR.LT.4.0) TUN2 = (U2STAR)/(SIN(U2STAR)) 
      TUNCOR=1.0
      IF(TUN2.NE.0.0) TUNCOR= TUN1/TUN2 
      VP=T1/B1
      EXC=T2/B2 
      ZPE=EXP(T3/2.)/EXP(B3/2.) 
      FREQ=WT1(1)/WT2(1)
      HRR=FREQ*VP*EXC*ZPE 
      AKIE = HRR*TUNCOR 
      IF (KT.GT.1) GO TO 80 
      WRITE(6,123) name 
 123  format (a)
      WRITE (6,160) 
      DIFFSR=0. 
      DIFFST=0. 
      WRITE (6,180) 
      DO 70 I=1,N6
      DIFFR = WR1(I)-WR2(I) 
      DIFFT = WT1(I)-WT2(I) 
      DIFFSR = DIFFSR+DIFFR 
      DIFFST = DIFFST+DIFFT 
      IF(I.EQ.1)WRITE (6,171) WR1(I),WR2(I),DIFFR,WT1(I),WT2(I),DIFFT 
   70 IF(I.NE.1)WRITE (6,170) WR1(I),WR2(I),DIFFR,WT1(I),WT2(I),DIFFT 
      WRITE (6,195) DIFFSR, DIFFST
      WRITE (6,190) 
      WRITE (6,200) 
   80 WRITE (6,210) T,FREQ,VP,EXC,ZPE,HRR
C  80 WRITE (6,210) T,FREQ,VP,EXC,ZPE,HRR,TUN1,TUN2,TUNCOR,AKIE 
   90 CONTINUE
      GO TO 50
  95  GO TO 59
 194  CALL EXIT
 198  write (6,199)
 199  format('Error in number of vibrations read in')
      call exit
C 
  100 FORMAT (1H1)
  110 FORMAT (A) 
  120 FORMAT (//25X,8A10) 
  140 FORMAT (I5/10F8.1)
  150 FORMAT (7(F10.5,1X))
  160 FORMAT (  //10X,32HVIBRATIONAL FREQUENCIES  .  .  .,//,10X,2(8HREA
     1CTANT,6X),4X,10HDIFFERENCE,4X,2(8HT. STATE,6X),10HDIFFERENCE,/11X,
     25HLIGHT,9X,5HHEAVY,27X,5HLIGHT,9X,5HHEAVY,/10X,2(7HISOTOPE,7X),18X
     3,2(7HISOTOPE,7X)//) 
  170 FORMAT (11X,3(F7.2,7X),4X,3(F7.2,7X)) 
  171 FORMAT (11X,3(F7.2,7X),4X,3(F8.2,6X)) 
  180 FORMAT (64X,1H*,13X,1H*)
  190 FORMAT (10X,1H*,73H   NORMAL VIBRATIONAL FREQUENCIES ASSOCIATED WI
     1TH THE REACTION COORDINATE,  ////)
  195 FORMAT(13X,"SUM OF ISOTOPIC SHIFTS ", F10.2,36X,F10.2,//) 
  200 FORMAT ( 7X, 40HComponents of partition function ratios:,// 7X, 64
     1HTEMP.    V1(*)/V2(*)      VP          EXC        ZPE         HRR)
C    27X,"U*(1)     U*(2)   TUN COR       KIE"/)
  210 FORMAT ( 7X,F8.4,F11.5,F12.5,1X,F10.5,2X,F10.5,F12.6)
C 210 FORMAT ( 7X,F8.4,F11.5,F12.5,1X,F10.5,2X,F10.5,F12.6,4F10.3)
  197 WRITE(6,2239) 
 2239 FORMAT(' TAPE4 CONTAINS YOUR INPUT DATA FOR SAVING')
      CALL EXIT 
C 
      END 
