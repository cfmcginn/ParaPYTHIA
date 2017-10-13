C...A simple skeleton program, illustrating a typical Pythia run:
C...Z0 production at LEP 1. 
C...Toy task: compare multiplicity distribution with matrix elements
C...and with parton showers (using same fragmentation parameters).

C-----------------------------------------------------------------

C...Preamble: declarations.
 
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
      INTEGER*4 n, time
C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYDATR/MRPY(6),RRPY(100)


C-----------------------------------------------------------------

C...First section: initialization.
 
C...Main parameters of run: c.m. energy and number of events.
      ECM=5020.0D0
      NEV=10000
      CKIN(3)=500D0
      CKIN(4)=-1D0
      MSTP(81)=0
C...Select gamma*/Z0 production process.
      MSEL=1
C...      MSUB(11)=1
C...      MSUB(12)=1
C...      MSUB(13)=1
C...      MSUB(28)=1
C...      MSUB(53)=1
C...      MSUB(68)=1
C...      MSUB(96)=1

C...      MSEL=0

C...Only allow Z0 decay to quarks (i.e. no leptonic final states).
C...      DO 100 IDC=MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
C...        IF(IABS(KFDP(IDC,1)).GE.6) MDME(IDC,1)=MIN(0,MDME(IDC,1))
C..  100 CONTINUE

C...Initialize.
      CALL PYINIT('CMS','p+','p+',ECM)

C...Check that Z0 decay channels set correctly.
C...      CALL PYSTAT(2)

C-----------------------------------------------------------------

C...Second section: event loop.
      n = time()
      MRPY(1)=n

C...Outer loop over ME and PS options.
      DO 300 ICA=1
        MSTP(81)=0
         
C...Begin event loop.
        DO 200 IEV=1,NEV
          CALL PYEVNT

C...List first few events.
C...          IF(IEV.LE.2) CALL PYLIST(1)
C...          print *, "pre-event: pthat"
C...          print *, PARI(17)
C...          CALL PYLIST(1)

C...Extract and fill event properties.
C...          CALL PYEDIT(3)
C...          CALL PYFILL(ICA,DBLE(N),1D0)

C...End event loop.
  200   CONTINUE

C...End outer loop.
  300 CONTINUE

C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
C...      CALL PYHIST

      END
