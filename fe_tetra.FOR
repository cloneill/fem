C     ******************************************************************
C     *                                                                *
C     * PROGRAM THAT CALCULATES DISPLACEMENTS, REACTIONS, STRESSES     *
C     * AND STRENGTHS IN A FINITE ELEMENT MESH OF LINEAR TETRAHEDRONS  *    
C     *                                                                *
C     * USED REFERENCES FOR IMPLEMENTATION                             *
C     * http://www.yorku.ca/gzhu/ENG3350/Lecture_notes/Chapter_6.pdf   *
C     *                                                                *
C     ******************************************************************
      PROGRAM MAIN
      PARAMETER (NDOFN=3,NNODE=4,NELFN=NDOFN*NNODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ITER
      ALLOCATABLE ::COORD(:,:),FLOAD(:),FRCTN(:),GSTIF(:,:),GSTIF0(:,:),
     +              IDFIX(:,:),INDX(:),LINKN(:,:),MATYP(:),PROPS(:,:)    
     		
      OPEN(15,FILE='INPUT.DAT',STATUS='OLD')
      OPEN(16,FILE='OUTPUT.txt')
C
C---- INPUT CONTROL DATA
C
      CALL DIMEN(NDOFN,NELEM,NEQNS,NOMAT,NPOIN)
C
C --- DYNAMIC ALLOCATE ADDRESS FOR DYNAMIC ARRAYS
C
      ALLOCATE (COORD(NDOFN,NPOIN),FLOAD(NEQNS),GSTIF0(NEQNS,NEQNS),
     +          FRCTN(NEQNS),GSTIF(NEQNS,NEQNS),IDFIX(NDOFN,NPOIN),
     +          INDX(NEQNS),LINKN(NNODE,NELEM),MATYP(NELEM),
     +          PROPS(NOMAT,3))
C
C---- INPUT DATA
C
      CALL INPUT(NDOFN,NELEM,NEQNS,NNODE,NOMAT,NPOIN,COORD,
     +	       FLOAD,IDFIX,LINKN,MATYP,PROPS)
      CLOSE(15)
C
C---- FORM ELEMENT STIFFNESS MATRIX AND ASSEMBLE IT TO GLOBAL MASS MATRIX
C
      CALL FORMK(NDOFN,NELEM,NELFN,NEQNS,NNODE,NOMAT,NPOIN,
     +           COORD,GSTIF,LINKN,MATYP,PROPS)
      GSTIF0=GSTIF
C
C---- APPLY BOUNDARY CONDITIONS
C
      CALL BOUNDARY(NDOFN,NEQNS,NPOIN,GSTIF,IDFIX)
C
C---- SOLVE EQUATIONS: DISPLACEMENT SOLUTION IS STORED IN FLOAD
C
      CALL LUDCMP(GSTIF,NEQNS,INDX,D)
      CALL LUBKSB(GSTIF,NEQNS,INDX,FLOAD)
C
C---  CALCULATE NODAL FORCE FROM DISPLACEMENT AND REACTION FORCES FROM K*U
C
      CALL NSOLU(NDOFN,NELEM,NEQNS,NNODE,NOMAT,NPOIN,COORD,
     +           FLOAD,LINKN,MATYP,PROPS)
      CALL MULTIPLY(GSTIF0,FLOAD,FRCTN,NEQNS,NEQNS,1)
C
C---- OUTPUT DATA
C
      CALL OUTPUT(NDOFN,NEQNS,NPOIN,FLOAD,FRCTN,IDFIX)
      CLOSE(16)
      STOP
      END
C
      SUBROUTINE ASSEMK(IELEM,NDOFN,NELEM,NELFN,NNODE,NEQNS,LINKN,ESTIF,
     +                  GSTIF)                    
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE ASSEMBLES GLOBAL STIFFNESS MATRIX             *
C     *                                                                *
C     *  CALLED BY FORMGK                                              *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ESTIF(NELFN,NELFN),GSTIF(NEQNS,NEQNS),L2G(NELFN),
     +          LINKN(NNODE,NELEM)
C
C ***  SETUP LINK ARRAY BETWEEN LOCAL AND GLOBAL VARIABLES
C
      DO INODE=1,NNODE
        IDOFG=(LINKN(INODE,IELEM)-1)*NDOFN
        IDOFL=(INODE-1)*NDOFN
        DO IDOFN=1,NDOFN
          L2G(IDOFL+IDOFN)=IDOFG+IDOFN
        END DO
      END DO
C
C *** ASSEMBLE WHOLE LOCAL MATRIX INTO GLOBAL MATRIX
C
      DO ICOLE=1,NELFN
        ICOLG=L2G(ICOLE)
        DO IROWE=1,NELFN
          IROWG=L2G(IROWE)
          GSTIF(IROWG,ICOLG)=GSTIF(IROWG,ICOLG)+ESTIF(IROWE,ICOLE)
        END DO
      END DO

      RETURN
      END
C
      SUBROUTINE BOUNDARY(NDOFN,NEQNS,NPOIN,GSTIF,IDFIX)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE DEALS WITH BOUNDARY CONDITIONS AND ECHOS      *
C     *  THEM INTO GLOBAL LOAD ARRAY AND STIFFNESS MATRIX              *
C     *                                                                *
C     *  CALLED BY FRAME                                               *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GSTIF(NEQNS,NEQNS),IDFIX(NDOFN,NPOIN)
      FACTOR=1.D40
C
C---- ECHO BOUNDARY CONDITIONS INTO GLOBAL STIFFNESS MATRIX-------------
C
      DO IPOIN=1,NPOIN
        II=(IPOIN-1)*NDOFN
        DO JDOFN=1,NDOFN
          IF(IDFIX(JDOFN,IPOIN).NE.0) THEN
            KK=II+JDOFN
            GSTIF(KK,KK)=GSTIF(KK,KK)+FACTOR
          END IF
        END DO
      END DO
      RETURN
      END
C
      SUBROUTINE DIMEN(NDOFN,NELEM,NEQNS,NOMAT,NPOIN)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE READS CONTROL DATA DEFINING DIMENSIONS OF     *
C     *  THE PROBLEM.                                                  *
C     *                                                                *
C     *  CALLED BY FRAME                                               *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      READ(15,1000) TITLE
 1000 FORMAT(A80)
      WRITE(16,1100) TITLE
 1100 FORMAT(A80)
C
C---- READ CONTROL DATA ------------------------------------------------
C
      READ(15,*) NPOIN,NELEM,NOMAT
      NEQNS=NDOFN*NPOIN
      WRITE(16,1200) NPOIN,NELEM,NOMAT,NEQNS
 1200 FORMAT(//,5X,'C O N T R O L    D A T A '//
     +       'NUMBER OF NODES IN THE STRUCTURE              =',I10/
     +       'NUMBER OF ELEMENTS IN THE STRUCTURE           =',I10/
     +       'NUMBER OF MATERIALS IN THE STRUCTURE          =',I10/
     +       'NUMBER OF TOTAL EQUATIONS                     =',I10)   
      RETURN
      END
C
      SUBROUTINE FORMK(NDOFN,NELEM,NELFN,NEQNS,NNODE,NOMAT,NPOIN,
     +                 COORD,GSTIF,LINKN,MATYP,PROPS)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE FORMS ELEMENT STIFFNESS MATRIX AND ASSEMBLES  *
C     *  GLOBAL STIFFNESS MATRIX                                       *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION COORD(NDOFN,NPOIN),GSTIF(NEQNS,NEQNS),
     +          LINKN(NNODE,NELEM),MATYP(NELEM),PROPS(NOMAT,3)
      DIMENSION AKGE(NELFN,NELFN)
C
C---- INITIALIZE GLOBAL STIFFNESS MATRIX--------------------------------
C
      GSTIF(:,:)=0.0D0
C
C---- LOOP OVER ELEMENT-------------------------------------------------
C
      DO IELEM=1,NELEM
C
C---- RETRIEVE MATERIAL PROPERTY AND ELEMENT GEOMETRY-------------------
C
        E=PROPS(MATYP(IELEM),1)
        ANU=PROPS(MATYP(IELEM),2)
C
C---- CALCULATE ELEMENT STIFFENSS MATRIX IN GLOBAL REFERENCE SYSTEM-----
C
      CALL FORMKL(AKGE,NDOFN,NPOIN,COORD,NNODE,
     +                  NELEM,LINKN,IELEM,E,ANU)
C
C    ASSEMBLE GLOBAL STIFFNESS MATRIX
C
        CALL ASSEMK(IELEM,NDOFN,NELEM,NELFN,NNODE,NEQNS,LINKN,AKGE,
     +              GSTIF)                    
      END DO
      RETURN
      END
C
      SUBROUTINE INPUT(NDOFN,NELEM,NEQNS,NNODE,NOMAT,NPOIN,COORD,
     +	             FLOAD,IDFIX,LINKN,MATYP,PROPS)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE INPUTS DATA DEFINING GEOMETRY, BOUNDARY       *
C     *  CONDITIONS, MATERIAL PROPERTY AND APPLIED LOADING.            *
C     *                                                                *
C     *  CALLED BY FRAME                                               *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(NDOFN,NPOIN),FLOAD(NEQNS),
     +          IDFIX(NDOFN,NPOIN),LINKN(NNODE,NELEM),MATYP(NELEM),
     +          PROPS(NOMAT,3)
C                                                                       
C---- READ NODAL POINT DATA--------------------------------
C                                                                       
      WRITE (16,1000)
 1000 FORMAT(//'NODAL POINT INPUT DATA'/
     +       'NODE',3X,'BOUNDARY CONDITION CODES',5X,
     +       'NODAL POINT COORDINATES',15X,
     +        /'NUMBER',10X,'X',4X,'Y',4X,'Z',17X,'-X-',12X,'-Y-',
     +        12X,'-Z-',12X)
      DO J=1,NPOIN
        READ(15,*)     N,(IDFIX(I,N),I=1,NDOFN),(COORD(I,N),I=1,NDOFN)
        WRITE(16,1100) N,(IDFIX(I,N),I=1,NDOFN),(COORD(I,N),I=1,NDOFN)
      END DO
 1100 FORMAT (I4,8X,3I5,7X,3(2X,F13.3))
C
C---- READ ELEMENT DATA------------------------------------
C
      WRITE(16,1800)
 1800 FORMAT(//'F R A M E    E L E M E N T    D A T A'//'ELE. NO.',5X,
     +       'N--I',5X,'N--J',5X,'N--K',5X,'N--L',5X,'MAT. NO.',10X)
      DO I=1,NELEM
        READ(15,*) NUM,(LINKN(INODE,NUM),INODE=1,NNODE),
     +             MATYP(NUM)
        WRITE(16,1900) NUM,(LINKN(INODE,NUM),INODE=1,NNODE),
     +                 MATYP(NUM)
      END DO
C 1900 FORMAT(I4,9X,I4,5X,I4,5X,I4,5X,I4,5X,I4,4X,1PE18.4,I5)
 1900 FORMAT(I4,9X,I4,5X,I4,5X,I4,5X,I4,5X,I4,4X)
C
C---- READ MATERIAL PROPERTY PLUS NEWMARK INTEGRATION PARAMETERS -------
C                                                                       
      WRITE(16,2100)
 2100 FORMAT(//'MATERIAL PROPERTY TABLE   '//
     +       'NUMBER',13X,'YOUNGS MODULAR',13X,'POISSON RATIO',13X,
     +        'MATERIAL DENSITY')
      DO I=1,NOMAT
        READ(15,*) NMAT,(PROPS(NMAT,J),J=1,3)
        WRITE(16,2200) NMAT,(PROPS(NMAT,J),J=1,3)
      END DO
 2200 FORMAT(I4,12X,1PE15.4,12X,1PE15.4,12X,1PE15.4)
C
C---- READ NODAL POINT LOADING------------------------------------------
C
      FLOAD(:)=0.0
      WRITE(16,2350)
 2350 FORMAT(/'NODE   ',16X,'INPUT CONCENTRATE NODAL LOAD DATA'/
     +       'NUMBER',17X,'-Fx-',15X,'-Fy-'15X,'-Fz-')
      READ(15,*) NLOAD
      DO 800 I=1,NLOAD
		READ(15,*)     NP,FLOAD((NP-1)*NDOFN+1),FLOAD((NP-1)*NDOFN+2),
     +       FLOAD((NP-1)*NDOFN+3)
		WRITE(16,2400) NP,FLOAD((NP-1)*NDOFN+1),FLOAD((NP-1)*NDOFN+2),
     +       FLOAD((NP-1)*NDOFN+3)
  800 CONTINUE
 2400 FORMAT(I4,5X,1PE21.4,1PE19.4,1PE18.4)
      RETURN
      END
C
      SUBROUTINE ludcmp(a,n,indx,d)
C     ******************************************************************
C     *  THIS SUBROUTINE SOLVES THE LINEAR SIMULTANEOUS EQUATIONS BY   *
C     *  L*D*TP(L) METHOD                                              *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INDX(N)
      DOUBLE PRECISION d,a(n,n),TINY
      PARAMETER (TINY=1.0d-30)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(n)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C
      SUBROUTINE lubksb(a,n,indx,b)
C     *********************************************************
C     *  THIS SUBROUTINE DOES THE BACK-SUBSTITUTION FOR THE   *
C     *  L*D*TP(L) METHOD                                     *
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION INDX(N)
      DOUBLE PRECISION a(n,n),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C
      SUBROUTINE OUTPUT(NDOFN,NEQNS,NPOIN,GDISP,FRCTN,IDFIX)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE OUTPUTS NODAL DISPLACEMENT AND MEMBER FORCE   *
C     *                                                                *
C     *  CALLED BY FRAME                                               *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION GDISP(NEQNS),FRCTN(NEQNS),IDFIX(NDOFN,NPOIN)
C
C---- OUTPUT NODAL DISPLACEMENT-----------------------------------------
C
      WRITE(16,1000)
 1000 FORMAT(//'O U T P U T    D A T A'//'NODAL DISPLACEMENTS'/
     +       'NODE NO.',10X,'-U-',10X,'-V-',10X,'-W-')
      DO 100 IPOIN=1,NPOIN
		K=(IPOIN-1)*NDOFN
		WRITE(16,1100) IPOIN,(GDISP(K+I),I=1,NDOFN)
  100 CONTINUE
 1100 FORMAT(I4,6X,1PE13.4,1PE13.4,1PE13.4)
C
C---- OUTPUT REACTION FORCES-----------------------------------------
C
      WRITE(16,2000)
 2000 FORMAT(/'REACTION FORCES'/
     +       'NODE NO.',9X,'-RX-',9X,'-RY-',9X,'-RZ-')
      DO IPOIN=1,NPOIN
        K=0
        DO IDOFN=1,NDOFN
          K=K+IDFIX(IDOFN,IPOIN)
        END DO
C        IF(K.NE.0) THEN
		K=(IPOIN-1)*NDOFN
		WRITE(16,1100) IPOIN,(FRCTN(K+I),I=1,NDOFN)
C        END IF
      END DO
C      RETURN
      END
C
      SUBROUTINE NSOLU(NDOFN,NELEM,NEQNS,NNODE,NOMAT,NPOIN,COORD,
     +                 GDISP,LINKN,MATYP,PROPS)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE CALCULATES REACTION FORCES                    *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION COORD(NDOFN,NPOIN),PROPS(NOMAT,3),
     +          GDISP(NEQNS),LINKN(NNODE,NELEM),MATYP(NELEM),
     +          BLMTX(6,12),DMATX(6,6),VMAT(4,4),
     +          DISPL(12),STRAN(6),STRAS(6)
C
C---- LOOP OVER ELEMENT-------------------------------------------------
C
      DO IELEM=1,NELEM
C
C---- RETRIEVE MATERIAL PROPERTY AND ELEMENT GEOMETRY-------------------
C
        E=PROPS(MATYP(IELEM),1)
        ANU=PROPS(MATYP(IELEM),2)
C
C	RETRIEVE NODAL DISPLACEMENT ASSOCIATED WITH ELEMENT NODES
C
        NP1=LINKN(1,IELEM)
        NP2=LINKN(2,IELEM)
        NP3=LINKN(3,IELEM)
        NP4=LINKN(4,IELEM)

        DISPL(1)=GDISP(NP1*3-2)
        DISPL(2)=GDISP(NP1*3-1)
        DISPL(3)=GDISP(NP1*3)
        DISPL(4)=GDISP(NP2*3-2)
        DISPL(5)=GDISP(NP2*3-1)
        DISPL(6)=GDISP(NP2*3)
        DISPL(7)=GDISP(NP3*3-2)
        DISPL(8)=GDISP(NP3*3-1)
        DISPL(9)=GDISP(NP3*3)
        DISPL(10)=GDISP(NP4*3-2)
        DISPL(11)=GDISP(NP4*3-1)
        DISPL(12)=GDISP(NP4*3)
C
C     CALCULATE NODAL SOLUTIONS AT ELEMENT LEVEL
C
		VMAT(1,1)=1.0d0
		VMAT(1,2)=1.0d0
		VMAT(1,3)=1.0d0
		VMAT(1,4)=1.0d0
		
		VMAT(2,1)=COORD(1,LINKN(1,IELEM))
		VMAT(2,2)=COORD(1,LINKN(2,IELEM))
		VMAT(2,3)=COORD(1,LINKN(3,IELEM))
		VMAT(2,4)=COORD(1,LINKN(4,IELEM))

		VMAT(3,1)=COORD(2,LINKN(1,IELEM))
		VMAT(3,2)=COORD(2,LINKN(2,IELEM))
		VMAT(3,3)=COORD(2,LINKN(3,IELEM))
		VMAT(3,4)=COORD(2,LINKN(4,IELEM))

		VMAT(4,1)=COORD(3,LINKN(1,IELEM))
		VMAT(4,2)=COORD(3,LINKN(2,IELEM))
		VMAT(4,3)=COORD(3,LINKN(3,IELEM))
		VMAT(4,4)=COORD(3,LINKN(4,IELEM))

        CALL M44DET(VMAT,V6)
		
        VOLU=V6/6.0d0
        CALL FORMB(NDOFN,NPOIN,COORD,NNODE,NELEM,LINKN,IELEM,ANU,
     +             BLMTX,VOLU)
        CALL MULTIPLY(BLMTX,DISPL,STRAN,6,12,1)
        CALL FORMD(ANU,DMATX,E)
        CALL MULTIPLY(DMATX,STRAN,STRAS,6,6,1)
      END DO

      RETURN
      END
C
      SUBROUTINE M44DET(A,DET)
C     ************************************************
C     * Compute the determinant of a 4x4 matrix      *
C     ************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(4,4)

      DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*
     +       (A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-
     +       A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*
     +       A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+
     +       A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*
     +       (A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-
     +       A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-
     +       A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+
     +       A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*
     +       A(4,2)-A(3,2)*A(4,1)))

      RETURN
      END
C
      SUBROUTINE TMULTIPLY(A,B,C,L,M,N)
C     ************************************************
C     *  THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT *
C     *                 C=TRANSPOSE(A)*B             *
C     ************************************************
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N
      DOUBLE PRECISION A,B,C
      REAL*16 S
      DIMENSION A(M,L),B(M,N),C(L,N)
      DO I=1,L
        DO J=1,N
          C(I,J)=0.D0
          S=0.D0
          DO K=1,M
            S=S+A(K,I)*B(K,J)
          END DO
        C(I,J)=S
        END DO
      END DO
      RETURN
      END
C
      SUBROUTINE MULTIPLY(A,B,C,L,M,N)
C     ************************************************
C     *  THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT *
C     *                 C=A*B                        *
C     ************************************************
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N
      DOUBLE PRECISION A,B,C
      DIMENSION A(L,M),B(M,N),C(L,N)
      DO I=1,L
        DO J=1,N
          C(I,J)=0.D0
C          S=0.D0
          DO K=1,M
			C(I,J)=C(I,J)+A(I,K)*B(K,J)
C			PRINT *,A(I,K)
C			PRINT *,B(K,J)
C			PRINT *,T
C			PRINT *,S
C 			PRINT *,C(I,J)
          END DO
C        C(I,J)=S
        END DO
      END DO
      RETURN
      END
C
      SUBROUTINE FORMD(ANU,DMATX,E)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE DEALS WITH DETERMINING THE D MATIRX            *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION DMATX(6,6)
C
      DMATX=0.0

      CONST=E/((1.0+ANU)*(1.0-2.0*ANU))
      DMATX(1,1)=CONST*(1.0-ANU)
      DMATX(1,2)=CONST*ANU
      DMATX(1,3)=CONST*ANU
	  
      DMATX(2,1)=CONST*ANU
      DMATX(2,2)=CONST*(1.0-ANU)
      DMATX(2,3)=CONST*ANU
	  
      DMATX(3,1)=CONST*ANU
      DMATX(3,2)=CONST*ANU
      DMATX(3,3)=CONST*(1.0-ANU)
	  
      DMATX(4,4)=CONST*0.5*(1.0-2*ANU)
      DMATX(5,5)=CONST*0.5*(1.0-2*ANU)
      DMATX(6,6)=CONST*0.5*(1.0-2*ANU)
C
      RETURN
      END
C
      SUBROUTINE FORMB(NDOFN,NPOIN,COORD,NNODE,NELEM,LINKN,IELEM,ANU,
     +                 BLMTX,VOLU)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE DEALS WITH DETERMINING THE LINEAR B MATRIX    *
C     *                                                                *
C     ******************************************************************	
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION BLMTX(6,12),COORD(NDOFN,NPOIN),LINKN(NNODE,NELEM),
     +          A(4),B(4),C(4)
C
        A(1)=-(COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(2,IELEM))+
     +       COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(4,IELEM)))

        A(2)=COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(4,IELEM))
	
        A(3)=-(COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(2,LINKN(4,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(4,IELEM)))
	
        A(4)=COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(3,IELEM))+
     +       COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(2,LINKN(2,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(2,LINKN(3,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(2,LINKN(1,IELEM))*COORD(3,LINKN(3,IELEM))

        B(1)=COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(2,IELEM))+
     +       COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(4,IELEM))

        B(2)=-(COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(3,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(4,IELEM)))

        B(3)=COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(4,IELEM))

        B(4)=-(COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(3,IELEM))+
     +       COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(3,LINKN(1,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(3,LINKN(2,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(3,LINKN(3,IELEM)))

        C(1)=-(COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(2,IELEM))+
     +       COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(3,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(2,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(3,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(4,IELEM)))

        C(2)=COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(3,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(1,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(3,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(4,IELEM))

        C(3)=-(COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(4,IELEM))+
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(2,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(1,IELEM))-
     +       COORD(1,LINKN(4,IELEM))*COORD(2,LINKN(2,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(4,IELEM)))

        C(4)=COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(3,IELEM))+
     +       COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(1,IELEM))+
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(2,IELEM))-
     +       COORD(1,LINKN(2,IELEM))*COORD(2,LINKN(1,IELEM))-
     +       COORD(1,LINKN(3,IELEM))*COORD(2,LINKN(2,IELEM))-
     +       COORD(1,LINKN(1,IELEM))*COORD(2,LINKN(3,IELEM))
C
C	FORMING LINEAR B MATRIX
C
      BLMTX=0.0
C	
      DO I=1,4
		BLMTX(1,1+(I-1)*3)=A(I)/(6*VOLU)
		BLMTX(2,2+(I-1)*3)=B(I)/(6*VOLU)
		BLMTX(3,3+(I-1)*3)=C(I)/(6*VOLU)
	
		BLMTX(4,1+(I-1)*3)=B(I)/(6*VOLU)
		BLMTX(4,2+(I-1)*3)=A(I)/(6*VOLU)

		BLMTX(5,2+(I-1)*3)=C(I)/(6*VOLU)
		BLMTX(5,3+(I-1)*3)=B(I)/(6*VOLU)

		BLMTX(6,1+(I-1)*3)=C(I)/(6*VOLU)
		BLMTX(6,3+(I-1)*3)=A(I)/(6*VOLU)
      END DO
C
      RETURN
      END
C
      SUBROUTINE FORMKL(ESTIF,NDOFN,NPOIN,COORD,NNODE,NELEM,LINKN,
     +                IELEM,E,ANU)
C     ******************************************************************
C     *                                                                *
C     *  THIS SUBROUTINE FORMS ELEMENT STIFFNESS MATRIX                *
C     *                                                                *
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ESTIF(12,12),BLMTX(6,12),DMATX(6,6),VMAT(4,4),
     +          COORD(NDOFN,NPOIN),LINKN(NNODE,NELEM)
C
      CALL FORMD(ANU,DMATX,E)
		
      VMAT(1,1)=1.0d0
      VMAT(1,2)=1.0d0
      VMAT(1,3)=1.0d0
      VMAT(1,4)=1.0d0
		
      VMAT(2,1)=COORD(1,LINKN(1,IELEM))
      VMAT(2,2)=COORD(1,LINKN(2,IELEM))
      VMAT(2,3)=COORD(1,LINKN(3,IELEM))
      VMAT(2,4)=COORD(1,LINKN(4,IELEM))

      VMAT(3,1)=COORD(2,LINKN(1,IELEM))
      VMAT(3,2)=COORD(2,LINKN(2,IELEM))
      VMAT(3,3)=COORD(2,LINKN(3,IELEM))
      VMAT(3,4)=COORD(2,LINKN(4,IELEM))

      VMAT(4,1)=COORD(3,LINKN(1,IELEM))
      VMAT(4,2)=COORD(3,LINKN(2,IELEM))
      VMAT(4,3)=COORD(3,LINKN(3,IELEM))
      VMAT(4,4)=COORD(3,LINKN(4,IELEM))

      CALL M44DET(VMAT,V6)
		
      VOLU=V6/6.0d0

      CALL FORMB(NDOFN,NPOIN,COORD,NNODE,NELEM,LINKN,IELEM,ANU,
     +           BLMTX,VOLU)
	  CALL BTAB(ESTIF,DMATX,BLMTX,6,12)
	  ESTIF=ESTIF*VOLU
C	
      RETURN
      END
C
      SUBROUTINE BTAB(C,A,B,N,M)
C     **************************************************
C     *   THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT  *
C     *   C=TR(B)*A*B  WHERE  A  IS A SQUARE MATRIX    *
C     *   A(N,N), B(N,M), C(M,M), BB(N,M)              *
C     *   BB:AUXILIARY MATRIX                          *
C     **************************************************
      IMPLICIT NONE
      INTEGER M,N,J,K,KK
      DOUBLE PRECISION A,B,C,BB
      DIMENSION C(M,M),A(N,N),B(N,M),BB(N,M)
      DO J=1,N
        DO K=1,M
          BB(J,K)=0.0D0
          DO KK=1,N
            BB(J,K)=BB(J,K)+A(J,KK)*B(KK,K)
          END DO
        END DO
      END DO
      DO KK=1,M
        DO J=1,M
          C(KK,J)=0.0D0
          DO K=1,N
            C(KK,J)=C(KK,J)+B(K,KK)*BB(K,J)
          END DO
        END DO
      END DO
      RETURN
      END
C
