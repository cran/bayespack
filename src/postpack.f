*********************************************************
* Orignal BAYESPACK source code (author: Alan Genz)
*********************************************************
*
*
*    This file contains a collection of Log Posterior example densities.
*    The densities below are given in alphabetical order:
*
*    BWLPST: Birthweight Data, 7-d 
*    BXLPST: Bio Oxygen Data, 2-d 
*    CTLPST: Contingency Table Data, 9-d 
*    EXLPST: Extreme Value Distribution, 5-d
*    FLLPST: Fuller Data, 3-d
*    HRLPST: Stanford Heart Transplant Data, 3-d 
*    K3LPST: Econometric Data, 3-d
*    K5LPST: Econometric Data, 5-d
*    K8LPST: Econometric Data, 8-d
*    LBLPST: Lubricant Data, 10-d
*    LGLPST: Multivariate Logistic Distribution, 5-d 
*    LMLPST: Linear Model Data, 10-d
*    MTLPST: Motorette Data, 3-d
*    NLLPST: Nonlinear Regression Model, 3-d 
*    NRLPST: Multivariate Normal Distribution, 6-d 
*    PGLPST: Porgi Data, 14-d
*    PHLPST: Photocarcinogen Data, 5-d
*    PRLPST: Proportional Hazards Data, 7-d
*    PSLPST: Pearson Type IV, 1-d
*    RDLPST: Radiotherapy Data, 2-d
*    TNLPST: Tornado Data, 11-d
*
*
      SUBROUTINE BWLPST(Z, OUT)
*
*     Log Posterior for Birthweight Data example from
*     Venables, W. N., and Ripley, B. D. (1994),
*     "Modern Applied Statistics with S-Plus", Springer-Verlag, New York, 
*     pp. 193-195. 
*     The model implemented is a seven variable model,
*     and only 1/IS of the data is used.
*
*     Approximate mode full data model (IS = 1):
*      DOUBLE PRECISION MODE(7)
*      DATA MODE/  
*     &   1.49955,-0.05232,-0.01404, 0.48526, 1.31542, 1.86334, 0.68343/
*
*     Approximate mode reduced data model (IS = 5):
*      DOUBLE PRECISION MODE(7)
*      DATA MODE/  
*     &  -1.80561, 0.01919,-0.00560, 1.89159, 2.65261, 0.49013,-1.05303/
*
      DOUBLE PRECISION OUT
      INTEGER I, N, IC, IS
      PARAMETER ( N = 189 )
      DOUBLE PRECISION FVALUE, Z(7), P
      INTEGER low(N), age(N), lwt(N), race(N), smoke(N), ptl(N), 
     &        ht(N), ui(N), ftv(N), bwt(N), ptd(N), ftv1(N), ftv2(N)
      SAVE low, age, lwt, race, smoke, ptl, ht, ui, ftv, bwt, 
     &     ptd, ftv1, ftv2, IC
      COMMON /STPBLK/ IS
      DATA IC / 0 /
      DATA (low(I), age(I), lwt(I), race(I), smoke(I), ptl(I), 
     &      ht(I), ui(I), ftv(I), bwt(I), I = 1,40)/
     & 0,19,182,2,0,0,0,1,0,2523,  0,33,155,3,0,0,0,0,3,2551,
     & 0,20,105,1,1,0,0,0,1,2557,  0,21,108,1,1,0,0,1,2,2594,
     & 0,18,107,1,1,0,0,1,0,2600,  0,21,124,3,0,0,0,0,0,2622,
     & 0,22,118,1,0,0,0,0,1,2637,  0,17,103,3,0,0,0,0,1,2637,
     & 0,29,123,1,1,0,0,0,1,2663,  0,26,113,1,1,0,0,0,0,2665,
     & 0,19, 95,3,0,0,0,0,0,2722,  0,19,150,3,0,0,0,0,1,2733,
     & 0,22, 95,3,0,0,1,0,0,2751,  0,30,107,3,0,1,0,1,2,2750,
     & 0,18,100,1,1,0,0,0,0,2769,  0,18,100,1,1,0,0,0,0,2769,
     & 0,15, 98,2,0,0,0,0,0,2778,  0,25,118,1,1,0,0,0,3,2782,
     & 0,20,120,3,0,0,0,1,0,2807,  0,28,120,1,1,0,0,0,1,2821,
     & 0,32,121,3,0,0,0,0,2,2835,  0,31,100,1,0,0,0,1,3,2835,
     & 0,36,202,1,0,0,0,0,1,2836,  0,28,120,3,0,0,0,0,0,2863,
     & 0,25,120,3,0,0,0,1,2,2877,  0,28,167,1,0,0,0,0,0,2877,
     & 0,17,122,1,1,0,0,0,0,2906,  0,29,150,1,0,0,0,0,2,2920,
     & 0,26,168,2,1,0,0,0,0,2920,  0,17,113,2,0,0,0,0,1,2920,
     & 0,17,113,2,0,0,0,0,1,2920,  0,24, 90,1,1,1,0,0,1,2948,
     & 0,35,121,2,1,1,0,0,1,2948,  0,25,155,1,0,0,0,0,1,2977,
     & 0,25,125,2,0,0,0,0,0,2977,  0,29,140,1,1,0,0,0,2,2977,
     & 0,19,138,1,1,0,0,0,2,2977,  0,27,124,1,1,0,0,0,0,2922,
     & 0,31,215,1,1,0,0,0,2,3005,  0,33,109,1,1,0,0,0,1,3033/
      DATA (low(I), age(I), lwt(I), race(I), smoke(I), ptl(I), 
     &      ht(I), ui(I), ftv(I), bwt(I), I = 41,80)/
     & 0,21,185,2,1,0,0,0,2,3042,  0,19,189,1,0,0,0,0,2,3062,
     & 0,23,130,2,0,0,0,0,1,3062,  0,21,160,1,0,0,0,0,0,3062,
     & 0,18, 90,1,1,0,0,1,0,3062,  0,18, 90,1,1,0,0,1,0,3062,
     & 0,32,132,1,0,0,0,0,4,3080,  0,19,132,3,0,0,0,0,0,3090,
     & 0,24,115,1,0,0,0,0,2,3090,  0,22, 85,3,1,0,0,0,0,3090,
     & 0,22,120,1,0,0,1,0,1,3100,  0,23,128,3,0,0,0,0,0,3104,
     & 0,22,130,1,1,0,0,0,0,3132,  0,30, 95,1,1,0,0,0,2,3147,
     & 0,19,115,3,0,0,0,0,0,3175,  0,16,110,3,0,0,0,0,0,3175,
     & 0,21,110,3,1,0,0,1,0,3203,  0,30,153,3,0,0,0,0,0,3203,
     & 0,20,103,3,0,0,0,0,0,3203,  0,17,119,3,0,0,0,0,0,3225,
     & 0,17,119,3,0,0,0,0,0,3225,  0,23,119,3,0,0,0,0,2,3232,
     & 0,24,110,3,0,0,0,0,0,3232,  0,28,140,1,0,0,0,0,0,3234,
     & 0,26,133,3,1,2,0,0,0,3260,  0,20,169,3,0,1,0,1,1,3274,
     & 0,24,115,3,0,0,0,0,2,3274,  0,28,250,3,1,0,0,0,6,3303,
     & 0,20,141,1,0,2,0,1,1,3317,  0,22,158,2,0,1,0,0,2,3317,
     & 0,22,112,1,1,2,0,0,0,3317,  0,31,150,3,1,0,0,0,2,3321,
     & 0,23,115,3,1,0,0,0,1,3331,  0,16,112,2,0,0,0,0,0,3374,
     & 0,16,135,1,1,0,0,0,0,3374,  0,18,229,2,0,0,0,0,0,3402,
     & 0,25,140,1,0,0,0,0,1,3416,  0,32,134,1,1,1,0,0,4,3430,
     & 0,20,121,2,1,0,0,0,0,3444,  0,23,190,1,0,0,0,0,0,3459/
      DATA (low(I), age(I), lwt(I), race(I), smoke(I), ptl(I), 
     &      ht(I), ui(I), ftv(I), bwt(I), I = 81,120)/
     & 0,22,131,1,0,0,0,0,1,3460,  0,32,170,1,0,0,0,0,0,3473,
     & 0,30,110,3,0,0,0,0,0,3544,  0,20,127,3,0,0,0,0,0,3487,
     & 0,23,123,3,0,0,0,0,0,3544,  0,17,120,3,1,0,0,0,0,3572,
     & 0,19,105,3,0,0,0,0,0,3572,  0,23,130,1,0,0,0,0,0,3586,
     & 0,36,175,1,0,0,0,0,0,3600,  0,22,125,1,0,0,0,0,1,3614,
     & 0,24,133,1,0,0,0,0,0,3614,  0,21,134,3,0,0,0,0,2,3629,
     & 0,19,235,1,1,0,1,0,0,3629,  0,25, 95,1,1,3,0,1,0,3637,
     & 0,16,135,1,1,0,0,0,0,3643,  0,29,135,1,0,0,0,0,1,3651,
     & 0,29,154,1,0,0,0,0,1,3651,  0,19,147,1,1,0,0,0,0,3651,
     & 0,19,147,1,1,0,0,0,0,3651,  0,30,137,1,0,0,0,0,1,3699,
     & 0,24,110,1,0,0,0,0,1,3728,  0,19,184,1,1,0,1,0,0,3756,
     & 0,24,110,3,0,1,0,0,0,3770,  0,23,110,1,0,0,0,0,1,3770,
     & 0,20,120,3,0,0,0,0,0,3770,  0,25,241,2,0,0,1,0,0,3790,
     & 0,30,112,1,0,0,0,0,1,3799,  0,22,169,1,0,0,0,0,0,3827,
     & 0,18,120,1,1,0,0,0,2,3856,  0,16,170,2,0,0,0,0,4,3860,
     & 0,32,186,1,0,0,0,0,2,3860,  0,18,120,3,0,0,0,0,1,3884,
     & 0,29,130,1,1,0,0,0,2,3884,  0,33,117,1,0,0,0,1,1,3912,
     & 0,20,170,1,1,0,0,0,0,3940,  0,28,134,3,0,0,0,0,1,3941,
     & 0,14,135,1,0,0,0,0,0,3941,  0,28,130,3,0,0,0,0,0,3969,
     & 0,25,120,1,0,0,0,0,2,3983,  0,16, 95,3,0,0,0,0,1,3997/
      DATA (low(I), age(I), lwt(I), race(I), smoke(I), ptl(I), 
     &      ht(I), ui(I), ftv(I), bwt(I), I = 121,160)/
     & 0,20,158,1,0,0,0,0,1,3997,  0,26,160,3,0,0,0,0,0,4054,
     & 0,21,115,1,0,0,0,0,1,4054,  0,22,129,1,0,0,0,0,0,4111,
     & 0,25,130,1,0,0,0,0,2,4153,  0,31,120,1,0,0,0,0,2,4167,
     & 0,35,170,1,0,1,0,0,1,4174,  0,19,120,1,1,0,0,0,0,4238,
     & 0,24,116,1,0,0,0,0,1,4593,  0,45,123,1,0,0,0,0,1,4990,
     & 1,28,120,3,1,1,0,1,0, 709,  1,29,130,1,0,0,0,1,2,1021,
     & 1,34,187,2,1,0,1,0,0,1135,  1,25,105,3,0,1,1,0,0,1330,
     & 1,25, 85,3,0,0,0,1,0,1474,  1,27,150,3,0,0,0,0,0,1588,
     & 1,23, 97,3,0,0,0,1,1,1588,  1,24,128,2,0,1,0,0,1,1701,
     & 1,24,132,3,0,0,1,0,0,1729,  1,21,165,1,1,0,1,0,1,1790,
     & 1,32,105,1,1,0,0,0,0,1818,  1,19, 91,1,1,2,0,1,0,1885,
     & 1,25,115,3,0,0,0,0,0,1893,  1,16,130,3,0,0,0,0,1,1899,
     & 1,25, 92,1,1,0,0,0,0,1928,  1,20,150,1,1,0,0,0,2,1928,
     & 1,21,200,2,0,0,0,1,2,1928,  1,24,155,1,1,1,0,0,0,1936,
     & 1,21,103,3,0,0,0,0,0,1970,  1,20,125,3,0,0,0,1,0,2055,
     & 1,25, 89,3,0,2,0,0,1,2055,  1,19,102,1,0,0,0,0,2,2082,
     & 1,19,112,1,1,0,0,1,0,2084,  1,26,117,1,1,1,0,0,0,2084,
     & 1,24,138,1,0,0,0,0,0,2100,  1,17,130,3,1,1,0,1,0,2125,
     & 1,20,120,2,1,0,0,0,3,2126,  1,22,130,1,1,1,0,1,1,2187,
     & 1,27,130,2,0,0,0,1,0,2187,  1,20, 80,3,1,0,0,1,0,2211/
      DATA (low(I), age(I), lwt(I), race(I), smoke(I), ptl(I), 
     &      ht(I), ui(I), ftv(I), bwt(I), I = 161,N)/
     & 1,17,110,1,1,0,0,0,0,2225,  1,25,105,3,0,1,0,0,1,2240,
     & 1,20,109,3,0,0,0,0,0,2240,  1,18,148,3,0,0,0,0,0,2282,
     & 1,18,110,2,1,1,0,0,0,2296,  1,20,121,1,1,1,0,1,0,2296,
     & 1,21,100,3,0,1,0,0,4,2301,  1,26, 96,3,0,0,0,0,0,2325,
     & 1,31,102,1,1,1,0,0,1,2353,  1,15,110,1,0,0,0,0,0,2353,
     & 1,23,187,2,1,0,0,0,1,2367,  1,20,122,2,1,0,0,0,0,2381,
     & 1,24,105,2,1,0,0,0,0,2381,  1,15,115,3,0,0,0,1,0,2381,
     & 1,23,120,3,0,0,0,0,0,2410,  1,30,142,1,1,1,0,0,0,2410,
     & 1,22,130,1,1,0,0,0,1,2410,  1,17,120,1,1,0,0,0,3,2414,
     & 1,23,110,1,1,1,0,0,0,2424,  1,17,120,2,0,0,0,0,2,2438,
     & 1,26,154,3,0,1,1,0,1,2442,  1,20,105,3,0,0,0,0,3,2450,
     & 1,26,190,1,1,0,0,0,0,2466,  1,14,101,3,1,1,0,0,0,2466,
     & 1,28, 95,1,1,0,0,0,2,2466,  1,14,100,3,0,0,0,0,2,2495,
     & 1,23, 94,3,1,0,0,0,0,2495,  1,17,142,2,0,0,1,0,0,2495,
     & 1,21,130,1,1,0,1,0,3,2495/ 
      IF ( IC .EQ. 0 ) THEN
         DO I = 1, N
            ptd(I) = MIN( ptl(I), 1 )
            ftv1(I) = 1 - MIN( ABS( ftv(I) - 1 ), 1 )
            ftv2(I) = MIN( ftv(I)/2, 1 )
         END DO
         IC = 1
      END IF      
      FVALUE = 0
      IS = 1
      DO I = 1, N
         P = Z(1) + age(I)*Z(2) + lwt(I)*Z(3) + smoke(I)*Z(4) 
     &            + ptd(I)*Z(5) +  ht(I)*Z(6) +    ui(I)*Z(7) 
         FVALUE = FVALUE + low(I)*P
         IF ( ABS(P) .LT. 50 ) THEN
            FVALUE = FVALUE - LOG( 1 + EXP(P) )
         ELSE IF ( P .GT. 50 ) THEN
            FVALUE = FVALUE - P
         END IF
      END DO
      OUT = FVALUE
      END
*
      SUBROUTINE BXLPST(Z, OUT)
*
*     To compute log posterior for Bio Oxygen Data problem taken from
*     Bates, D. M., and Watts, D. G. (1988),
*     Nonlinear Regression and Its Applications, Wiley, New York,
*     pp. 40-43, 270. 
*     
*     Original variables are transformed using
*     \theta_1  = 60exp(Z(1))/( 1 + exp(Z(1) )
*     \theta_2  = 60exp(Z(2))/( 1 + exp(Z(2) )
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(2)
*      DATA MODE/ -0.7861D0, -2.25D0 /
*
      DOUBLE PRECISION OUT
      INTEGER I, N
      PARAMETER ( N = 6 )
      INTEGER X(N)
      DOUBLE PRECISION Z(2), E1, E2, THETA(2), SUM, Y(N)
      SAVE Y, X
      DATA ( X(I), Y(I), I = 1, N ) /
     &     1, 8.3D0, 2, 10.3D0, 3, 19D0, 4, 16D0, 5, 15.6D0, 7, 19.8D0/
      E1 = EXP( Z(1) )
      E2 = EXP( Z(2) )
      THETA(1) = 60*E1/( 1 + E1 )
      THETA(2) =  6*E2/( 1 + E2 )
      OUT = LOG( THETA(1)*THETA(2)/( ( 1 + E1 )*( 1 + E2 ) ) )
      SUM = 0
      DO I = 1, N
         SUM = SUM + ( Y(I) - THETA(1)*( 1-EXP(-THETA(2)*X(I)) ) )**2
      END DO 
      OUT = OUT - 3*LOG(SUM) - 7.24085214731852D0 
      END
*
      SUBROUTINE CTLPST(X, OUT)
*
*     To compute posterior integrand for Contingency Table Data Problem
*              Example 2 from M. Evans and T. Swartz 
*     "Methods for Approxmating Integrals with Applications to Statistics"
*      with Special Emphasis on Bayesian Integration Problems"
*      Statistical Science, 10 (1995), 254-272. 
*
*     Original variables are tranformed using 
*     \theta = exp(x(1))/( 2 + 2*exp(x(1)) ),
*     \alpha_{1,1} = exp(x(2))/( 1 + exp(x(2)) + exp(x(3)) ),
*     \alpha_{2,1} = exp(x(3))/( 1 + exp(x(2)) + exp(x(3)) ),
*     \alpha_{1,2} = exp(x(4))/( 1 + exp(x(4)) + exp(x(5)) ),
*     \alpha_{2,2} = exp(x(5))/( 1 + exp(x(4)) + exp(x(5)) ),
*     \beta_{1,1} = exp(x(6))/( 1 + exp(x(6)) + exp(x(8)) ),
*     \beta_{2,1} = exp(x(7))/( 1 + exp(x(6)) + exp(x(8)) ),
*     \beta_{1,2} = exp(x(8))/( 1 + exp(x(8)) + exp(x(9)) ) and
*     \beta_{2,2} = exp(x(9))/( 1 + exp(x(8)) + exp(x(9)) ).
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(9)
*      DATA MODE/ 
*     &   1.40560,-2.14270,-0.47991, 1.63775,-0.33202,-1.73717,-0.03699,
*     &   2.63396, 1.66428/
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION X(9)
      INTEGER I, J, N(3,3)
      DATA N/ 43, 6, 9, 16, 11, 18, 3, 10, 16/
      DOUBLE PRECISION THETA, P, CSUM,  S, E1, E2, AL(3,2), BT(3,2)
      E1 = EXP( X(1) )
      THETA = E1/( 2*( 1 + E1 ) )
      CSUM = LOG( THETA/( 1 + E1 ) )
      DO J = 1,2
         E1 = EXP( X( 2*J     ) )
         E2 = EXP( X( 2*J + 1 ) )
         S = 1 + E1 + E2
         AL(1,J) = E1/S   
         AL(2,J) = E2/S
         AL(3,J) = 1/S
         E1 = EXP( X( 2*J + 4  ) )
         E2 = EXP( X( 2*J + 5 ) )
         S = 1 + E1 + E2
         BT(1,J) = E1/S   
         BT(2,J) = E2/S
         BT(3,J) = 1/S
      END DO
      DO I = 1,3
         DO J = 1,3
            P = THETA*AL(I,1)*BT(J,1) + ( 1 - THETA )*AL(I,2)*BT(J,2)
            CSUM = CSUM + N(I,J)*LOG(P)
         END DO
         CSUM  = CSUM + LOG( AL(I,1)*AL(I,2)*BT(I,1)*BT(I,2) )
      END DO
      OUT = CSUM 
      END
*
      SUBROUTINE EXLPST(X, OUT)
*
*     To compute log posterior for Extreme Value Distribution
*
*     Mode:
*      DOUBLE PRECISION MODE(M)
*      DATA MODE / M*0 /
*
*     Exact means should all be approximately -.57721
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION X(*), SUM
      INTEGER I, M
      COMMON /EXBLCK/M
      M = 5
      SUM = 0
      DO I = 1,M
         SUM = SUM + X(I) + EXP(-X(I))
      END DO
      OUT = - SUM
      END
*
      SUBROUTINE FLLPST(THETA, OUT)
*
*     To compute log posterior for Fuller Data
*     Reference: Fuller, W.A. (1976),
*       {\it Introduction to Statistical Time Series}, 
*       John Wiley and Sons, New York, p. 228.
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(3)
*      DATA MODE / 141.6062, -83.9727, 1.34778 /
*
      DOUBLE PRECISION OUT
      INTEGER I, N
      PARAMETER ( N = 12 )
      DOUBLE PRECISION THETA(3), SUM, Y(N), X(N)
      SAVE Y, X
      DATA ( Y(I), X(I), I = 1, N ) /
     &     47.3D0, 0D0,  87D0, 0.4D0, 120.1D0, 0.8D0, 130.4D0, 1.6D0, 
     &     58.8D0, 0D0, 111.9D0, 1D0, 136.5D0, 2D0, 132D0, 4D0, 68.8D0, 
     &     0D0, 138.1D0, 1.5D0, 145.7D0, 3D0, 143D0, 5.9D0 /
      SUM = 0
      DO I = 1, N
         SUM = SUM + ( Y(I)-THETA(1)-THETA(2)*EXP(-THETA(3)*X(I)) )**2
      END DO 
      OUT = -6.1D0*LOG( 10 + SUM/2 )
      END
*
      SUBROUTINE HRLPST(X, OUT)
*
*     To compute log posterior for Stanford Heart Transplant Data
*      Reference: 
*       Naylor, J. C. and Smith, A. F. M. (1988),  "Econometric 
*       Illustrations of Novel Numerical Integration Strategies for 
*       Bayesian Inference", J. Econometrics, 38, pp. 103-125.
*
*     Original variables are tranformed using \lambda = exp(x(1)),
*      \tau = exp(x(2)), and p = exp(x(3)). 
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(3)
*      DATA MODE/ 3.38498, -0.09179, -0.72319 /
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION X(3)
      DOUBLE PRECISION LAMBDA, TAU, P, LGL, LGT, LGP, FVAL, ESUM, ETSUM
      INTEGER I, NCASES
      PARAMETER ( NCASES = 82 )
      INTEGER WAIT(NCASES),SURVIV(NCASES),EXACT(NCASES),TREATD(NCASES)
      SAVE TREATD, WAIT, EXACT, SURVIV
      DATA TREATD/ 30*0, 52*1/
      DATA WAIT/
     *   49,    5,   17,    2,   39,   84,    7,    0,
     *   35,   36, 1400,    5,   34,   15,   11,    2,
     *    1,   39,    8,  101,    2,  148,    1,   68,
     *   31,    1,   20,  118,   91,  427,    0,   35,
     *   50,   11,   25,   16,   36,   27,   19,   17,
     *    7,   11,    2,   82,   24,   70,   15,   16,
     *   50,   22,   45,   18,    4,    1,   40,   57,
     *    0,    1,   20,   35,   82,   31,   40,    9,
     *   66,   20,   77,    2,   26,   32,   13,   56,
     *    2,    9,    4,   30,    3,   26,    4,   45,
     *   25,    5/
      DATA EXACT/ 10*1, 0, 16*1, 3*0, 14*1, 0, 3*1, 2*0, 
     *     4*1, 2*0, 1, 0, 1, 2*0, 3*1, 0, 1, 0, 1, 2*0,  
     *     3*1, 0, 1, 2*0, 2*1, 3*0/
      DATA SURVIV/ 30*0,
     *   15,    3,
     *  624,   46,  127,   61, 1350,  312,   24,   10,
     * 1024,   39,  730,  136, 1379,    1,  836,   60,
     * 1140, 1153,   54,   47,    0,   43,  971,  868,
     *   44,  780,   51,  710,  663,  253,  147,   51,
     *  479,  322,  442,   65,  419,  362,   64,  228,
     *   65,  264,   25,  193,  196,   63,   12,  103,
     *   60,   43/
      LGL = X(1)
      LAMBDA = EXP(LGL)
      LGT = X(2)
      TAU = EXP(LGT)
      LGP = X(3)
      P = EXP(LGP)
      FVAL = P*LGL*NCASES
      ESUM = 0
      ETSUM = 0
      DO I = 1, NCASES
         ESUM = ESUM + EXACT(I) 
         ETSUM = ETSUM + EXACT(I)*TREATD(I)
         FVAL = FVAL - (P+EXACT(I))*LOG( LAMBDA+WAIT(I)+TAU*SURVIV(I) )
      END DO
      OUT = X(1)+X(2)+X(3) + FVAL + ESUM*LGP + ETSUM*LGT
      END
*
      SUBROUTINE K3LPST(X, OUT)
      DOUBLE PRECISION X(3)
*
*     To compute log posterior for Econometric problem
*      Reference: 
*           Kloek, T., and van Dijk, H. K. (1978), 
*           ``Bayesian Estimates of Equation System Parameters: 
*           An  Application of Integration by Monte Carlo,''
*           Econometrica, 46, 1-19.
*
*     Original variables are tranformed using 
*      t_1 = exp(x(1))/( 1 + exp(x(1)) + exp(x(2)) ),
*      t_2 = exp(x(2))/( 1 + exp(x(1)) + exp(x(2)) ) and
*      t_3 = exp(x(3))/( 1 + exp(x(3)) ).
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(3)
*      DATA MODE/ 0.01221, -1.61596, -0.53693 /
*
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION T(3), PI2(2,2), OMEGN(3), ZPB1Z(3), BD(2,2)
      DOUBLE PRECISION S11, S12, S22, LTST, E
      SAVE OMEGN, ZPB1Z, PI2
*                   DON'T NEED INDIVIDUAL DATA
*                   OMEGN IS SSE MATRIX OF Y'S (C I) ON X'S (Z LAGI)
      DATA OMEGN/ .89754012354D0, .23746932393D0, .073946591778D0 / 
*                   ZPB1Z IS COV MX OF (Z LAGI)*DF OR CENTERED X'X
      DATA ZPB1Z/ 1.4176836D0, .1507912D0, .5799304D0 / 
*                   PI2 IS BOTTOM PART OF REDUCED FORM ESTIMATES
      DATA PI2/ .424716D0, 2.757635D0, .049935D0, 1.034737D0 /
*                   CHECK PARAMETER CONSTRAINTS 
         E = EXP( MAX( -4D1, MIN( 4D1, X(1) ) ) )
      T(2) = EXP( MAX( -4D1, MIN( 4D1, X(2) ) ) )
      T(1) =    E/( 1 + E + T(2) )
      T(2) = T(2)/( 1 + E + T(2) )
      T(3) = EXP( MAX( -4D1, MIN( 4D1, X(3) ) ) )
      T(3) = T(3)/( 1 + T(3) )      
*
*                   GAMMA MATRIX IS  ( T(1)-1    T(2)  )
*                                    (  T(1)    T(2)-1 )
*
*                   PUT N*GAMMA'*OMEGA*GAMMA INTO SS
*      
*
      S11 = ( T(1) - 1 )*( T(1) - 1 )*OMEGN(1) + 
     &        2*( T(1) - 1 )*T(1)*OMEGN(2) + T(1)*T(1)*OMEGN(3)
      S12 = ( T(1) - 1 )*T(2)*OMEGN(1) + ( T(1)*T(2) + 
     & ( T(1) - 1 )*( T(2) - 1 ) )*OMEGN(2) + T(1)*( T(2) - 1 )*OMEGN(3)
      S22 = ( T(2) - 1 )*( T(2) - 1 )*OMEGN(3) + 
     &        2*( T(2) - 1 )*T(2)*OMEGN(2) + T(2)*T(2)*OMEGN(1)
*
*                   NOW DIFFERENCE IN BETA FROM ESTIMATE
*                    GET BETA ESTIMATE FROM REDUCED FORM PI2 AND GAMMA
*
      BD(1,1) = T(1) - PI2(1,1) + T(1)*( PI2(1,1) + PI2(1,2) )
      BD(2,1) =  0   - PI2(2,1) + T(1)*( PI2(2,1) + PI2(2,2) )
      BD(1,2) = T(2) - PI2(1,2) + T(2)*( PI2(1,1) + PI2(1,2) )
      BD(2,2) = T(3) - PI2(2,2) + T(2)*( PI2(2,1) + PI2(2,2) )
*
*                   NOW GET SS -- ADD N*GAM*OMEG*GAM TO QF IN BD AND Z'Z
*
      S11 = S11 + BD(1,1)*BD(1,1)*ZPB1Z(1) + 2*BD(1,1)*BD(2,1)*ZPB1Z(2) 
     &          + BD(2,1)*BD(2,1)*ZPB1Z(3)
      S12 = S12 + BD(1,1)*BD(1,2)*ZPB1Z(1) + BD(1,1)*BD(2,2)*ZPB1Z(2) +
     &            BD(2,1)*BD(1,2)*ZPB1Z(2) + BD(2,2)*BD(2,1)*ZPB1Z(3)
      S22 = S22 + BD(1,2)*BD(1,2)*ZPB1Z(1) + 2*BD(1,2)*BD(2,2)*ZPB1Z(2)
     &          + BD(2,2)*BD(2,2)*ZPB1Z(3)
      OUT = -9*LOG( S11*S22 - S12*S12 )/2 
      LTST = ( 1 - T(1) - T(2) )**11*T(1)*T(2)*( 1 - T(3) )*T(3) 
      IF ( LTST .GT. 1D-40  ) THEN
         OUT = OUT + LOG( LTST )
      ELSE
         OUT = OUT - 100
      END IF
      END
*
      SUBROUTINE K5LPST(X, OUT)
      DOUBLE PRECISION X(5)
*
*     To compute log posterior for Econometric problem
*      Reference: 
*           Kloek, T., and van Dijk, H. K. (1978), 
*           ``Bayesian Estimates of Equation System Parameters: 
*           An  Application of Integration by Monte Carlo,''
*           Econometrica, 46, 1-19.
*
*     Original variables are tranformed using 
*      t_3 = exp(x(3))/( 1 + exp(x(3)) + exp(x(4)) ),
*      t_4 = exp(x(4))/( 1 + exp(x(3)) + exp(x(4)) ) and
*      t_5 = exp(x(5))/( 1 + exp(x(5)) ).
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(5)
*      DATA MODE/ 4.0859, -0.7206,  0.1162, -1.4866, -0.5329 /
*
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION  T(5), YCIZ(11,4), U1,U2, S11,S12,S22, LTST, E
      INTEGER I, J, IN
      SAVE IN, YCIZ
      DATA IN, ( ( YCIZ(I,J), J = 1,4 ), I = 1,11 ) / 0, 
     &     13895D0,10706D0,1024D0,2165D0,14377D0,10940D0,1078D0,2359D0,
     &     14843D0,11250D0,1123D0,2470D0,15307D0,11089D0,1052D0,3166D0,
     &     15360D0,11023D0, 980D0,3357D0,15951D0,11474D0,1073D0,3404D0,
     &     16680D0,12023D0,1281D0,3376D0,17237D0,12443D0,1474D0,3320D0,
     &     17547D0,12548D0,1591D0,3408D0,17788D0,12802D0,1668D0,3318D0,      
     &     17699D0,13096D0,1709D0,2894D0 /
      IF ( IN .EQ. 0 ) THEN
         DO J = 1, 4
            DO I = 1, 11
               YCIZ(I,J) = YCIZ(I,J)/1000
            END DO
         END DO
         IN = 5
      END IF
      T(1) = X(1)
      T(2) = X(2)
         E = EXP( MAX( -4D1, MIN( 4D1, X(3) ) ) )
      T(4) = EXP( MAX( -4D1, MIN( 4D1, X(4) ) ) )
      T(3) =    E/( 1 + E + T(4) )
      T(4) = T(4)/( 1 + E + T(4) )
      T(5) = EXP( MAX( -4D1, MIN( 4D1, X(5) ) ) )
      T(5) = T(5)/( 1 + T(5) )      
*      
      S11 = 0
      S12 = 0
      S22 = 0
*     NOTE INDEX STARTS AT 2 BECAUSE OF LAG INVESTMENT
      DO I = 2,11
*     FORM RESIDUALS
         U1 = YCIZ(I,2) - T(1) - T(3)*YCIZ(I,1)
         U2 = YCIZ(I,3) - T(2) - T(4)*YCIZ(I,1)
     &                         - T(5)*YCIZ(I-1,3)
*     NOW GET SS
         S11 = S11 + U1*U1
         S12 = S12 + U1*U2
         S22 = S22 + U2*U2
      END DO
      OUT = -5*LOG( S11*S22 - S12*S12 ) 
      LTST = ( 1 - T(3) - T(4) )**11*T(3)*T(4)*( 1 - T(5) )*T(5) 
      IF ( LTST .GT. 1D-40  ) THEN
         OUT = OUT + LOG( LTST )
      ELSE
         OUT = OUT - 100
      END IF
      END
*
      SUBROUTINE K8LPST(X, OUT)
      DOUBLE PRECISION X(*)
*
*     To compute log posterior for Econometric problem
*      Reference: 
*           Kloek, T., and van Dijk, H. K. (1978), 
*           ``Bayesian Estimates of Equation System Parameters: 
*           An  Application of Integration by Monte Carlo,''
*           Econometrica, 46, 1-19.
*
*     Original variables are tranformed using 
*      t_3 = exp(x(3))/( 1 + exp(x(3)) + exp(x(4)) ),
*      t_4 = exp(x(4))/( 1 + exp(x(3)) + exp(x(4)) ),
*      t_5 = exp(x(5))/( 1 + exp(x(5)) ),
*      t_6 = exp(x(6)) and
*      t_8 = exp(x(8)).
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(8)
*      DATA MODE/ 3.6524, -0.8612,  0.2569, -1.3141, -0.5332, -1.5009,
*     &  0.0709, -3.4271 /
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION  T(8), YCIZ(11,4), U1, U2, S11, S22, LTST, E
      INTEGER I, J, IN
      SAVE IN, YCIZ
      DATA IN, ( ( YCIZ(I,J), J = 1,4 ), I = 1,11 ) / 0, 
     &     13895D0,10706D0,1024D0,2165D0,14377D0,10940D0,1078D0,2359D0,
     &     14843D0,11250D0,1123D0,2470D0,15307D0,11089D0,1052D0,3166D0,
     &     15360D0,11023D0, 980D0,3357D0,15951D0,11474D0,1073D0,3404D0,
     &     16680D0,12023D0,1281D0,3376D0,17237D0,12443D0,1474D0,3320D0,
     &     17547D0,12548D0,1591D0,3408D0,17788D0,12802D0,1668D0,3318D0,      
     &     17699D0,13096D0,1709D0,2894D0 /
      IF ( IN .EQ. 0 ) THEN
         DO J = 1, 4
            DO I = 1, 11
               YCIZ(I,J) = YCIZ(I,J)/1000
            END DO
         END DO
         IN = 8
      END IF
      T(1) = X(1)
      T(2) = X(2)
         E = EXP( MAX( -4D1, MIN( 4D1, X(3) ) ) )
      T(4) = EXP( MAX( -4D1, MIN( 4D1, X(4) ) ) )
      T(3) =    E/( 1 + E + T(4) )
      T(4) = T(4)/( 1 + E + T(4) )
      T(5) = EXP( MAX( -4D1, MIN( 4D1, X(5) ) ) )
      T(5) = T(5)/( 1 + T(5) )      
      T(6) = EXP( MAX( -4D1, MIN( 4D1, X(6) ) ) )
      T(7) = X(7)
      T(8) = EXP( MAX( -4D1, MIN( 4D1, X(8) ) ) )
*      
      S11 = 0
      S22 = 0
*     NOTE INDEX STARTS AT 2 BECAUSE OF LAG INVESTMENT
      DO I = 2,11
*     FORM RESIDUALS
         U1 = YCIZ(I,2) - T(1) - T(3)*YCIZ(I,1)
         U2 = YCIZ(I,3) - T(2) - T(4)*YCIZ(I,1)
     &                         - T(5)*YCIZ(I-1,3)
         U1 = U1/T(6)
         U2 = ( U2 - T(7)*U1 )/T(8)
*     NOW GET SS
         S11 = S11 + U1*U1
         S22 = S22 + U2*U2
      END DO
      OUT = - ( S11 + S22 )/2 
      LTST = (1-T(3)-T(4))**11*T(3)*T(4)*(1-T(5))*T(5)/(T(6)*T(8))**12 
      IF ( LTST .GT. 1D-40  ) THEN
         OUT = OUT + LOG( LTST )
      ELSE
         OUT = OUT - 100
      END IF
      END
*
      SUBROUTINE LBLPST(THETA, OUT)
*
*     To compute log posterior for Lubricant Data problem taken from
*     Bates, D. M., and Watts, D. G. (1988),
*     Nonlinear Regression and Its Applications, Wiley, New York,
*     pp. 87-89, 275. 
*     
*     Original variable \sigma = exp(theta(10)) and original variables
*     \theta_1 and \theta_2 are scaled by 1/100 and 1/10, respectively.
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(10)
*      DATA MODE/ 
*     &      10.55D0, 20.66D0, 1.46D0, -0.259D0, 0.0225D0, 
*     &      0.402D0, .0352D0, 57.4D0, -0.476D0, -3.198D0 /
*
      DOUBLE PRECISION OUT
      INTEGER I, N, IS
      PARAMETER ( N = 53 )
      DOUBLE PRECISION THETA(*),SUM,DIF,SIGMA,W, T(N),P(N),V(N),TH1,TH2
      SAVE T, P, V
      DATA ( T(I), I = 1, N ), ( P(I), V(I), I = 1, N ) /
     &     11*0D0, 12*25D0, 15*37.8D0, 15*98.9D0, 
     &        1.000,  5.10595,  740.803,  6.38705, 1407.470,  7.38511, 
     &      363.166,  5.79057,    1.0  ,  5.10716,  805.500,  6.36113,
     &     1868.090,  7.97720, 3285.100, 10.47250, 3907.470, 11.92720,
     &     4125.470, 12.42620, 2572.030,  9.15630, 
     &        1.000,  4.54223,  805.500,  5.82452, 1505.920,  6.70515,
     &     2339.960,  7.71659,  422.941,  5.29782, 1168.370,  6.22654,
     &     2237.290,  7.57338, 4216.890, 10.35400, 5064.290, 11.98440, 
     &     5280.880, 12.44350, 3647.270,  9.52333, 2813.940,  8.34496, 
     &      516.822,  5.17275, 1737.990,  6.64963, 1008.730,  5.80754, 
     &     2749.240,  7.74101, 1375.820,  6.23206,  191.084,  4.66060,
     &        1.000,  4.29865, 2922.940,  7.96731, 4044.600,  9.34225, 
     &     4849.800, 10.51090, 5605.780, 11.82150, 6273.850, 13.06800,
     &     3636.720,  8.80445, 1948.960,  6.85530, 1298.470,  6.11898,
     &        1.000,  3.38099,  685.950,  4.45783, 1423.640,  5.20675, 
     &     2791.430,  6.29101, 4213.370,  7.32719, 2103.670,  5.76988,
     &      402.195,  4.08766,    1.000,  3.37417, 2219.700,  5.83919,
     &     3534.750,  6.72635, 4937.710,  7.76883, 6344.170,  8.91362,
     &     7469.350,  9.98334, 5640.940,  8.32329, 4107.890,  7.13210/
      SIGMA = EXP(THETA(10))
      SUM = 0
      TH1 = THETA(1)*100
      TH2 = THETA(2)*10
      DO I = 1, N
         W = P(I)/1000
         DIF = V(I) - TH1/( TH2 + T(I) )
     &       - W*( THETA(3) + W*( THETA(4) + W*THETA(5) ) )
         IF ( T(I) .LT. ABS( THETA(8)+THETA(9)*W**2 )*100 )
     &        DIF = DIF - W*( THETA(6) + W**2*THETA(7) )
     &                     *EXP( -T(I)/( THETA(8)+THETA(9)*W**2 ) )
         SUM = SUM + DIF**2
      END DO 
      OUT = -( N - 1 )*THETA(10) - SUM/(2*SIGMA**2) 
      END
*
*
      SUBROUTINE LGLPST(Z, OUT)
      INTEGER M, I
      DOUBLE PRECISION OUT
*
*     To compute log posterior for multivariate logistic distribution
*
      DOUBLE PRECISION Z(*), SUM, SUME
      COMMON /LGLBLK/M
      M = 5
      SUM = 0
      SUME = 0
      DO I = 1,M
         SUM = SUM + Z(I)
         SUME = SUME + LOG( 1 + EXP(-Z(I)) )
      END DO
      OUT = -SUM - 2*SUME
      END
*
      SUBROUTINE LMLPST(THETA, OUT)
*
*
*     To compute log posterior for Linear Model Data Problem from
*     Evans, M. and Swartz, T. (1995), 
*     "Methods for Apprixmating Integrals with Applications to Statistics
*      with Special Emphasis on Bayesian Integration Problems"
*      Statistical Science, 10, 254-272. 
*
*     Original variable \sigma is transformed using \sigma = exp(theta(10)).
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(10)
*      DATA MODE/
*     &   2.01819, 0.11640, 0.88300, 0.02861,-0.40753,-0.17198,-0.27059,
*     &  -0.62908, 0.36961,-0.23244/
*
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION THETA(*)
      DOUBLE PRECISION LAMBDA, SIGMA, Z, INRPRD, GAMCON
      INTEGER I, IC, J, K, M, N
      PARAMETER ( M = 10, N = 5 )
      PARAMETER ( LAMBDA = 3, GAMCON =  -0.451582705289499D0 )
*      GAMCON = LogGamma( ( LAMBDA + 1 )/2 ) - LogGamma( LAMBDA/2 )  
*     &          - Log( LAMBDA - 2 )/2 - Log(sqrt(\pi))
*      GAMCON = GAMMLN( ( LAMBDA + 1 )/2 ) - GAMMLN( LAMBDA/2 )  
*     &          - LOG( LAMBDA - 2 )/2 - 0.5723649429247000D0
      DOUBLE PRECISION Y(5,9)
      SAVE IC, Y
      DATA IC, Y / 0, 
     & 0.151754965,1.085632419, 0.591451482, -.068423422, -.378073973,
     & -.638321684,0.186636475, -.009514517, 0.495339037, 0.261643198,
     & 2.833250710,1.098441287, -.476050335, -.383252862, 1.128895153,
     & 0.041946195,0.174042986, -.581710481, 0.447193674, -.071116204,
     & -.616061349,-.513835377,-1.321292735, 0.043171116, 0.072504593,
     & -.816614843,0.934854744, -.156455228, -.172299695, -.174379640,
     & -.878534226,0.028068692, -.110434509,-1.193255325, 0.973952349,
     & -.610598648,-.548582027, -.544598812, -.402896237,-1.299835892,
     & 0.289558508,0.755147497, -.149990572, 0.746211189,-3.006485474 /
      IF ( IC .EQ. 0 ) THEN
         IC = 1
         DO J = 1,N
            Y(J,1) = Y(J,1) + 1.8D0
         END DO
      END IF
      SIGMA = EXP( THETA(M) )
      INRPRD = 0
      DO I = 1, M-1
         DO J = 1, N
            Z = (  Y(J,I) - THETA(I) )/SIGMA
            INRPRD = INRPRD + LOG( 1 + Z*Z/( LAMBDA - 2 ) )
         END DO
      END DO
      OUT = N*(M-1)*( GAMCON - THETA(M) ) - ( LAMBDA + 1 )*INRPRD/2 
      END
*
      SUBROUTINE MTLPST(BETA, OUT)
*
*     To compute log posterior for Motorette Data from
*      Tanner, M. A. (1993), Tools for Statistical Inference, 2nd Ed.,
*      Springer-Verlag, New York, p. 41.
*
*     Original variables are transformed using
*     \sigma  = exp( BETA(3) )
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(3)
*      DATA MODE/ -4.9306D0, 3.7471D0, -1.8631D0/
*
*
      DOUBLE PRECISION OUT
      INTEGER I, N, M, IC
      PARAMETER ( M = 17, N = 40 )
      DOUBLE PRECISION BETA(*), SUMSQR, SIGMA, T(N), V(N)
      INTEGER TEMP(N), TIME(N)
      SAVE TEMP, TIME, T, V, IC
      DATA IC, ( TEMP(I), I = 1, N ), ( TIME(I), I = 1, N ) / 0, 
     &     7*170, 5*190, 5*220, 10*150, 3*170, 5*190, 5*220, 
     &     1764, 2772, 3444, 3542, 3780, 4860, 5196, 2*408, 2*1344,
     &     1440, 2*408, 3*504, 10*8064, 3*5448, 5*1680, 5*528 /
      IF ( IC .EQ. 0 ) THEN
         DO I = 1,N
            V(I) = 1D3/( TEMP(I) + 273.2D0 )
            T(I) = LOG10( DBLE(TIME(I)) )
         END DO
         IC = 1
      END IF
      SIGMA = EXP(BETA(3))
      SUMSQR = 0
      DO I = 1, N
         SUMSQR = SUMSQR + ( T(I) - BETA(1) - BETA(2)*V(I) )**2
      END DO 
      OUT = -(N-1)*BETA(3) - SUMSQR/( 2*SIGMA**2 ) 
      END
*
      SUBROUTINE NLLPST(Z, OUT)
*
*     To compute log posterior for Nonlinear Regression Model from
*       P. M. Reilly, "The Numerical Computation of Posterior
*       Distributions iin Bayesian Statistical Inference",
*       Appl Statist. 25 (1976), pp. 201-209.
*
*     
*     Original variables are transformed using 
*     \alpha = exp(Z(1))
*     \beta  = exp(Z(2))
*     \sigma = exp(Z(3))
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(3)
*      DATA MODE / 1.3862D0, 0.8676D0, -5.969D0 /
*
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION Z(*)
      DOUBLE PRECISION ALPHA, BETA, SIGMA, FVAL
      INTEGER I, N
      PARAMETER ( N = 5 )
      INTEGER X(0:N)
      DOUBLE PRECISION  Y(0:N)
      SAVE X, Y
      DATA X/ 0, 1, 2, 3, 4, 5/
      DATA Y/ 4.11D0, 6.32D0, 8.21D0, 10.43D0, 14.29D0, 16.78D0 /
      ALPHA = EXP( Z(1) )
      BETA = EXP( Z(2) )
      SIGMA = EXP( Z(3) ) 
      FVAL = 0
      DO I = 0, N
         FVAL = FVAL - LOG( Y(I)/( ALPHA + BETA*X(I) ) )**2
      END DO
      OUT = Z(1) + Z(2) - 3*Z(3) + FVAL/(2*SIGMA)
      END
*
      SUBROUTINE NRLPST(X, OUT)
*
*     To compute log posterior for Normal Density
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION X(*), LGCON, SUMSQR
      INTEGER I, M
      PARAMETER ( M = 6, LGCON = .91893853320467274177D0 )
      SUMSQR = 0
      DO I = 1,M
         SUMSQR = SUMSQR + X(I)**2
      END DO
      OUT = - SUMSQR/2 - M*LGCON
      END
*
      SUBROUTINE PGLPST(THETA, OUT)
*
*     Log posterior for Porgi Example, taken from
*
*     Shaw, J. E. H. (1988), Aspects of Numerical Integration and 
*     Summarisation, in Bayesian Statistics 3, J. Bernardo, H. H. Degroot, 
*     D. V. Lindley and A. F. M. Smith (Eds.), Oxford University Press,   
*     pp. 411-428.
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(14)
*      DATA MODE/ 
*     & 0.0172,-0.2626,-0.3833, 1.4509, 0.1196, 0.0934, 1.4998, 
*     & 0.3177, 0.5976, 1.2672, 0.2683, 0.6650, 1.2261, 0.3795 /
*
      DOUBLE PRECISION OUT
      INTEGER I, J
      DOUBLE PRECISION THETA(14)
*
      DOUBLE PRECISION LOGSUM, LPRIOR, TWO, PTWO, HALF, NSUM
      PARAMETER ( TWO = 2, PTWO = TWO/10, HALF = 1/TWO )
      DOUBLE PRECISION MU(5), SIGMA(5), P(5), FF, FO, FI, PHID
      INTEGER FREQS(21)
      SAVE FREQS
      DATA ( FREQS(I),I = 1,21 )/
     &     509, 2240, 2341, 623, 476, 1230, 1439, 921, 448, 512, 
     &     719, 673, 445, 341, 310, 228, 168, 140, 114, 64, 22/
*
*     Computation of Prior 
*
      MU(1) = 11 + THETA(1)
      SIGMA(1) = EXP(THETA(2))
      P(1) =                    EXP(THETA(3 ))/(1+EXP(THETA(3)))
      P(2) = (1-P(1))          *EXP(THETA(6 ))/(1+EXP(THETA(6)))
      P(3) = (1-P(1)-P(2))     *EXP(THETA(9 ))/(1+EXP(THETA(9)))
      P(4) = (1-P(1)-P(2)-P(3))*EXP(THETA(12))/(1+EXP(THETA(12)))
      P(5) = (1-P(1)-P(2)-P(3)-P(4))
      LPRIOR = 0
      DO I = 1,4
         LPRIOR = LPRIOR - ( LOG( P(I+1)/P(I) ) + HALF )**2/2
         MU(I+1) = MU(I) + EXP(THETA(3*I+1))
         SIGMA(I+1) = EXP(THETA(3*I+2))
         LPRIOR = LPRIOR - (THETA(3*I+2) - THETA(3*I-1) + PTWO)**2/2
         LPRIOR = LPRIOR + THETA(3*I) + THETA(3*I+1)
     &          - ( 6 - I )*LOG( 1 + EXP(THETA(3*I) ) )
      END DO
      LPRIOR = LPRIOR - 4*LOG( MU(5) - MU(1) )
      DO I = 2,4
         LPRIOR = LPRIOR + THETA(3*I+1) + THETA(3*I-2)
     &          - 2*LOG( MU(I+1) - MU(I-1) )
      END DO
*
*     Compute LogLikelihood
*
      LOGSUM = LPRIOR
      NSUM = 0
      DO J = 1,22
         FF = 0
         DO I = 1,5
            FF = FF + P(I)*PHID( ( 8 + J - MU(I) )/SIGMA(I) )
         END DO
         IF ( J .GT. 1 ) THEN
            LOGSUM = LOGSUM + FREQS(J-1)*LOG( FF -FO )
            NSUM = NSUM + FREQS(J-1)
         ELSE
            FI = FF
         ENDIF
         FO = FF
      END DO
      OUT = LOGSUM - NSUM*LOG( FF - FI )
      END
      DOUBLE PRECISION FUNCTION PHID(Z)
*
*       Normal distribution probabilities accurate to 1.e-15.
*       Z = no. of standard deviations from the mean.
*
*       Based upon algorithm 5666 for the error function, from:
*       Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*
*       Programmer: Alan Miller
*
*       Latest revision - 30 March 1986
*
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, 
     &     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     &     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(  P0 = 220.20 68679 12376 1D0,
     &            P1 = 221.21 35961 69931 1D0, 
     &            P2 = 112.07 92914 97870 9D0,
     &            P3 = 33.912 86607 83830 0D0,
     &            P4 = 6.3739 62203 53165 0D0,
     &            P5 = .70038 30644 43688 1D0, 
     &            P6 = .035262 49659 98910 9D0 )
      PARAMETER(  Q0 = 440.41 37358 24752 2D0,
     &            Q1 = 793.82 65125 19948 4D0, 
     &            Q2 = 637.33 36333 78831 1D0,
     &            Q3 = 296.56 42487 79673 7D0, 
     &            Q4 = 86.780 73220 29460 8D0,
     &            Q5 = 16.064 17757 92069 5D0, 
     &            Q6 = 1.7556 67163 18264 2D0,
     &            Q7 = .088388 34764 83184 4D0 )
      PARAMETER(  ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER(  CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37
*     
      IF ( ZABS .GT. 37 ) THEN
         P = 0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP(-ZABS**2/2)
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL
     &           *( ( ( ( ( ( P6*ZABS + P5 )*ZABS + P4 )*ZABS + P3 )
     &                          *ZABS + P2 )*ZABS + P1 )*ZABS + P0 )
     &           /( ( ( ( ( ( ( Q7*ZABS + Q6 )*ZABS + Q5 )*ZABS + Q4 )
     &                            *ZABS + Q3 )*ZABS + Q2 )
     &                            *ZABS + Q1 )*ZABS + Q0 )
*
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/( ZABS + 1/( ZABS + 2/( ZABS + 3/( ZABS
     &                        + 4/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0 ) P = 1 - P
      PHID = P
      END
*
      SUBROUTINE PHLPST(X, OUT)
*
*     To compute log posterior for Photocarcinogen Data problem taken from
*     Dellaportas, P. amd Smith, A. F. M. (1993), Bayesian Inference for 
*     Generalized Linear and Proportional Hazards Models via Gibbs
*     Sampling, Appl. Statist. 42, 443-459.
*     
*     Original variable p is transformed p  = exp(Z(1))
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(5)
*      DATA MODE/ 1.1903, -10.8469, -1.1791, -0.3517, 0.4020 /
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION X(5)
      DOUBLE PRECISION BETA(0:3), P, SUM, ZS
      INTEGER I, J, N, M, NM, IC
      PARAMETER ( N = 65, M = 15, NM = 80 )
      INTEGER Z(3,NM), CEN(NM), WEEK(NM)
      DOUBLE PRECISION LNWEEK(NM)
      SAVE LNWEEK, Z
      DATA IC, WEEK/ 0, 
     &     12, 17, 21, 25, 11, 26, 27, 30, 13, 12,
     &     21, 20, 23, 25, 23, 29, 35, 40, 31, 36,
     &     32, 27, 23, 12, 18, 40, 40, 38, 29, 30,
     &     40, 32, 40, 40, 40, 40, 25, 30, 37, 27,
     &     22, 26, 10, 28, 19, 15, 12, 35, 35, 10,
     &     22, 18, 24, 12, 40, 40, 31, 24, 37, 29,
     &     27, 18, 22, 13, 18, 29, 28, 20, 16, 22, 
     &     26, 19, 29, 10, 17, 28, 26, 12, 17, 26/
      DATA ((Z(J,I),I=1,NM),J=1,3)/20*0, 20*1, 80*0, 20*1, 80*0, 20*1/
      DATA CEN/ 
     &     17*0, 1, 7*0, 2*1, 3*0, 1, 0, 4*1, 6*0,
     &     1, 9*0, 1, 0, 2*1, 11*0, 1, 4*0, 2*1, 6*0/
      IF ( IC .EQ. 0 ) THEN
         IC = 1
         DO I = 1,NM
            LNWEEK(I) = LOG( DBLE(WEEK(I)) )
         END DO
      END IF
      P = EXP(X(1))
      BETA(0) = X(2) 
      BETA(1) = X(3)
      BETA(2) = X(4)
      BETA(3) = X(5)
      SUM = (N+1)*X(1)
      DO I = 1,NM
         ZS = BETA(0) + BETA(1)*Z(1,I)+ BETA(2)*Z(2,I)+ BETA(3)*Z(3,I)
         IF ( CEN(I) .EQ. 0 ) SUM = SUM + ZS + (P-1)*LNWEEK(I) 
         SUM = SUM - EXP( ZS + P*LNWEEK(I) )
      END DO
      OUT = SUM
      END

*
      SUBROUTINE PRLPST(X, OUT)
      DOUBLE PRECISION X(7)
*
*     To compute log posterior for Proportional Hazards Problem
*      Reference: 
*           Dellaportas, P., and Wright, D. (1992), 
*           A Numerical Integration Strategy in Bayesian Analysis,
*           in Bayesian Statistics 4, J. M. Bernardo, et al. eds, 
*           Oxford University Press, 601-606.
*
*     Original variables p is tranformed using p = exp(x(1)). 
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(7)
*      DATA MODE/ 
*      &   0.13732,-3.90576, 1.87334,-0.13476,-0.02137,-0.03869, 0.13793/
*
      DOUBLE PRECISION OUT
      INTEGER I, J, N, M
      PARAMETER ( N = 48, M = 17 )
      DOUBLE PRECISION P, BETA(0:5), SUM, ZS, LNT, Z(5,M+N), T(M+N)
      SAVE T, Z
      DATA ( T(I), (Z(J,I), J = 1,5 ), I = 1, 20 ) /
     &  1.0,  0.8262308, -0.8015385,  6.8461538, 0.0,-0.1230769,
     &  1.0,  0.5482308,  1.7984615,-22.1538462, 0.0, 7.8769231,
     &  2.0,  0.1272308, -0.4015385, 20.8461538, 0.0, 4.8769231,
     &  2.0,  0.3562308,  1.0984615, 14.8461538, 0.0, 1.8769231,
     &  2.0, -0.0907692, -5.1015385, -3.1538462, 0.0,-1.1230769,
     &  3.0,  0.1522308, -3.5015385,-14.1538462, 1.0,-0.1230769,
     &  5.0,  0.8442308, -0.1015385,-10.1538462, 1.0,-1.1230769,
     &  5.0,  0.2892308, -3.7015385, 13.8461538, 0.0,-1.1230769,
     &  6.0, -0.0297692, -1.2015385, 16.8461538, 0.0,-2.1230769,
     &  6.0,  0.7222308, -0.0015385,  9.8461538, 1.0,-2.1230769,
     &  6.0, -0.2777692, -0.5015385, -0.1538462, 0.0,-0.1230769,
     &  6.0,  0.0232308,  0.1984615,  6.8461538, 1.0,-2.1230769,
     &  7.0,  0.5862308, -0.7015385,-12.1538462, 0.0,-0.1230769,
     &  7.0, -0.3507692, -5.1015385,  0.8461538, 1.0,-0.1230769,
     &  7.0, -0.2157692,  1.1984615, -7.1538462, 1.0, 2.8769231,
     &  9.0,  0.3322308, -2.0015385, -5.1538462, 0.0, 1.8769231,
     & 11.0, -0.2777692,  3.7984615,  0.8461538, 0.0,-0.1230769,
     & 11.0, -0.1617692,  1.7984615,-17.1538462, 0.0,-1.1230769,
     & 11.0, -0.0907692,  2.9984615,  4.8461538, 0.0,-0.1230769,
     & 11.0,  0.1162308, -2.7015385,  9.8461538, 0.0, 1.8769231/
      DATA ( T(I), (Z(J,I), J = 1,5 ), I = 21, 40 ) /
     & 11.0, -0.3127692, -0.6015385, -9.1538462, 1.0,-1.1230769,
     & 13.0, -0.6137692, -4.7015385, -0.1538462, 1.0,-0.1230769,
     & 14.0,  0.0062308,  4.3984615,  5.8461538, 0.0,-0.1230769,
     & 15.0,  0.2102308,  0.3984615,  9.8461538, 0.0, 0.8769231,
     & 16.0, -0.0497692, -1.2015385,-12.1538462, 0.0,-0.1230769,
     & 16.0, -0.0697692, -1.4015385,  1.8461538, 1.0,-0.1230769,
     & 17.0, -0.1617692, -0.2015385, -7.1538462, 0.0,-1.1230769,
     & 17.0,  0.1992308,  0.9984615,  7.8461538, 0.0,-0.1230769,
     & 18.0,  0.0552308, -2.7015385,  4.8461538, 1.0,-2.1230769,
     & 19.0, -0.3127692,  4.1984615, -9.1538462, 0.0, 4.8769231,
     & 19.0, -0.1367692, -2.7015385, -0.1538462, 1.0,-1.1230769,
     & 24.0, -0.0907692,  4.3984615, -4.1538462, 1.0,-1.1230769,
     & 25.0, -0.3917692,  2.1984615,  6.8461538, 0.0,-0.1230769,
     & 26.0, -0.1617692,  0.9984615,-11.1538462, 1.0, 0.8769231,
     & 32.0, -0.0697692,  0.3984615,-14.1538462, 0.0,-1.1230769,
     & 35.0, -0.2777692, -3.2015385,-12.1538462, 0.0,-0.1230769,
     & 37.0,  0.2102308,  0.7984615,  2.8461538, 0.0,-1.1230769,
     & 41.0, -0.3917692, -0.0015385,  8.8461538, 0.0,-0.1230769,
     & 42.0, -0.2457692, -5.2015385,  9.8461538, 1.0,-1.1230769,
     & 51.0,  0.1762308, -2.5015385, 13.8461538, 0.0, 2.8769231/
      DATA ( T(I), (Z(J,I), J = 1,5 ), I = 41, 60 ) /
     & 52.0, -0.3917692, -0.1015385, -0.1538462, 1.0,-0.1230769,
     & 54.0, -0.1367692, -1.2015385,-11.1538462, 0.0,-0.1230769,
     & 58.0, -0.1877692,  1.8984615,-18.1538462, 1.0,-0.1230769,
     & 66.0,  0.0552308, -3.6015385, -1.1538462, 0.0,-1.1230769,
     & 67.0, -0.0697692,  2.5984615, -8.1538462, 0.0,-0.1230769,
     & 88.0, -0.2157692,  0.3984615,-13.1538462, 1.0,-1.1230769,
     & 89.0, -0.0697692,  3.7984615,  2.8461538, 0.0,-1.1230769,
     & 92.0,  0.0392308,  0.7984615, -2.1538462, 1.0, 0.8769231,
     &  4.0,  0.5532308, -0.0015385, -1.1538462, 0.0,-0.1230769,
     &  4.0,  0.5322308, -0.2015385,-11.1538462, 1.0, 2.8769231,
     &  7.0, -0.2777692,  2.1984615,-12.1538462, 1.0,-0.1230769,
     &  7.0,  0.1402308, -0.0015385, 20.8461538, 0.0, 0.8769231,
     &  8.0, -0.3127692, -0.3015385, -3.1538462, 1.0,-2.1230769,
     & 12.0, -0.2457692,  1.3984615,-14.1538462, 1.0,-3.1230769,
     & 11.0,  0.2212308,  3.7984615, -0.1538462, 0.0,-1.1230769,
     & 12.0,  0.0062308, -1.4015385,  5.8461538, 1.0,-1.1230769,
     & 13.0,  0.2712308, -5.3015385, 10.8461538, 1.0,-1.1230769,
     & 16.0, -0.2457692,  2.7984615, -5.1538462, 0.0,-1.1230769,
     & 19.0, -0.0697692,  2.7984615, -1.1538462, 1.0,-0.1230769,
     & 19.0, -0.0697692,  0.5984615,  8.8461538, 1.0,-0.1230769/
      DATA ( T(I), (Z(J,I), J = 1,5 ), I = 61, N+M ) /
     & 28.0, -0.1617692, -2.9015385, 21.8461538, 1.0,-1.1230769,
     & 41.0,  0.3642308,  2.5984615, 11.8461538, 0.0,-1.1230769,
     & 53.0, -0.2777692,  1.7984615,  5.8461538, 0.0, 0.8769231,
     & 57.0, -0.1367692,  2.2984615,  5.8461538, 0.0, 0.8769231,
     & 77.0, -0.3127692,  3.7984615, -0.1538462, 0.0, 1.8769231/
      P = EXP( X(1) )
      BETA(0) = X(2)
      BETA(1) = X(3)
      BETA(2) = X(4)
      BETA(3) = X(5)
      BETA(4) = X(6)
      BETA(5) = X(7)
      SUM = 0 
      DO I = 1, N+M
         ZS = BETA(0) + BETA(1)*Z(1,I) + BETA(2)*Z(2,I) + BETA(3)*Z(3,I)
     &                + BETA(4)*Z(4,I) + BETA(5)*Z(5,I)
         LNT = LOG( T(I) )
         IF ( I .LE. N ) SUM = SUM + ZS + (P-1)*LNT
         IF ( ZS + P*LNT .GT. -100 ) SUM = SUM - EXP( ZS + P*LNT )
      END DO
      OUT = SUM + (N+1)*X(1)
      END
*
      SUBROUTINE PSLPST(X, OUT)
*
*     To compute log posterior for Pearson Type IV, example from
*     Genz, A., and Kass, R. (1997),
*     Subregion Adaptive Integration of Functions Having a Dominant Peak,
*     to appear in J. Comp. Graph. Stat.
*
*     Mode at X = 32, mean at X = 160/3.
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION PITWO, X, T, LAMBDA, OMEGA, RHO, NU, SQNU
      PARAMETER ( PITWO = 1.57079 63267 94897 ) 
      PARAMETER ( LAMBDA = 0, OMEGA = 1, RHO = 20, NU = 4, SQNU = 2 ) 
      T = ( X - LAMBDA )/( OMEGA*SQNU )
      OUT = RHO*NU*( ATAN(T) - PITWO ) - ( NU + 1 )*LOG( 1 + T**2 )/2
      END
*
      SUBROUTINE RDLPST(Z, OUT)
*
*     To compute log posterior for Radiotherapy Data example from
*      Tanner, M. A. (1993), Tools for Statistical Inference, 2nd Ed.,
*      Springer-Verlag, New York, p. 14.
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(2)
*      DATA MODE/ 3.81944D0, -0.08648D0 /
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION ALPHA, BETA, P, FVAL, Z(2)
      INTEGER I, N
      PARAMETER ( N = 24 )
      INTEGER DAYS(N), RESPNS(N)
      SAVE DAYS, RESPNS
      DATA DAYS/ 21,24,25,26,28,31,33,34,35,37,43,49,
     &           51,55,25,29,43,44,46,46,51,55,56,58 /,
     &     RESPNS/ 14*1,10*0 /
      ALPHA = Z(1)
      BETA = Z(2)
      FVAL = 0
      DO 10 I = 1,N
         P = ALPHA + BETA*DAYS(I)
         IF ( P .GT. 100 ) THEN
            FVAL = FVAL + ( RESPNS(I) - 1 )*P
         ELSE IF ( P .LT. -100 ) THEN
            FVAL = FVAL + RESPNS(I)*P
         ELSE 
            FVAL = FVAL + RESPNS(I)*P - LOG( 1+EXP(P) )
         END IF
 10   CONTINUE
      OUT = FVAL
      END

*
      SUBROUTINE TNLPST( ALPHBT, OUT )
*
*     To compute log posterior for Tornado Data Problem
*     Monahan, J. F., Schrab, K. J., and Anderson, C. E. (1993), 
*     ``Statistical Methods for Forecasting Tornado Intensity,'' 
*      in Statistical Sciences and Data Analysis, 
*      K. Matushita, et al., eds., VSP, Utrecht, The Netherlands, 13-24.
*
*     Approximate mode:
*      DOUBLE PRECISION MODE(11)
*      DATA MODE/
*     &    0.6584,  0.2414,  0.1538, -0.8338,  0.1092, 17.8495, 18.3402,
*     &   20.0268, 21.2174, 23.2027, 28.2810 /
*
      DOUBLE PRECISION OUT
      DOUBLE PRECISION ALPHBT(11) 
      DOUBLE PRECISION LGSUM, LGFACT, FLOGIT, BETAX
      INTEGER ND, I, J, KI
      PARAMETER ( ND = 157 )
      INTEGER UMAX(ND), MDA(ND), INTNS(ND), VORT(ND) 
      DOUBLE PRECISION SRAW(ND)
      SAVE UMAX, MDA, INTNS, VORT, SRAW
      DATA (  UMAX(I), I = 1, ND )/
     &  10, 15, 20, 10, 13, 18,  8,  5,  4, 18, 12,  5,  8, 13, 10, 15,
     &  15, 12, 10, 13, 15, 16, 33, 30, 27, 20, 15, 16, 30,  0,  8, 12,
     &  22, 25, 27, 50, 12, 12, 12, 12, 15, 24, 20,  8, 15,  7,  5,  5,
     &  10, 20, 30, 10,  2,  2,  4,  8, 20, 15, 38,  2, 10,  3, 25, 25,
     &  26, 31, 20, 10, 20, 20,  5,  6, 10, 15, 19, 10,  2, 15, 22,  8,
     &  10, 13,  9, 13,  5, 12, 20, 10, 10, 10, 18,  7,  7, 15, 23,  8,
     &   5, 20, 15,  8,  8, 15, 17, 13, 15, 10, 16, 20, 15, 13,  2, 15,
     &  10,  8, 30, 15, 18, 13, 17, 13, 10, 15, 13, 10, 15,  5,  0,  5,
     &   0,  5,  6, 10,  0, 10, 10,  3, 10, 10,  3, 12, 15, 15,  8, 13,
     &   5, 15,  8, 10, 15,  7,  5, 18, 13, 10, 10, 15, 23/
      DATA (   MDA(I), I = 1, ND )/
     &  -4, 15, 15, 11, 14, 11, -5,  5, -4, 13, 17,  7,  7,  0,  5, 21,
     &  16, 16, 15, 10, 30, 22, 32,  7, 35, 16, 22, 34, 51, 17, 10, 21,
     &  11, 13, 15, 15,  8,  0,  5,  8, -1, 13, 12,  5, 21, -9, -3, 13,
     & -13, 12, 18, -3, -2, 10, -1,  2,  9, 17, 38,  9, -5,  2, 36, 54,
     &  29, 16, 25, 14, 29, 28,  9, 15, -7, 18, 15, 22, 18, 25, 11, -8,
     &   5,  0,  6, 17, 12,  8, 23, 24, -6, -8, 23, 13, -3, -6, 15, 18,
     &  13, -1, -3,  2, -7, -8, 10,  1, 24, -4, 11, 16, 17, 10, 18, 11,
     &  18, 21, 21,  4, 14,  9, 13,  1,  7, 14,  7, 42, 13, 11, -1, -5,
     &  13, -6, 15, -4,  2,  4,  3,  6,  8,  8, 11, 13,  7, -4,  2, -6,
     &  11, 12, 28, 15, 14, -3, -1, 19, 18, 16, 14, 15,  9/
      DATA ( INTNS(I), I = 1, ND )/
     &   0,  0,  3,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2,
     &   0,  0,  0,  0,  1,  0,  4,  0,  5,  0,  0,  2,  5,  0,  0,  0,
     &   3,  4,  5,  4,  0,  0,  0,  0,  0,  3,  3,  0,  4,  0,  0,  0,
     &   0,  4,  5,  0,  0,  0,  0,  0,  2,  0,  5,  0,  0,  0,  5,  6,
     &   4,  4,  2,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   0,  0,  0,  0,  0,  0,  3,  3,  0,  0,  5,  0,  0,  0,  2,  0,
     &   0,  0,  0,  0,  0,  0,  4,  0,  3,  0,  3,  4,  2,  0,  0,  0,
     &   0,  0,  4,  0,  2,  0,  0,  0,  0,  0,  0,  5,  2,  0,  0,  0,
     &   0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &   0,  0,  0,  1,  5,  0,  0,  4,  2,  0,  1,  1,  3/
      DATA (  VORT(I), I = 1, ND )/
     &  18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
     &  19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
     &  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 54, 54, 54, 54,
     &  54, 54, 54, 54, 54, 54, 54, 54, 54,  9,  9,  9,  9,  9,  9,  9,
     &   9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 53,
     &  53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 48, 48, 48, 48, 48,
     &  48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 24, 24, 24, 24, 24,
     &  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 26, 26, 26, 26, 26,
     &  26, 26, 26, 26, 26, 26, 26, 26, 45, 45, 45, 45, 45, 45, 45, 45,
     &  45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45/
      DATA (  SRAW(I), I = 1, ND )/
     &  19.2, 11.1, 19.2, 14.8, 18.4, 15.1, 23.7, 16.2, 22.6, 12.3,
     &  12.9, 22.4, 18.6, 15.0, 17.7,  3.7,  9.0, 11.0, 15.0, 16.0,
     &  13.0, 11.0, 15.0,  8.0,  7.0, 12.0, 10.0,  8.0,  9.0, 14.0,
     &  17.0,  8.0, 20.0, 21.0, 21.0, 19.0, 18.0, 21.0, 14.0, 17.0,
     &  16.0, 28.0, 31.0,  3.0, 15.7, 10.5,  5.7,  9.7,  8.0, 13.8,
     &  16.1, 11.6, 10.5, 11.5, 11.3, 10.7,  7.2, 26.0, 23.0, 15.0,
     &  11.0,  9.0, 17.0, 13.0, 14.0, 21.0, 22.0, 18.0, 21.0, 22.0,
     &  17.0, 10.0,  9.0, 13.0, 14.0,  4.0, 14.0,  8.0, 20.0, 16.2,
     &   1.9, 15.1, 10.5, 12.7, 20.1, 24.6, 14.1, 14.7, 18.0, 16.9,
     &  25.2,  1.3,  2.6, 14.9, 24.7,  6.1,  3.9,  8.7, 21.3, 22.6,
     &  10.3, 12.6, 16.8, 15.6, 18.3, 16.4, 18.6, 21.4, 20.8, 15.2,
     &  14.8, 20.2, 20.1, 16.2, 16.5, 20.1, 23.8, 11.8, 13.0,  6.8,
     &  17.7, 10.1, 22.5, 23.9, 33.3, 14.6, 23.8, 21.5, 22.6, 14.9,
     &  20.6, 25.5, 18.2, 15.4, 19.5, 19.7, 10.2, 14.8, 10.9,  9.8,
     &  19.7, 15.9, 10.0, 11.6, 16.0,  7.7, 16.2, 16.5, 16.7,  5.9,
     &  20.5, 10.6, 18.8, 14.4, 22.5, 14.9, 17.0/
*
*     Compute Product for Likelihood
*
      LGSUM = 0
      DO I = 1,ND
         BETAX = UMAX(I)*ALPHBT(1)+MDA(I)*ALPHBT(2)+SRAW(I)*ALPHBT(3)      
     *         + UMAX(I)**2*ALPHBT(4)/100 + VORT(I)*ALPHBT(5)
         KI = INTNS(I)
         IF ( KI .EQ. 0 ) THEN
            LGFACT = FLOGIT( ALPHBT(6) - BETAX )
         ELSE IF ( KI .EQ. 6 ) THEN
            LGFACT = FLOGIT( BETAX - ALPHBT(11) )
         ELSE
            LGFACT = FLOGIT( ALPHBT(KI+6) - BETAX ) 
     &             - FLOGIT( ALPHBT(KI+5) - BETAX )
         ENDIF
         LGFACT = MAX( 1D-15, LGFACT )
         LGSUM = LGSUM + LOG(LGFACT)
      END DO
      OUT = LGSUM
      END
      DOUBLE PRECISION FUNCTION FLOGIT(S)
      DOUBLE PRECISION S, E
      IF ( S .LT. -4D1 ) THEN 
         FLOGIT = 2D-17
      ELSE IF ( S .GT. 4D1 ) THEN
         FLOGIT = 1
      ELSE
         E = EXP( S )
         FLOGIT = E/( 1 + E )
      END IF
      END
*
