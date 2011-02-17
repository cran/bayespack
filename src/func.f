      DOUBLE PRECISION FUNCTION USRLGP (X)
      DOUBLE PRECISION X(*), usrlgpC
      USRLGP = usrlgpc(X)
      END
    
      SUBROUTINE USRMNS (X, GFUNS )
      DOUBLE PRECISION X(*), GFUNS(*)
      CALL usrmnsC(X, GFUNS)
      END
