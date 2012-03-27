module Math.LinearAlgebra.Sparse.Algorithms.Diagonal
where

import Data.Monoid
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Matrix 
import Math.LinearAlgebra.Sparse.Vector
import Math.LinearAlgebra.Sparse.Algorithms.Staircase 

isDiag :: SparseMatrix α -> Bool
isDiag m = M.foldlWithKey' (\true i row -> true && (keys row == [i])) True (mx m)

toDiag m = toDiag_ m (idMx (height m))
                     (idMx (width  m))

toDiag_ m t u | isZeroMx m =  (m,t,u)
              | height m /= height t = error "height(mM) /= height(mT)"
              | width m /= width u   = error "width(mM) /= width(mU)"
              | otherwise = 
              let (s,t') = staircase' m t
              in  dm True s t' u

  -- Here  s  is a staircase matrix.
  -- If it is not diagonal, then transp(s) it brought to staircase
  -- form - this corresponds to the column elementary 
  -- transformations of s -  and so on,  until the diagonal matrix 
  -- is obtained (even number of `transp' to be applied).

dm even s t u = case (even, isDiag s, trans s) of
    (True , True , _ ) -> (s , t, u)
    (False, True , s') -> (s', t, u)
    (True , False, s') -> dm False s'' t  u' where (s'',u') = staircase' s' u
    (False, False, s') -> dm True  s'' t' u  where (s'',t') = staircase' s' t
