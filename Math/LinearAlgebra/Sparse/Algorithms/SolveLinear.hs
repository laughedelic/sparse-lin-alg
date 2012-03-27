module Math.LinearAlgebra.Sparse.Algorithms.SolveLinear
where

import Data.Monoid
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Matrix 
import Math.LinearAlgebra.Sparse.Vector
import Math.LinearAlgebra.Sparse.Algorithms.Staircase 
import Math.LinearAlgebra.Sparse.Algorithms.Diagonal

solveLinSystems :: Integral α => SparseMatrix α -> SparseMatrix α -> SparseMatrix α
solveLinSystems m bs = 
    let as = trans bs
        (d,t,u) = toDiag m
        solve = fmap vec . solveDiagonal (d,t,u) . SV (height d)
    in SM (height as, width m) $ M.mapMaybe solve (mx as)

solveLinear :: Integral α =>SparseMatrix α -> SparseVector α -> Maybe (SparseVector α)
solveLinear m b = solveDiagonal (toDiag m) b 

solveDiagonal :: Integral α =>(SparseMatrix α, SparseMatrix α, SparseMatrix α)-> SparseVector α -> Maybe (SparseVector α)
solveDiagonal (d,t,u) b
    | height d /= dim b = error "width m /= dim b"
    | isZeroVec b = Just $ zeroVec (width d)
    | isZeroMx  d = Nothing
    | otherwise   =
    let (dd,a)    = (mainDiag d, t ×· b)
        (bad,sol) = solveDiag dd a
    in if dim dd < dim a then Nothing
          else if not (M.null bad) then Nothing
                  else Just $ (SV (width d) sol) ·× u

solveDiag :: Integral c =>SparseVector c -> SparseVector c -> (IntMap Index, IntMap c)
solveDiag d a = M.mapEitherWithKey solveOne (vec d)
    where solveOne i 0 = Left i
          solveOne i _ | a!i == 0 = Right 0
          solveOne i x = if (a!i) `mod` x /= 0 
                            then Left i 
                            else Right ((a!i) `div` x)
