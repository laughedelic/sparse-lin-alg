module Math.LinearAlgebra.Sparse.Algorithms.SolveLinear
(
solveLinear, solveLinSystems
)
where

import Data.Maybe
import Data.Monoid
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Matrix
import Math.LinearAlgebra.Sparse.Vector
import Math.LinearAlgebra.Sparse.Algorithms.Staircase
import Math.LinearAlgebra.Sparse.Algorithms.Diagonal

-- | Solves system for matrix in diagonal form
solveDiagonal :: Integral α
              => (SparseMatrix α, SparseMatrix α, SparseMatrix α)   -- ^ return triple-value of `toDiag`
              -> SparseVector α                                     -- ^ right-hand-side vector
              -> Either String (SparseVector α)                     -- ^ either error message or solution
solveDiagonal (d,t,u) b
    | height d /= dim b = Left "width m /= dim b"
    | isZeroVec b = Right $ zeroVec (width d)
    | isZeroMx  d = Left "left side is zero-matrix, while right side is not zero-vector"
    | otherwise   =
    let (dd,a)    = (mainDiag d, t ×· b)
        (bad,sol) = solveDiag dd a
    in if size (vec dd) < size (vec a) then Left   "zeroes in lhs, while non-zeroes in rhs"
          else if not (M.null bad)     then Left $ "integral fraction error at lines "++(show (elems bad))
                  else Right $ (SV (width d) sol) ·× u

solveDiag :: Integral α => SparseVector α -> SparseVector α -> (IntMap Index, IntMap α)
solveDiag dd a = M.mapEitherWithKey solveOne (vec a)
    where solveOne i r = let l = dd!i
                             (x,e) = r `divMod` l
                         in if (l == 0 && r /= 0) || e /= 0
                               then Left i else Right x

-- | Just solves system of linear equations in matrix form
--   for given left-side matrix and right-side vector
solveLinear :: Integral α =>SparseMatrix α -> SparseVector α -> Maybe (SparseVector α)
solveLinear m b = case solveDiagonal (toDiag m) b of
                       Left msg -> error msg
                       Right s -> Just s

--  Solves a system in form:
--
-- >>> 4×3       3×5         4×5
-- >>> [     ]   [       ]   [       ]
-- >>> [  m  ] × [   x   ] = [   b   ]
-- >>> [     ]   [       ]   [       ]
-- >>> [     ]               [       ]
--
-- @solveLinSystems m b@ returns solution matrix @x@
--

-- | Solves a set of systems for given left-side matrix and each right-side vector of given set (sparse vector)
solveLinSystems :: Integral α
                => SparseMatrix α                        -- ^ left-side matrix
                -> SparseVector (SparseVector α)         -- ^ right-side vectors
                -> Maybe (SparseVector (SparseVector α)) -- ^ if one of systems has no solution `Nothing` is returned
solveLinSystems m bs = if ok then Just (SV (dim bs) sols) else Nothing
    where (ok,sols)  = M.mapAccum solve True (vec bs)
          solve ok b = case solveLinear m b of
                            Nothing -> (False, emptyVec)
                            Just s  -> (ok, s)
--solveLinSystems m bs = if Nothing `elem` sols
--                          then error "system is not solvable"
--                          else catMaybes sols
--    where sols = fmap (solveLinear m) bs
