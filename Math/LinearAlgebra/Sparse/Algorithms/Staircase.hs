module Math.LinearAlgebra.Sparse.Algorithms.Staircase 
    (staircase, clearColumn, extGCD)
where

import Math.LinearAlgebra.Sparse.Matrix 
import Math.LinearAlgebra.Sparse.Vector

-- | Staircase Form of matrix.
--   Using of `divMod` causes `Integral` context. (TODO: eliminate it)
--   Method:
--   Gauss method applied to the rows of matrix. Though α may be not
--   a field, we repeat the remainder division to obtain zeroes down
--   in the column. 
staircase :: (Num α, Integral α) => SparseMatrix α -> SparseMatrix α
staircase m | height m <= 1 = m  -- m is either zero matrix or one-row
            | otherwise    =     -- Main loop. m is non-zero.
    let m'   = clearColumn m
        row1 = m' `row` 1
    in if dim row1 == 1 then m'  -- we reached last column
          else if row1!1 == 0
                  -- this means that the first column is zero —+
                  then addZeroCol 1  -- × place it back        |
                     $ staircase     -- ↑ apply recursion      |
                     $ delCol 1 m'   -- ↑ cut this column <————'
                  -- else m'(1,1) == gcd(column) /= 0
                  -- and  m'(i,1) == 0 for i>1        —————————————————————+
                  else addRow row1 1     -- × and first row back           |
                     $ addZeroCol 1      -- ↑ return zero column           |
                     $ staircase         -- ↑ apply recursion              |
                     $ delRowCol 1 1 m'  -- ↑ so, we take smaller matrix <—'

-- | clearColumn m --> m'  
--   From the start, length(m) > 1.
--   m'(1,1) = gcd(firstColumn(m)), m'(i,1)==0  for i>1.
--   m'(1,1) = 0 means that column was zero.
clearColumn :: (Num α, Integral α) => SparseMatrix α -> SparseMatrix α
clearColumn m = c nz rest
    where (nz, rest) = partitionMx (\r -> r!1 /= 0) m

          -- Each ai = nz # (i,1) is non-zero,  
          -- The subcolumn (a1,a2) reduces to the form (a,0) by
          -- the Euclidean gcd algorithm, and the transformation  
          -- 2x2 matrix tt is accumulated, then it is applied to
          -- nz' without 1st column and nz(2) moves to
          -- rest. This continues while (heigth nz) > 1.
          c nz rest             -- nz are the rows with the non-zero head
              | height nz == 0 = rest                        -- zero column
              | height nz == 1 = addRow (nz `row` 1) 1 rest  -- single non-zero
              | otherwise     =
                let (r1,r2)   = (nz `row` 1, nz `row` 2)
                    nz'       = delRow 1 $ delRow 1 $ nz
                    (a,tt)    = extGCD (r1!1) (r2!1)
                    tr        = tt × delCol 1 (fromRows [r1,r2])
                    (r1',r2') = (tr `row` 1, tr `row` 2)
                in c (addRow (a.>r1') 1 nz')
                     (addRow (0.>r2') 1 rest)

-- | extGCD a b --> (gcd(a,b), tt)  
--   a,b are divided repeatedly with remainder, like in 
--   extended gcd method. tt is a protocol 2x2 matrix
--   so, [a,b] ·× tt = [gcd(a,b),0]
extGCD :: (Num α, Integral α) => α -> α -> (α, SparseMatrix α)
extGCD a b = egcd a b (idMx 2)
    where egcd a b tt =
            let (q,r) = divMod b a  -- quotRem ???
                (row1, row2) = (tt `row` 1, tt `row` 2)
                row2' = if q == 0 then row2 
                                 else row2 - (fmap (q*) row1)
            in if r /= 0
                  then egcd r a (fromRows [row2', row1])
                  else (a, fromRows [row1, row2'])
