module Math.LinearAlgebra.Sparse.Matrix 
-- TODO: explicit export list
where

import Data.Functor
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))
import Data.Monoid

import Math.LinearAlgebra.Sparse.Vector

--------------------------------------------------------------------------------
-- SPARSE MATRIX DATATYPE --
----------------------------

-- | Internal storage of matrix
type SMx  α = SVec (SVec α)

-- | Sparse matrix is indexed map of non-zero rows, 
data SparseMatrix α = SM 
     { dims :: (Int,Int) -- ^ real height and width of filled matrix
     , mx   :: SMx α     -- ^ IntMap (IntMap α) representing non-zero values
     } deriving Eq

-- | `fmap` applies given function on all non-zero values
instance Functor SparseMatrix where
    fmap f m = m {mx = fmap (fmap f) (mx m)}

-- instance Foldable SparseMatrix where
--     foldr f d m = F.foldr f d $ fmap (F.foldr f d) (mx m)

instance (Eq α, Num α) => Num (SparseMatrix α) where
    (SM (h1,w1) m) + (SM (h2,w2) n) 
        = SM (max h1 h2, max w1 w2) $ M.filter (not . M.null)
             (unionWith (unionWith (+)) m n)
    (SM (h1,w1) m) * (SM (h2,w2) n)
        = SM (max h1 h2, max w1 w2) $ M.filter (not . M.null)
             (intersectionWith (intersectionWith (*)) m n)
    negate         = fmap negate
    fromInteger x  = diagonalMx [fromInteger x]
    abs            = fmap abs
    signum         = fmap signum

-- | Matrix real height and width
height, width :: SparseMatrix α -> Int
height = fst . dims
width  = snd . dims

-- | Sets height and width of matrix
setSize :: (Num α) => (Int,Int) -> SparseMatrix α -> SparseMatrix α
setSize s m = m { dims = s }

-- | Matrix of zero size with no values
emptyMx ::  SparseMatrix α
emptyMx = SM (0,0) M.empty

-- | Zero matrix of given size
zeroMx ::  Num α => (Int, Int) -> SparseMatrix α
zeroMx (h,w) = setSize (h,w) emptyMx

-- | Checks if vector has no non-zero values (i.e. is empty)
isZeroVec, isNotZeroVec :: SparseVector α -> Bool
isZeroVec = M.null . vec
isNotZeroVec = not . isZeroVec

-- | Checks if matrix has no non-zero values (i.e. is empty)
isZeroMx, isNotZeroMx  :: SparseMatrix α -> Bool
isZeroMx = M.null . mx
isNotZeroMx = not . isZeroMx

-- | Identity matrix of given size
idMx :: (Num α, Eq α) => Int -> SparseMatrix α
idMx n = diagonalMx (L.replicate n 1)

--------------------------------------------------------------------------------
-- ADDING/DELETING ROW/COLUMNS --
---------------------------------

-- | Adds row at given index, increasing matrix height by 1 
--   and shifting indexes after it
addRow :: (Num α) => SparseVector α -> Index -> SparseMatrix α -> SparseMatrix α
addRow v i m 
    | isZeroMx m = SM (1, dim v)
                      (singleton 1 (vec v))
    | otherwise  = SM (height m + 1, max (width m) (dim v))
                      (addElem mbv i (mx m))
    where mbv = if isZeroVec v then Nothing else Just (vec v)

-- | Adds column at given index, increasing matrix width by 1 
--   and shifting indexes after it
addCol :: (Num α) => SparseVector α -> Index -> SparseMatrix α -> SparseMatrix α
addCol v j m 
    | isZeroMx m = SM (dim v, 1)
                      (M.map (singleton 1) (vec v))
    | otherwise  = SM (max (height m) (dim v), width m + 1)
                      (M.mapWithKey insCol (mx m))
    where insCol i row = addElem (M.lookup i (vec v)) j row

-- | Just adds zero row at given index
addZeroRow ::  Num α => Index -> SparseMatrix α -> SparseMatrix α
addZeroRow i m = addRow (zeroVec (width m)) i m

-- | Just adds zero column at given index
addZeroCol ::  Num α => Index -> SparseMatrix α -> SparseMatrix α
addZeroCol i m = addCol (zeroVec (height m)) i m

-- | Deletes row at given index, decreasing matrix height by 1 
--   and shifting indexes after it
delRow :: (Num α) => Index -> SparseMatrix α -> SparseMatrix α
delRow i m | isZeroMx m = m
           | otherwise  = SM (height m - 1, width m)
                             (delElem i (mx m))

-- | Deletes column at given index, decreasing matrix width by 1 
--   and shifting indexes after it
delCol :: (Num α) => Index -> SparseMatrix α -> SparseMatrix α
delCol j m | isZeroMx m = m
           | otherwise  = SM (height m, width m - 1)
                             (M.map (delElem j) (mx m))

-- | Deletes row and column at given indexes
delRowCol :: Num α => Index -> Index -> SparseMatrix α -> SparseMatrix α
delRowCol i j m = delCol j (delRow i m)

-- | Partitions matrix, using pedicate on rows and returns two /new/ matrices,
--   one constructed from rows satisfying predicate, and another from rows that don't
partitionMx :: (Num α) => (SparseVector α -> Bool) -> SparseMatrix α -> (SparseMatrix α, SparseMatrix α)
partitionMx p (SM (h,w) m) = (SM (st,w) t, SM (h-st,w) f)
    where (t,f) = partitionMap (p . SV w) m
          st = size t

--------------------------------------------------------------------------------
-- LOOKUP/UPDATE --
-------------------

-- | Looks up an element in the matrix (if not found, zero is returned)
(#) :: (Num α) => SparseMatrix α -> (Index,Index) -> α
m # (i,j) = maybe 0 (findWithDefault 0 j) (M.lookup i (mx m))

-- | Returns row of matrix at given index
row :: (Num α) => SparseMatrix α -> Index -> SparseVector α
m `row` i = SV (width m) (findWithDefault M.empty i (mx m))

-- | Returns column of matrix at given index
col :: (Num α, Eq α) => SparseMatrix α -> Index -> SparseVector α
m `col` i = (trans m) `row` i
-- TODO: transpositioning whole matrix is not effective

-- | Updates values in row using given function
updRow :: (Num α) => SparseMatrix α -> (SparseVector α -> SparseVector α) -> Index -> SparseMatrix α
updRow m f i = m { mx = M.adjust f' i (mx m) }
    where f' = vec . f . SV (width m)

-- | Fills row with zeroes (i.e. deletes it, but size of matrix doesn't change)
eraseRow :: (Num α) => SparseMatrix α -> Index -> SparseMatrix α
m `eraseRow` i = m { mx = M.delete i (mx m) }

-- | Erases matrix element at given index
erase :: (Num α) => SparseMatrix α -> (Index,Index) -> SparseMatrix α 
m `erase` (i,j) = if isZeroVec (m' `row` i)     -- if that was the last element
                     then m' `eraseRow` i       -- delete this row
                     else m'
    where m' = updRow m (`eraseInVec` j) i

-- | Inserts new element to the sparse matrix (replaces old value)
ins :: (Num α, Eq α) => SparseMatrix α -> ((Index,Index), α) -> SparseMatrix α 
m `ins` ((i,j),0) = m `erase` (i,j)
m `ins` ((i,j),x) = m { mx = newMx }
    where newMx = M.insertWith' M.union i (M.singleton j x) (mx m)

--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

------------
-- Matrices:

-- | Constructs square matrix with given elements on diagonal
diagonalMx :: (Num α, Eq α) => [α] -> SparseMatrix α
diagonalMx = L.foldl add emptyMx
    where add m x = let i = height m + 1
                    in setSize (i,i) (m `ins` ((i,i),x))

-- | Constructs matrix from a list of rows
fromRows :: (Num α) => [SparseVector α] -> SparseMatrix α
fromRows = L.foldl (\m r -> addRow r (height m + 1) m) emptyMx

-- | Converts sparse matrix to associative list,
--   adding fake zero element, to save real size for inverse conversion
toAssocList :: (Num α, Eq α) => SparseMatrix α -> [ ((Index,Index), α) ]
toAssocList (SM s m) = (s, 0) : 
    [ ((i,j), x) | (i,row) <- M.toAscList m, (j,x) <- M.toAscList row, x /= 0 ] 

-- | Converts associative list to sparse matrix,
--   using maximal index as matrix size
fromAssocList :: (Num α, Eq α) => [ ((Index,Index), α) ] -> SparseMatrix α 
fromAssocList l = let size = L.maximum $ fmap fst l
                      m    = L.foldl ins emptyMx l
                  in m { dims = size }

-- | Converts sparse matrix to plain list-matrix with all zeroes restored
fillMx :: (Num α) => SparseMatrix α -> [[α]]
fillMx m = [ [ m # (i,j) | j <- [1 .. width  m] ]
                         | i <- [1 .. height m] ]

-- | Converts plain list-matrix to sparse matrix, throwing out all zeroes
sparseMx :: (Num α, Eq α) => [[α]] -> SparseMatrix α
sparseMx [] = emptyMx
sparseMx m@(r:_) = SM (length m, length r) $ M.fromList
    [ (i,row) | (i,row) <- zipWith pair [1..] m, not (M.null row) ]
    where pair i r = (i, vec (sparseList r))

-- | Shows size and filled matrix (but without zeroes)
instance (Show α, Eq α, Num α) => Show (SparseMatrix α) where
    show = showSparseMatrix . fillMx
    -- show (dims m) ++ show (mx m) -- 

showSparseMatrix :: (Show α, Eq α, Num α) => [[α]] -> String
showSparseMatrix [] = "[]"
showSparseMatrix m = show (length m, length (head m))++": \n"++
    (unlines $ L.map (("["++) . (++"]") . intercalate "|") 
             $ transpose $ L.map column $ transpose m)

column :: (Show α, Eq α, Num α) => [α] -> [String]
column c = let c'       = L.map showNonZero c
               width    = L.maximum $ L.map length c'
               offset x = replicate (width - (length x)) ' ' ++ x
           in L.map offset c'

--------------------------------------------------------------------------------
-- TRANSPOSITION --
-------------------

-- | Transposes matrix (rows became columns)
trans :: (Num α, Eq α) => SparseMatrix α -> SparseMatrix α
trans m = let indexes = [ (i,j) | i <- [1 .. height m], j <- [1 .. width m] ]
              add acc (i,j) = acc `ins` ((j,i), m # (i,j))
              mt = F.foldl' add emptyMx indexes
          in mt { dims = (width m, height m) }

--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

-- | Matrix-by-vector multiplication
mulMV :: (Num α, Eq α) => SparseMatrix α -> SparseVector α -> SparseVector α
mulMV = (×·)
-- | Unicode alias for `mulMV`
(×·)  :: (Num α, Eq α) => SparseMatrix α -> SparseVector α -> SparseVector α
(SM (h,_) m) ×· (SV _ v) = SV h (M.filter (0/=) (M.map (v··) m))  -- dot-p v with each row

-- | Vector-by-matrix multiplication
mulVM :: (Num α, Eq α) => SparseVector α -> SparseMatrix α -> SparseVector α
mulVM = (·×)
-- | Unicode alias for `mulVM`
(·×)  :: (Num α, Eq α) => SparseVector α -> SparseMatrix α -> SparseVector α
v ·× m = (trans m) ×· v

-- | Sparse matrices multiplication
mul :: (Num α, Eq α) => SparseMatrix α -> SparseMatrix α -> SparseMatrix α
mul = (×)
-- | Unicode alias for `mul`
(×) :: (Num α, Eq α) => SparseMatrix α -> SparseMatrix α -> SparseMatrix α
a × b = let d  = (height a, width b)    -- size of result
            bt = mx (trans b)           -- columns of b
            m  = M.filter (not . M.null)
               $ M.map (\aRow -> M.filter (0/=) (M.map (aRow··) bt)) (mx a)
            -- each row of a should be dot-multiplied on b columns
        in SM d m
