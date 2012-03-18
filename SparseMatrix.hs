module SparseMatrix 
-- TODO: explicit export list
where

import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))

--------------------------------------------------------------------------------
-- SPARSE VECTOR DATATYPE --
----------------------------

type Index = Int

type SVec α = IntMap α

-- | Sparse vector is just indexed map of non-zero values
data SparseVector α = SV
     { dim :: Int    -- ^ real size of filled vector
     , vec :: SVec α -- ^ IntMap repsesenting non-zero values
     } deriving Eq

instance Functor SparseVector where
    fmap f v = v {vec = fmap f (vec v)}

instance Foldable SparseVector where
    foldr f d v = F.foldr f d (vec v)

--------------------------------------------------------------------------------
-- SPARSE MATRIX DATATYPE --
----------------------------

type SMx  α = SVec (SVec α)

-- | Sparse matrix is indexed map of non-zero rows, 
data SparseMatrix α = SM 
     { dims :: (Int,Int) -- ^ real height and width of filled matrix
     , mx   :: SMx α     -- ^ IntMap (IntMap α) representing non-zero values
     } deriving Eq

instance Functor SparseMatrix where
    fmap f m = m {mx = fmap (fmap f) (mx m)}

-- | Matrix real height and width
height, width :: SparseMatrix α -> Int
height = fst . dims
width  = snd . dims

-- for internal use
emptyMx = SM (0,0) M.empty

--------------------------------------------------------------------------------
-- LOOKUP/UPDATE --
-------------------

-- | Looks up an element in the vector (if not found, zero is returned)
(!) :: (Num α) => SparseVector α -> Index -> α
v ! i = findWithDefault 0 i (vec v)

-- | Looks up an element in the vector (if not found, zero is returned)
(#) :: (Num α) => SparseMatrix α -> (Index,Index) -> α
m # (i,j) = maybe 0 (findWithDefault 0 j) (M.lookup i (mx m))


-- | Erases matrix element at given index
erase :: SparseMatrix α -> (Index,Index) -> SparseMatrix α 
m `erase` (i,j) = m { mx = newMx }
    where newMx = M.adjust (M.delete j) i (mx m)

-- | Inserts new element to the sparse matrix (replaces old value)
ins :: (Num α) => SparseMatrix α -> ((Index,Index), α) -> SparseMatrix α 
m `ins` ((i,j),0) = m `erase` (i,j)
m `ins` ((i,j),x) = m { mx = newMx }
    where newMx = M.insertWith' M.union i (M.singleton j x) (mx m)

--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

-----------
-- Vectors:

-- | Returns plain list with all zeroes restored
fillVec :: (Num α) => SparseVector α -> [α]
fillVec v = [ v ! i | i <- [1 .. dim v] ]

-- | Converts plain list to sparse vector, throwing out all zeroes
sparseList :: (Num α) => [α] -> SparseVector α
sparseList l = SV (length l) $ M.fromList [ (i,x) | (i,x) <- zip [1..] l, x /= 0 ]

instance (Show α, Num α) => Show (SparseVector α) where
    show = showSparseList . fillVec

showSparseList :: (Show α, Num α) => [α] -> String
showSparseList l = show (length l)++":  ["++
    (intercalate "|" $ L.map showNonZero l)++"]"

showNonZero x  = if x == 0 then " " else show x
               
------------
-- Matrices:

-- ? loosing information about matrix size
--
-- | Converts sparse matrix to associative list,
--   adding fake zero element, to save real size for inverse conversion
toAssocList :: (Num α) => SparseMatrix α -> [ ((Index,Index), α) ]
toAssocList (SM s m) = (s, 0) : 
    [ ((i,j), x) | (i,row) <- M.toList m, (j,x) <- M.toList row, x /= 0 ] 

-- | Converts associative list to sparse matrix,
--   using maximal index as matrix size
fromAssocList :: (Num α) => [ ((Index,Index), α) ] -> SparseMatrix α 
fromAssocList l = let size = L.maximum $ fmap fst l
                      m    = L.foldl ins emptyMx l
                  in m { dims = size }

-- | Converts sparse matrix to plain list-matrix with all zeroes restored
fillMx :: (Num α) => SparseMatrix α -> [[α]]
fillMx m = [ [ m # (i,j) | j <- [1 .. width  m] ]
                         | i <- [1 .. height m] ]

-- | Converts plain list-matrix to sparse matrix, throwing out all zeroes
sparseMx :: (Num α) => [[α]] -> SparseMatrix α
sparseMx [] = emptyMx
sparseMx m@(r:_) = SM (length m, length r) $ M.fromList
    [ (i,row) | (i,row) <- zipWith pair [1..] m, not (M.null row) ]
    where pair i r = (i, vec (sparseList r))

instance (Show α, Num α) => Show (SparseMatrix α) where
    show = showSparseMatrix . fillMx
    -- show (dims m) ++ show (mx m) -- 

showSparseMatrix :: (Show α, Num α) => [[α]] -> String
showSparseMatrix [] = "[]"
showSparseMatrix m = show (length m, length (head m))++": \n"++
    (unlines $ L.map (("["++) . (++"]") . intercalate "|") 
             $ transpose $ L.map column $ transpose m)

column :: (Show α, Num α) => [α] -> [String]
column c = let c'       = L.map showNonZero c
               width    = L.maximum $ L.map length c'
               offset x = replicate (width - (length x)) ' ' ++ x
           in L.map offset c'

--------------------------------------------------------------------------------
-- TRANSPOSITION --
-------------------

trans :: (Num α) => SparseMatrix α -> SparseMatrix α
trans m = let indexes = [ (i,j) | i <- [1 .. height m], j <- [1 .. width m] ]
              add acc (i,j) = acc `ins` ((j,i), m # (i,j))
              mt = F.foldl' add emptyMx indexes
          in mt { dims = (width m, height m) }

--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

-- | Dot product of two `IntMap`s (for internal use)
(··) :: (Num α) => SVec α -> SVec α -> α
v ·· w = M.foldlWithKey' f 0 v
    where f acc 0 _ = acc
          f acc i x = acc + ((findWithDefault 0 i w) * x)

-- | Dot product of two sparse vectors
dot :: (Num α) => SparseVector α -> SparseVector α -> α
dot = (·)
-- | Unicode alias for `dot`
(·) :: (Num α) => SparseVector α -> SparseVector α -> α
(SV n v) · (SV m w) | n < m     = v ·· w  -- uses shorter vector
                    | otherwise = w ·· v

-- | Matrix-by-vector multiplication
mulMV :: (Num α) => SparseMatrix α -> SparseVector α -> SparseVector α
mulMV = (×·)
-- | Unicode alias for `mulMV`
(×·)  :: (Num α) => SparseMatrix α -> SparseVector α -> SparseVector α
(SM (h,_) m) ×· (SV _ v) = SV h (M.map (v··) m)  -- dot-p v with each row

-- | Vector-by-matrix multiplication
mulVM :: (Num α) => SparseVector α -> SparseMatrix α -> SparseVector α
mulVM = (·×)
-- | Unicode alias for `mulVM`
(·×)  :: (Num α) => SparseVector α -> SparseMatrix α -> SparseVector α
v ·× m = (trans m) ×· v

-- | Sparse matrices multiplication
mul :: (Num α) => SparseMatrix α -> SparseMatrix α -> SparseMatrix α
mul = (×)
-- | Unicode alias for `mul`
(×) :: (Num α) => SparseMatrix α -> SparseMatrix α -> SparseMatrix α
a × b = let d  = (height a, width b)    -- size of result
            bt = mx (trans b)           -- columns of b
            m  = M.map (\aRow -> M.map (aRow··) bt) (mx a)
            -- each row of a should be dot-multiplied on b columns
        in SM d m
