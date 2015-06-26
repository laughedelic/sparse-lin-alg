module Math.LinearAlgebra.Sparse.Matrix
(

-- ** Sparse matrix datatype

SMx, SparseMatrix (..),

-- ** Basic functions

height, width, setSize, emptyMx, zeroMx, isZeroMx, isNotZeroMx , idMx,

-- ** Combining matrices

(//), hconcat, vconcat, sizedBlockMx, sizedBlockSMx, blockMx, blockSMx,

-- ** Adding\/deleting row\/columns

addRow, addCol, addZeroRow, addZeroCol, delRow, delCol, delRowCol, separateMx,

-- ** Lookup\/update

(#), row, col, updRow, eraseRow, erase, ins, findRowIndices, findRowIndicesR, popRow, (|>), (<|), replaceRow, exchangeRows, mapOnRows,
rows, columns,

-- ** To\/from list

diagonalMx, mainDiag, fromRows, toAssocList, fromAssocListWithSize, fromAssocList, fillMx, sparseMx,

-- ** Transposition

trans,

-- ** Multiplications

mulMV, (×·) , mulVM, (·×) , mul, (×),

)
where

import Math.LinearAlgebra.Sparse.IntMapUtilities

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

-- | All `Num` work on sparse matrices the same way as
--   on `SparseVector` (see documentation there)
instance (Eq α, Num α) => Num (SparseMatrix α) where
    (SM (h1,w1) m) + (SM (h2,w2) n) 
        = SM (max h1 h2, max w1 w2) $ M.filter (not . M.null)
             (unionWith (unionWith (+)) m n)
    (SM (h1,w1) m) * (SM (h2,w2) n)
        = SM (max h1 h2, max w1 w2) $ M.filter (not . M.null)
             (intersectionWith (intersectionWith (*)) m n)
    negate         = fmap negate
    fromInteger 0  = emptyMx
    fromInteger x  = diagonalMx [fromInteger x]
    abs            = fmap abs
    signum         = fmap signum

-- | `mempty` is just `emptyMx`
--
--   `mappend` is horisontal concatenation
instance Monoid (SparseMatrix α) where
    mempty = emptyMx
    (SM (h1,w1) m) `mappend` (SM (h2,w2) n)
        = SM (max h1 h2, w1 + w2) (M.unionWith M.union m (M.map (shiftKeys w1) n))

-- | Shows size and filled matrix (but without zeroes)
instance (Show α, Eq α, Num α) => Show (SparseMatrix α) where
    show = showSparseMatrix . fillMx

showSparseMatrix :: (Show α, Eq α, Num α) => [[α]] -> String
showSparseMatrix [] = "(0,0):\n[]\n"
showSparseMatrix m = show (length m, length (head m))++": \n"++
    (unlines $ L.map (("["++) . (++"]") . intercalate "|")
             $ transpose $ L.map column $ transpose m)

column :: (Show α, Eq α, Num α) => [α] -> [String]
column c = let c'       = L.map showNonZero c
               width    = L.maximum $ L.map length c'
               offset x = replicate (width - (length x)) ' ' ++ x
           in L.map offset c'

showNonZero x  = if x == 0 then " " else show x

--------------------------------------------------------------------------------
-- BASIC FUNCTIONS --
---------------------

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

-- | Checks if matrix has no non-zero values (i.e. is empty)
isZeroMx, isNotZeroMx  :: SparseMatrix α -> Bool
isZeroMx = M.null . mx
isNotZeroMx = not . isZeroMx

-- | Identity matrix of given size
idMx :: (Num α, Eq α) => Int -> SparseMatrix α
idMx n = diagonalMx (L.replicate n 1)

--------------------------------------------------------------------------------
-- COMBINING MATRICES --
------------------------

-- | Vertical concatenation
(//) :: SparseMatrix α -> SparseMatrix α -> SparseMatrix α
(SM (h1,w1) m) // (SM (h2,w2) n) =
    SM (h1 + h2, max w1 w2) (m `M.union` (shiftKeys h1 n))

-- | Batch horisontal\/vertical concatenation
hconcat, vconcat :: [SparseMatrix α] -> SparseMatrix α
hconcat = L.foldl' (<>) emptyMx
vconcat = L.foldl' (//) emptyMx

-- | Takes size of each block and matrix of sparse matrices
--   and constructs sparse matrix from this blocks
sizedBlockMx :: Num α => (Int, Int) -> [[SparseMatrix α]] -> SparseMatrix α
sizedBlockMx s = blockMx . fmap (fmap (setSize s))

-- | Fills sparse matrix of blocks and then applies `sizedBlockMx`
sizedBlockSMx :: (Eq α, Num α) =>(Int, Int) -> SparseMatrix (SparseMatrix α) -> SparseMatrix α
sizedBlockSMx s = sizedBlockMx s . fillMx

-- TODO: evaluate block sizes automaticaly
blockMx :: [[SparseMatrix α]] -> SparseMatrix α
blockMx = vconcat . fmap hconcat

blockSMx :: (Eq α, Num α) => SparseMatrix (SparseMatrix α) -> SparseMatrix α
blockSMx = blockMx . fillMx

--------------------------------------------------------------------------------
-- ADDING\/DELETING ROW\/COLUMNS --
---------------------------------

-- | Adds row at given index, increasing matrix height by 1
--   and shifting indexes after it
addRow :: (Num α) => SparseVector α -> Index -> SparseMatrix α -> SparseMatrix α
addRow v i m = SM (height m + 1, max (width m) (dim v))
                  (addElem mbv i (mx m))
    where mbv = if isZeroVec v then Nothing else Just (vec v)

-- | Adds column at given index, increasing matrix width by 1
--   and shifting indexes after it
addCol :: (Num α) => SparseVector α -> Index -> SparseMatrix α -> SparseMatrix α
addCol v j m = SM (max (height m) (dim v), width m + 1)
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
delRow i m | isZeroMx m = setSize (height m - 1, width m) m
           | otherwise  = SM (height m - 1, width m)
                             (delElem i (mx m))

-- | Deletes column at given index, decreasing matrix width by 1
--   and shifting indexes after it
delCol :: (Num α) => Index -> SparseMatrix α -> SparseMatrix α
delCol j m | isZeroMx m = setSize (height m, width m - 1) m
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
-- ^ WARNING: doesn't work with empty rows

-- | Separates matrix, using pedicate on rows and returns two matrices of the same size,
--   one only with rows satisfying predicate, and another with the rest rows
separateMx :: (Num α) => (SparseVector α -> Bool) -> SparseMatrix α -> (SparseMatrix α, SparseMatrix α)
separateMx p (SM (h,w) m) = (SM (h,w) t, SM (h,w) f)
    where (t,f) = M.partition (p . SV w) m
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
m `col` j = M.foldlWithKey' addElem (zeroVec (height m)) (mx m)
    where addElem acc i row = maybe acc (\x -> acc `vecIns` (i,x)) (M.lookup j row)
-- old, obvious variant: (trans m) `row` i
-- transpositioning the whole matrix is not effective

-- | Updates values in row using given function
updRow :: (Num α) => (SparseVector α -> SparseVector α) -> Index -> SparseMatrix α -> SparseMatrix α
updRow f i m = m { mx = M.adjust f' i (mx m) }
    where f' = vec . f . SV (width m)

-- | Fills row with zeroes (i.e. deletes it, but size of matrix doesn't change)
eraseRow :: (Num α) => Index -> SparseMatrix α -> SparseMatrix α
eraseRow i m = m { mx = M.delete i (mx m) }

-- | Erases matrix element at given index
erase :: (Num α) => SparseMatrix α -> (Index,Index) -> SparseMatrix α
m `erase` (i,j) = if isZeroVec (m' `row` i)     -- if that was the last element
                     then eraseRow i m'         -- delete this row
                     else m'
    where m' = updRow (`eraseInVec` j) i m

-- | Inserts new element to the sparse matrix (replaces old value)
ins :: (Num α, Eq α) => SparseMatrix α -> ((Index,Index), α) -> SparseMatrix α
m `ins` ((i,j),0) = m `erase` (i,j)
m `ins` ((i,j),x) = m { mx = newMx }
    where newMx = M.insertWith' M.union i (M.singleton j x) (mx m)

-- | Finds indices of rows, that satisfy given predicate. Searches from left to right (in ascending order of indices)
findRowIndices :: (SparseVector α -> Bool) -> SparseMatrix α -> [Int]
findRowIndices  p m = fst $ M.mapAccumRWithKey (\acc i x -> (if p (SV (width m) x) then i:acc else acc,x)) [] (mx m)

-- | Finds indices of rows, that satisfy given predicate. Searches from right to left (in descending order of indices)
findRowIndicesR :: (SparseVector α -> Bool) -> SparseMatrix α -> [Int]
findRowIndicesR p m = fst $ M.mapAccumWithKey  (\acc i x -> (if p (SV (width m) x) then i:acc else acc,x)) [] (mx m)

-- TODO: more effective implementation (is done in exchangeRows?)
-- moveRow i j m | i == j    = m
--               | otherwise = addRow r j $ delRow i m
--     where r = m `row` i

-- | Returns a row at given index and matrix without it
popRow :: Num α =>Index -> SparseMatrix α -> (SparseVector α, SparseMatrix α)
popRow i m = (m `row` i, delRow i m)

-- | Adds row to matrix at the top
(|>) ::  Num α => SparseVector α -> SparseMatrix α -> SparseMatrix α
r |> m = addRow r 1 m

-- | Adds row to matrix at the bottom
(<|) ::  Num α => SparseMatrix α -> SparseVector α -> SparseMatrix α
m <| r = addRow r (height m + 1) m

-- | Replaces row at given index with given vector
replaceRow :: Num α => SparseVector α -> Index -> SparseMatrix α -> SparseMatrix α
replaceRow r i m | isZeroVec r = eraseRow i m
                 | otherwise   = m { mx = M.insert i (vec r) (mx m) }

-- | Exchanges positions of two rows
exchangeRows :: Num α => Index -> Index -> SparseMatrix α -> SparseMatrix α
exchangeRows i j m | i == j    = m
                   | otherwise = replaceRow (m `row` i) j
                               $ replaceRow (m `row` j) i m

-- | Applies vector-function on matrix rows
mapOnRows :: (SparseVector α -> SparseVector β)-> SparseMatrix α -> SparseMatrix β
mapOnRows f m = m { mx = M.map (vec . f . (SV (width m))) (mx m) }

-- | Returns vector with matrix rows
rows :: SparseMatrix α -> SparseVector (SparseVector α)
rows (SM (h,w) m) = SV h (M.map (SV w) m)

-- | Returns vector with matrix columns (@rows . trans@)
columns :: (Eq α, Num α) => SparseMatrix α -> SparseVector (SparseVector α)
columns = rows . trans

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

-- | Collects main diagonal of matrix
mainDiag ::  (Eq α, Num α) => SparseMatrix α -> SparseVector α
mainDiag m = sparseList [ m#(i,i) | i <- [1 .. l] ]
    where l = min (height m) (width m)


-- | Constructs matrix from a set (list/vector/etc.) of rows
fromRows :: (Num α, Foldable φ) => φ (SparseVector α) -> SparseMatrix α
fromRows = F.foldl' (<|) emptyMx

-- | Converts sparse matrix to associative list,
--   adding fake zero element, to save real size for inverse conversion
toAssocList :: (Num α, Eq α) => SparseMatrix α -> [ ((Index,Index), α) ]
toAssocList (SM s m) = (s, 0) :
    [ ((i,j), x) | (i,row) <- M.toAscList m, (j,x) <- M.toAscList row, x /= 0 ]

-- | Converts associative list to sparse matrix,
--   of given size
fromAssocListWithSize :: (Num α, Eq α) => (Int,Int) -> [ ((Index,Index), α) ] -> SparseMatrix α
fromAssocListWithSize s l = L.foldl' ins (zeroMx s) l

-- | Converts associative list to sparse matrix,
--   using maximal index as matrix size
fromAssocList :: (Num α, Eq α) => [ ((Index,Index), α) ] -> SparseMatrix α
fromAssocList l = fromAssocListWithSize size l
    where size = L.foldl maxIndices (0, 0) l
          maxIndices (mX, mY) ((x, y), _) = (max mX x, max mY y)

-- | Converts sparse matrix to plain list-matrix with all zeroes restored
fillMx :: (Num α) => SparseMatrix α -> [[α]]
fillMx m = [ [ m # (i,j) | j <- [0 .. width  m] ]
                         | i <- [0 .. height m] ]

-- | Converts plain list-matrix to sparse matrix, throwing out all zeroes
sparseMx :: (Num α, Eq α) => [[α]] -> SparseMatrix α
sparseMx [] = emptyMx
sparseMx m@(r:_) = SM (length m, length r) $ M.fromList
    [ (i,row) | (i,row) <- zipWith pair [1..] m, not (M.null row) ]
    where pair i r = (i, vec (sparseList r))

--------------------------------------------------------------------------------
-- TRANSPOSITION --
-------------------

-- | Transposes matrix (rows become columns)
trans :: (Num α, Eq α) => SparseMatrix α -> SparseMatrix α
trans m = let mt                  = M.foldlWithKey'  accRow     emptyMx (mx m)
              accRow    acc i row = M.foldlWithKey' (accElem i) acc      row
              accElem i acc j x   = acc `ins` ((j,i),x)
          in setSize (width m, height m) mt

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
