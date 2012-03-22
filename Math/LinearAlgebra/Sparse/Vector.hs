module Math.LinearAlgebra.Sparse.Vector
-- TODO: explicit export list
where

import Data.Functor
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))
import Data.Monoid

--------------------------------------------------------------------------------
-- IntMap Utilities (for internal use) --
-----------------------------------------

-- | Dot product of two `IntMap`s (for internal use)
(··) :: (Num α) => IntMap α -> IntMap α -> α
v ·· w = M.foldlWithKey' f 0 v
    where f acc 0 _ = acc
          f acc i x = acc + ((findWithDefault 0 i w) * x)

-- | Shifts (re-enumerates) keys of IntMap by given number
shiftKeys :: Int -> IntMap α -> IntMap α
shiftKeys k m = M.fromAscList [ (i+k,x) | (i,x) <- M.toAscList m ]

-- | Adds element to the map at given index, shifting all keys after it
addElem :: Maybe α -> Int -> IntMap α -> IntMap α
addElem v i m = M.unions [a, maybeSingleton i v, maybeSingleton (i+1) x, shiftKeys 1 b]
     where (a, x, b) = M.splitLookup i m
           maybeSingleton k = maybe M.empty (singleton k)

-- | Deletes element of the map at given index, shifting all keys after it
delElem :: Int -> IntMap α -> IntMap α
delElem i m = a `M.union` (shiftKeys (-1) b)
    where (a,b) = M.split i m

-- | Splits map using predicate and returns a pair with filtered map and
--   re-enumereted second part (that doesn't satisfy predicate). For example:
--
-- >>> partitionMap (>0) (fromList [(1,1),(2,-1),(3,-2),(4,3),(5,-4)])
-- ( fromList [(1,1),(4,3)], fromList [(1,-1),(2,-2),(3,-4)] )
-- 
partitionMap :: (α -> Bool) -> IntMap α -> (IntMap α, IntMap α)
partitionMap p m = (m', f')
    where f  = M.filter (not . p) m
          f' = M.fromAscList $ zip [1..] (M.elems f)
          m' = L.foldl (\mm j -> delElem j mm) m $ L.reverse (M.keys f)

--------------------------------------------------------------------------------
-- SPARSE VECTOR DATATYPE --
----------------------------

type Index = Int

-- | Type of internal vector storage
type SVec α = IntMap α

-- | Sparse vector is just indexed map of non-zero values
data SparseVector α = SV
     { dim :: Int    -- ^ real size of vector (with zeroes)
     , vec :: SVec α -- ^ IntMap storing non-zero values
     } deriving Eq

-- | Sets vector's size
setLength ::  Int -> SparseVector α -> SparseVector α
setLength n v = v { dim = n }

-- | Vector of zero size with no values
emptyVec :: SparseVector α
emptyVec = SV 0 M.empty

-- | Vector of given size with no non-zero values
zeroVec ::  Int -> SparseVector α
zeroVec n = setLength n emptyVec

-- | Vector of length 1 with given value
singVec :: (Eq α, Num α) => α -> SparseVector α
singVec 0 = zeroVec 1
singVec x = SV 1 (singleton 1 x)

-- | `fmap` applies given function on vector non-zero values
instance Functor SparseVector where
    fmap f v = v {vec = fmap f (vec v)}

-- | fold functions are applied to non-zero values
instance Foldable SparseVector where
    foldr f d v = F.foldr f d (vec v)

-- | `Num` operators like @(*)@, @(+)@ and @(-)@ work on sparse vectors 
--   like @`zipWith` (…)@ works on lists, except size of result is maximum
--   of arguments sizes.
--
--   `signum`, `abs` and `negate` work through `fmap`, so usual `Num` laws
--   are satisfied (such as @(signum v)*(abs v) == v@.
--
--   `fromInteger` constructs a vector with single element. So, 
--
-- >>> 3 + (sparseList [0,2,1])
-- sparseList [3,2,1]
--
instance (Eq α, Num α) => Num (SparseVector α) where
    (SV n v) + (SV m w) = SV (max n m) (unionWith (+) v w)
    (SV n v) * (SV m w) = SV (max n m) (intersectionWith (*) v w)
    negate              = fmap negate
    fromInteger x       = singVec (fromInteger x)
    abs                 = fmap abs
    signum              = fmap signum

-- | Monoid `mappend` operation works like concatenation of two vectors 
--   (indexes of second are shifted by length of first)
instance Monoid (SparseVector α) where
    (SV n v) `mappend` (SV m w) = SV (n+m) (v `M.union` (shiftKeys n w))
    mempty = emptyVec

-- | This is like cons (`:`) operator for lists.
--
--   @x .> v = singVec x \<\> v@
--
(.>) :: (Eq α, Num α) => α -> SparseVector α -> SparseVector α
x .> v = singVec x <> v

--------------------------------------------------------------------------------
-- FILTER --
------------

-- | Splits vector using predicate and returns a pair with filtered values and
--   re-enumereted second part (that doesn't satisfy predicate). For example:
--
-- >>> partitionVec (>0) (sparseList [0,1,-1,2,3,0,-4,5,-6,0,7])
-- ( sparseList [0,1,0,2,3,0,0,5,0,0,7], sparseList [-1,-4,-6] )
-- 
partitionVec :: (Num α) => (α -> Bool) -> SparseVector α -> (SparseVector α, SparseVector α)
partitionVec p (SV d v) = (SV st t, SV (d-st) f)
    where (t,f) = partitionMap p v
          st = size t

--------------------------------------------------------------------------------
-- LOOKUP/UPDATE --
-------------------

-- | Looks up an element in the vector (if not found, zero is returned)
(!) :: (Num α) => SparseVector α -> Index -> α
v ! i = findWithDefault 0 i (vec v)

-- | Deletes element of vector at given index (size of vector doesn't change)
eraseInVec :: (Num α) => SparseVector α -> Index -> SparseVector α
v `eraseInVec` j = v { vec = M.delete j (vec v) }

--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

-----------
-- Vectors:

-- | Returns plain list with all zeroes restored
fillVec :: (Num α) => SparseVector α -> [α]
fillVec v = [ v ! i | i <- [1 .. dim v] ]

-- | Converts plain list to sparse vector, throwing out all zeroes
sparseList :: (Num α, Eq α) => [α] -> SparseVector α
sparseList l = SV (length l) $ M.fromList [ (i,x) | (i,x) <- zip [1..] l, x /= 0 ]

-- | Shows size and filled vector (but without zeroes)
instance (Show α, Eq α, Num α) => Show (SparseVector α) where
    show = showSparseList . fillVec

showSparseList :: (Show α, Eq α, Num α) => [α] -> String
showSparseList l = show (length l)++":  ["++
    (intercalate "|" $ L.map showNonZero l)++"]"

showNonZero x  = if x == 0 then " " else show x
               
--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

-- | Dot product of two sparse vectors
dot :: (Num α) => SparseVector α -> SparseVector α -> α
dot = (·)

-- | Unicode alias for `dot`
(·) :: (Num α) => SparseVector α -> SparseVector α -> α
(SV n v) · (SV m w) | n < m     = v ·· w  -- uses shorter vector
                    | otherwise = w ·· v
