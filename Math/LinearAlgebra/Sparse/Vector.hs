-- | This module provides common funtions for manipulating sparse vectors
module Math.LinearAlgebra.Sparse.Vector
(
-- ** Sparse vector datatype

Index, SVec, SparseVector (..),

-- ** Basic functions

setLength, emptyVec, zeroVec, isZeroVec, isNotZeroVec, singVec,

-- ** Filter

partitionVec,

-- ** Lookup/update

(!), eraseInVec, vecIns,

-- ** Vectors combining

unionVecsWith, intersectVecsWith,

-- ** To/from list

fillVec, sparseList, vecToAssocList, vecFromAssocListWithSize, vecFromAssocList,

-- ** Multiplications

dot, (·)

)
where

import Math.LinearAlgebra.Sparse.IntMapUtilities

import Data.Functor
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))
import Data.Monoid


--------------------------------------------------------------------------------
-- SPARSE VECTOR DATATYPE
----------------------------

type Index = Int

-- | Type of internal vector storage
type SVec α = IntMap α

-- | Sparse vector is just indexed map of non-zero values
data SparseVector α = SV
     { dim :: Int    -- ^ real size of vector (with zeroes)
     , vec :: SVec α -- ^ IntMap storing non-zero values
     } deriving Eq

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
--   `fromInteger` constructs a vector with single element (only @fromInteger 0 = `emptyVec`@). So, 
--
-- >>> 3 + (sparseList [0,2,1])
-- sparseList [3,2,1]
--
instance (Eq α, Num α) => Num (SparseVector α) where
    -- (SV n v) + (SV m w) = SV (max n m) (unionWith (+) v w)
    -- (SV n v) * (SV m w) = SV (max n m) (intersectionWith (*) v w)
    (+)           = unionVecsWith (+)
    (*)           = intersectVecsWith (*)
    negate        = fmap negate
    fromInteger 0 = emptyVec
    fromInteger x = singVec (fromInteger x)
    abs           = fmap abs
    signum        = fmap signum

-- | Monoid `mappend` operation works like concatenation of two vectors 
--   (indexes of second are shifted by length of first)
--
--   Examples:
--
-- >>> (sparseList [0,1,0,2]) <> (sparseList [3,0,4])
-- sparseList [0,1,0,2,3,0,4]
--
-- >>> 1 <> (sparseList [2,3])
-- sparseList [1,2,3]
--
instance Monoid (SparseVector α) where
    mempty = emptyVec
    (SV n v) `mappend` (SV m w) = SV (n+m) (v `M.union` (shiftKeys n w))

-- | Shows size and filled vector (but without zeroes)
instance (Show α, Eq α, Num α) => Show (SparseVector α) where
    show = showSparseList . fillVec

showSparseList :: (Show α, Eq α, Num α) => [α] -> String
showSparseList l = show (length l)++":  ["++
    (intercalate "|" $ L.map showNonZero l)++"]"

showNonZero x  = if x == 0 then " " else show x

--------------------------------------------------------------------------------
-- Basic functions --
---------------------

-- | Sets vector's size
setLength ::  Int -> SparseVector α -> SparseVector α
setLength n v = v { dim = n }

-- | Vector of zero size with no values
emptyVec :: SparseVector α
emptyVec = SV 0 M.empty

-- | Vector of given size with no non-zero values
zeroVec ::  Int -> SparseVector α
zeroVec n = setLength n emptyVec

-- | Checks if vector has no non-zero values (i.e. is empty)
isZeroVec, isNotZeroVec :: SparseVector α -> Bool
isZeroVec = M.null . vec
isNotZeroVec = not . isZeroVec

-- | Vector of length 1 with given value
singVec :: (Eq α, Num α) => α -> SparseVector α
singVec 0 = zeroVec 1
singVec x = SV 1 (singleton 1 x)

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

-- | Updates value at given index
vecIns :: (Eq α, Num α) => SparseVector α -> (Index, α) -> SparseVector α
v `vecIns` (i,0) = v `eraseInVec` i
v `vecIns` (i,x) = v { vec = M.insert i x (vec v) }

--------------------------------------------------------------------------------
-- VECTORS COMBINING --
-----------------------

-- | Unions non-zero values of vectors and applies given function on intersection
unionVecsWith :: (α -> α -> α) -> SparseVector α -> SparseVector α -> SparseVector α
unionVecsWith f (SV x v) (SV y w) = SV (max x y) $ M.unionWith f v w

-- | Intersects non-zero values of vectors and applies given function on them
intersectVecsWith :: (α -> α -> α) -> SparseVector α -> SparseVector α -> SparseVector α
intersectVecsWith f (SV x v) (SV y w) = SV (max x y) $ M.intersectionWith f v w

--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

-- | Returns plain list with all zeroes restored
fillVec :: (Num α) => SparseVector α -> [α]
fillVec v = [ v ! i | i <- [1 .. dim v] ]

-- | Converts plain list to sparse vector, throwing out all zeroes
sparseList :: (Num α, Eq α) => [α] -> SparseVector α
sparseList l = SV (length l) $ M.fromList [ (i,x) | (i,x) <- zip [1..] l, x /= 0 ]

-- | Converts sparse vector to an associative list,
--   adding fake zero element, to save real size for inverse conversion
vecToAssocList :: (Num α, Eq α) => SparseVector α -> [ (Index, α) ]
vecToAssocList v = (dim v, 0) : (M.toAscList (vec v))

-- | Converts associative list to sparse vector,
--   of given size
vecFromAssocListWithSize :: (Num α, Eq α) => Int -> [ (Index, α) ] -> SparseVector α
vecFromAssocListWithSize s l = L.foldl' vecIns (zeroVec s) l

-- | Converts associative list to sparse vector,
--   using maximal index as it's size
vecFromAssocList :: (Num α, Eq α) => [ (Index, α) ] -> SparseVector α
vecFromAssocList l = vecFromAssocListWithSize (L.maximum $ fmap fst l) l

--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

-- | Dot product of two sparse vectors
dot :: (Eq α, Num α) => SparseVector α -> SparseVector α -> α
dot = (·)

-- | Unicode alias for `dot`
(·) :: (Eq α, Num α) => SparseVector α -> SparseVector α -> α
v · w = (vec v) ·· (vec w)
