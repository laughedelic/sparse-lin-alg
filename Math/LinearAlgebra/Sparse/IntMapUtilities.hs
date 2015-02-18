-- | This module provides common funtions for manipulating IntMap datatype in sparse vectors and matrices
module Math.LinearAlgebra.Sparse.IntMapUtilities
where

import Data.List     as L
import Data.IntMap   as M

--------------------------------------------------------------------------------
-- IntMap Utilities (for internal use) --
-----------------------------------------

-- | Dot product of two `IntMap`s (for internal use)
(··) :: (Num α, Eq α) => IntMap α -> IntMap α -> Maybe α
v ·· w = case M.foldl' (+) 0 $ M.intersectionWith (*) v w of
  0 -> Nothing
  x -> Just x
--    M.foldlWithKey' f 0 v
--    where f acc 0 _ = acc
--          f acc i x = acc + ((findWithDefault 0 i w) * x)

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
