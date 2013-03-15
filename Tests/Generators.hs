module Tests.Generators where

import Test.QuickCheck

import Control.Monad (liftM)
import Data.Function (on)
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Vector
import Math.LinearAlgebra.Sparse.Matrix

----------------------------------------------------------------------

-- | Generates non negative number
genNonNegative :: (Arbitrary α, Ord α, Eq α, Num α) => Gen α
genNonNegative = do NonNegative x <- arbitrary
                    return x

genNonZero :: (Arbitrary α, Ord α, Eq α, Num α) => Gen α
genNonZero = do NonZero x <- arbitrary
                return x

-- | Combains two generators to pair
genPair g1 g2 = do
    x <- g1
    y <- g2
    return (x,y)

genSize (m,n) = choose (m,n)

genIndex (m,n) | m > n     = return 0
               | otherwise = choose (m,n)

-- | Generates IntMap of given size with indecies in given diapason 
--   and arbitrary values
genIntMap :: Gen α       -- ^ generator for values
          -> Int         -- ^ size of IntMap
          -> (Int,Int)   -- ^ diapason of indicies
          -> Gen (IntMap α) 
genIntMap g n (a,b) = liftM M.fromList 
                   $ vectorOf n $ genPair (choose (a,b)) g

--------------------------------------------------------------------------------
-- SPARSE VECTOR DATATYPE --
----------------------------

-- | Generates sparse vector of given size with indecies in given diapason 
--   and arbitrary values
genSVec :: (Arbitrary α, Ord α, Eq α, Num α) => Int -> (Int,Int) -> Gen (SVec α)
genSVec = genIntMap genNonZero

-- | Generates `SparseVector` of given length with arbitrary values
genSparseVector :: (Arbitrary α, Ord α, Eq α, Num α) => Int -> Gen (SparseVector α)
genSparseVector m = do
    n <- choose (0,m)
    v <- genSVec n (1,m)
    return $ SV m v

-- | Arbitrary `SparseVector` of random length
instance (Arbitrary α, Eq α, Ord α, Num α) => Arbitrary (SparseVector α) where
    arbitrary = genSparseVector =<< genNonNegative

-- | Generates list of given length with frequent zeroes
genSparseList :: (Arbitrary α, Num α) => Int -> Gen [α]
genSparseList n = vectorOf n $ frequency [(4, return 0), (1, arbitrary)]

genVecAssocList :: (Arbitrary α, Ord α, Eq α, Num α) => Int -> Int -> Gen [ (Index, α) ]
genVecAssocList n s = liftM ( ((s,0):) . nubBy ((==) `on` fst) ) $
    vectorOf n $ genPair (genIndex (1,s)) genNonZero

--------------------------------------------------------------------------------
-- SPARSE MATRIX DATATYPE --
----------------------------

genSparseMatrix :: (Arbitrary α, Ord α, Eq α, Num α) => (Int,Int) -> Gen (SparseMatrix α)
genSparseMatrix (h,w) = do
    n <- choose (1,w)
    let genRow = genSVec w (1,w)
    m <- choose (1,h)
    rows <- genIntMap genRow h (1,h)
    return $ SM (h,w) rows

instance (Arbitrary α, Ord α, Eq α, Num α) => Arbitrary (SparseMatrix α) where
    arbitrary = genSparseMatrix =<< genPair genNonNegative genNonNegative
        
genAssocList :: (Arbitrary α, Ord α, Eq α, Num α) => Int -> (Int,Int) -> Gen [ ((Index,Index), α) ]
genAssocList n (h,w) = liftM ( (((h,w),0) :) . nubBy ((==) `on` fst) ) $
    vectorOf n $ genPair (genPair (genIndex (1,h)) (genIndex (1,w))) genNonZero


