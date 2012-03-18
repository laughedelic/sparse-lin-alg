module Main where

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import Control.Monad (liftM)
import SparseMatrix

import Data.Function (on)
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))


----------------------------------------------------------------------
main = defaultMain tests

tests = [ testgr_lookup_update
        , testgr_to_from_list
        , testgr_trans
        , testgr_muls
        ]
----------------------------------------------------------------------


-- | Generates non negative number
genNonNegative = liftM abs arbitrary

genNonZero :: (Arbitrary α, Num α) => Gen α
genNonZero = do x <- arbitrary
                return $ if x == 0 then 1 else x

-- | Combains two generators to pair
genPair g1 g2 = do
    x <- g1
    y <- g2
    return (x,y)

genSize (m,n) = choose (m,n)

genIndex (m,n) | m > n     = return 0
               | otherwise = choose (m,n)

--------------------------------------------------------------------------------
-- SPARSE VECTOR DATATYPE --
----------------------------

-- | Generates IntMap of given size with indecies in given diapason 
--   and arbitrary values
genSVec :: (Arbitrary α, Num α) => Int -> (Int,Int) -> Gen (SVec α)
genSVec = genSVec' genNonZero

-- | Generates sparse vector
genSVec' :: Gen α       -- ^ generator for values
         -> Int         -- ^ size of IntMap
         -> (Int,Int)   -- ^ diapason of indicies
         -> Gen (SVec α) 
genSVec' g n (a,b) = liftM M.fromList 
                   $ vectorOf n $ genPair (choose (a,b)) g

-- | Generates `SparseVector` of given length with arbitrary values
genSparseVector :: (Arbitrary α, Num α) => Int -> Gen (SparseVector α)
genSparseVector m = do
    n <- choose (0,m)
    v <- genSVec n (1,m)
    return $ SV m v

-- | Arbitrary `SparseVector` of random length
instance (Arbitrary α, Num α) => Arbitrary (SparseVector α) where
    arbitrary = genSparseVector =<< genNonNegative

-- | Generates list of given length with frequent zeroes
genSparseList :: (Arbitrary α, Num α) => Int -> Gen [α]
genSparseList n = vectorOf n $ frequency [(4, return 0), (1, arbitrary)]

--------------------------------------------------------------------------------
-- SPARSE MATRIX DATATYPE --
----------------------------

genSparseMatrix :: (Arbitrary α, Num α) => (Int,Int) -> Gen (SparseMatrix α)
genSparseMatrix (h,w) = do
    n <- choose (1,w)
    let genRow = genSVec n (1,w)
    m <- choose (1,h)
    rows <- genSVec' genRow m (1,h)
    return $ SM (h,w) rows

instance (Arbitrary α, Num α) => Arbitrary (SparseMatrix α) where
    arbitrary = genSparseMatrix =<< genPair genNonNegative genNonNegative
        
genAssocList :: (Arbitrary α, Num α) => Int -> (Int,Int) -> Gen [ ((Index,Index), α) ]
genAssocList n (h,w) = liftM ( (((h,w),0) :) . nubBy ((==) `on` fst) ) $
    vectorOf n $ genPair (genPair (genIndex (1,h)) (genIndex (1,w))) genNonZero

--------------------------------------------------------------------------------
-- LOOKUP/UPDATE --
-------------------

testgr_lookup_update = testGroup "LOOKUP/UPDATE"
    [ testProperty "Insert/Erase" prop_ins ]


-- Test for `ins`, `(#)` and `erase`
prop_ins = 
    forAll (genPair (genSize (0,300)) (genSize (1,400))) $ \(h,w) ->
    forAll (genSparseMatrix (h,w))                       $ \m ->
    forAll (genPair (genIndex (1,h))   (genIndex (1,w))) $ \(i,j) ->
    forAll (arbitrary :: Gen Int)                        $ \x -> 
    -- collect (h,w) $
        let m' = m `ins` ((i,j),x)
        in  m' # (i,j) == x
        && (m' `erase` (i,j)) # (i,j) == 0

--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

testgr_to_from_list = testGroup "TO/FROM LIST"
    [ testgr_to_from_list_vectors
    , testgr_to_from_list_matrices
    ]


-----------
-- Vectors:

testgr_to_from_list_vectors = testGroup "Vectors"
    [ testProperty "fillVec . sparseList == id" prop_sparse_fill_vec
    , testProperty "sparseList . fillVec == id" prop_fill_sparse_vec
    ]


prop_fill_sparse_vec = 
    forAll (genSize (0,1000) >>= genSparseVector :: Gen (SparseVector Integer))
    $ \v -> sparseList (fillVec v) == v

prop_sparse_fill_vec = 
    forAll (genSize (0,1000) >>= genSparseList :: Gen [Integer]) 
    $ \l -> fillVec (sparseList l) == l

------------
-- Matrices:

testgr_to_from_list_matrices = testGroup "Matrices"
    [ testProperty "fillMx   . sparseMx == id" prop_fill_sparse_mx
    , testProperty "sparseMx . fillMx   == id" prop_sparse_fill_mx
    , testProperty "toAssocList   . fromAssocList == id" prop_to_from_AssocList
    , testProperty "fromAssocList . toAssocList   == id" prop_to_from_AssocList
    ]


prop_sparse_fill_mx =
    forAll (genPair (genSize (1,100)) (genSize (1,150))) $ \(h,w) ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        sparseMx (fillMx m) == m

prop_fill_sparse_mx = 
    forAll (genPair (genSize (1,100)) (genSize (1,150))) $ \(h,w) ->
    forAll (vectorOf h (genSparseList w) :: Gen [[Integer]]) $ \m ->
        fillMx (sparseMx m) == m

prop_to_from_AssocList =
    forAll (genPair (genSize (1,50)) (genSize (1,60))) $ \(h,w) ->
    forAll (choose (0,h*w)) $ \n ->
    forAll (genAssocList n (h,w) :: Gen [((Index,Index),Integer)]) $ \l ->
        sort (toAssocList (fromAssocList l)) == sort l

prop_from_to_AssocList =
    forAll (genPair (genSize (1,100)) (genSize (1,150))) $ \(h,w) ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        -- collect (toAssocList m) $
        fromAssocList (toAssocList m) == m

--------------------------------------------------------------------------------
-- TRANSPOSITION --
-------------------

testgr_trans = testGroup "TRANSPOSITION"
    [ testProperty "Transposition" prop_trans ]


prop_trans =
    forAll (genPair (genSize (1,100)) (genSize (1,150))) $ \(h,w) ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        trans (trans m) == m

--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

testgr_muls = testGroup "MULTIPLICATIONS"
    [ testProperty "Vector dot-production" prop_dot
    , testProperty "Matrix-by-vector multiplication" prop_mulMV
    , testProperty "Vector-by-matrix multiplication" prop_mulVM
    , testProperty "Sparse matrices multiplication"  prop_mul
    ]


prop_dot = 
    forAll (genSize (1,1000)) $ \n ->
    forAll (genSparseVector n :: Gen (SparseVector Integer)) $ \x ->
    forAll (genSparseVector n :: Gen (SparseVector Integer)) $ \y ->
        x·y == fillVec x `dot'` fillVec y

dot' :: (Num α) => [α] -> [α] -> α
dot' x y = L.sum $ zipWith (*) x y

prop_mulMV = 
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,200)) $ \w ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
    forAll (genSparseVector  h    :: Gen (SparseVector Integer)) $ \v ->
        fillVec (m ×· v) == [ row `dot'` fillVec v | row <- fillMx m ]

prop_mulVM = 
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,200)) $ \w ->
    forAll (genSparseVector  h    :: Gen (SparseVector Integer)) $ \v ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        fillVec (v ·× m) == [ fillVec v `dot'` col | col <- transpose (fillMx m) ]

mm :: (Num α) => [[α]] -> [[α]] -> [[α]]
a `mm` b = [ [ row `dot'` col | col <- transpose b ] | row <- a ]

prop_mul = 
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,150)) $ \n ->
    forAll (genSize (1,130)) $ \w ->
    forAll (genSparseMatrix (h,n) :: Gen (SparseMatrix Integer)) $ \a ->
    forAll (genSparseMatrix (n,w) :: Gen (SparseMatrix Integer)) $ \b ->
        fillMx  a `mm` fillMx  b == fillMx (a × b)
