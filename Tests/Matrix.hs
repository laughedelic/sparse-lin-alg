module Tests.Matrix where

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import Control.Monad (liftM)
import Data.Function (on)
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Matrix
import Math.LinearAlgebra.Sparse.Vector

import Tests.Generators

testgr_matrix = testGroup "MATRIX TESTS"
    [ testgr_lookup_update
    , testgr_to_from_list
    , testgr_trans
    , testgr_muls
    ]

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
    [ testProperty "fillMx   . sparseMx == id" prop_fill_sparse_mx
    , testProperty "sparseMx . fillMx   == id" prop_sparse_fill_mx
    , testProperty "toAssocList   . fromAssocList == id" prop_to_from_AssocList
    , testProperty "fromAssocList . toAssocList   == id" prop_from_to_AssocList
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
        collect (size (mx m), L.sum [ size r | r <- elems (mx m) ]) $
        trans (trans m) == m

--------------------------------------------------------------------------------
-- MULTIPLICATIONS --
---------------------

testgr_muls = testGroup "MULTIPLICATIONS"
    [ testProperty "Matrix-by-vector multiplication" prop_mulMV
    , testProperty "Vector-by-matrix multiplication" prop_mulVM
    , testProperty "Sparse matrices multiplication"  prop_mul
    ]

dotL :: (Num α) => [α] -> [α] -> α
dotL x y = L.sum $ zipWith (*) x y

prop_mulMV = 
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,200)) $ \w ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
    forAll (genSparseVector  h    :: Gen (SparseVector Integer)) $ \v ->
        fillVec (m ×· v) == [ row `dotL` fillVec v | row <- fillMx m ]

prop_mulVM = 
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,200)) $ \w ->
    forAll (genSparseVector  h    :: Gen (SparseVector Integer)) $ \v ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        fillVec (v ·× m) == [ fillVec v `dotL` col | col <- transpose (fillMx m) ]

mm :: (Num α) => [[α]] -> [[α]] -> [[α]]
a `mm` b = [ [ row `dotL` col | col <- transpose b ] | row <- a ]

prop_mul =
    forAll (genSize (1,100)) $ \h ->
    forAll (genSize (1,150)) $ \n ->
    forAll (genSize (1,130)) $ \w ->
    forAll (genSparseMatrix (h,n) :: Gen (SparseMatrix Integer)) $ \a ->
    forAll (genSparseMatrix (n,w) :: Gen (SparseMatrix Integer)) $ \b ->
        fillMx  a `mm` fillMx  b == fillMx (a × b)
