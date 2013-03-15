module Tests.Vector where

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import Control.Monad (liftM)
import Data.Function (on)
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse.Vector

import Tests.Generators

testgr_vector = testGroup "VECTOR TESTS"
    [ testgr_to_from_list
    , testgr_muls
    ]

--------------------------------------------------------------------------------
-- LOOKUP/UPDATE --
-------------------

testgr_lookup_update = testGroup "LOOKUP/UPDATE"
    [ ]


--------------------------------------------------------------------------------
-- TO/FROM LIST --
------------------

testgr_to_from_list = testGroup "TO/FROM LIST"
    [ testProperty "fillVec . sparseList == id" prop_sparse_fill_vec
    , testProperty "sparseList . fillVec == id" prop_fill_sparse_vec
    , testProperty "vecToAssocList   . vecFromAssocList == id" prop_vec_to_from_AssocList
    , testProperty "vecFromAssocList . vecToAssocList   == id" prop_vec_from_to_AssocList
    ]


prop_fill_sparse_vec = 
    forAll (genSize (0,1000) >>= genSparseVector :: Gen (SparseVector Integer))
    $ \v -> sparseList (fillVec v) == v

prop_sparse_fill_vec = 
    forAll (genSize (0,1000) >>= genSparseList :: Gen [Integer]) 
    $ \l -> fillVec (sparseList l) == l

prop_vec_to_from_AssocList =
    forAll (genSize (0,1000)) $ \s ->
    forAll (choose (0,s)) $ \n ->
    forAll (genVecAssocList n s :: Gen [(Index,Integer)]) $ \l ->
        sort (vecToAssocList (vecFromAssocList l)) == sort l

prop_vec_from_to_AssocList =
    forAll (genSize (0,1000) >>= genSparseVector :: Gen (SparseVector Integer))
    $ \v -> vecFromAssocList (vecToAssocList v) == v

--------------------------------------------------------------------------------
-- DOT PRODUCT --
-----------------

testgr_muls = testGroup "DOT PRODUCT"
    [ testProperty "Vector dot-product" prop_dot
    ]


prop_dot = 
    forAll (genSize (1,1000)) $ \n ->
    forAll (genSparseVector n :: Gen (SparseVector Integer)) $ \x ->
    forAll (genSparseVector n :: Gen (SparseVector Integer)) $ \y ->
        x·y == fillVec x `dotL` fillVec y
    where
        dotL :: (Num α) => [α] -> [α] -> α
        dotL x y = L.sum $ zipWith (*) x y
