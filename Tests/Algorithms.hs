module Tests.Algorithms where

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import Control.Monad (liftM)
import Data.Function (on)
import Data.Foldable as F
import Data.List     as L
import Data.IntMap   as M hiding ((!))

import Math.LinearAlgebra.Sparse

import Tests.Matrix
import Tests.Generators

testgr_algorithms = testGroup "ALGORITHMS TESTS"
    [ testProperty "Staircase form" prop_staircase
    , testProperty "Diagonal form"  prop_diagonal
    --, testProperty "Linear system solution" prop_solve
    ]

--------------------------------------------------------------------------------
-- ALGORITHMS --
----------------

prop_staircase =
    forAll (genSize (1,50)) $ \h ->
    forAll (genSize (1,70)) $ \w ->
    forAll (genSparseMatrix (h,w) :: Gen (SparseMatrix Integer)) $ \m ->
        let (s,t) = staircase m
            -- (sz,en) = (size (mx m), L.sum [ size r | r <- elems (mx m) ])
        in  -- collect (sz,en,show (en `div` 100)++"%") $ -- (h,w) $ 
            t × m == s

prop_diagonal =
    forAll (genSparseMatrix (30,50) :: Gen (SparseMatrix Integer)) $ \m ->
        let (d,t,u) = toDiag m
            (sz,en) = (size (mx m), L.sum [ size r | r <- elems (mx m) ])
        in  -- collect (sz,en,show (en `div` 100)++"%") $ -- (h,w) $ 
            t × m × (trans u) == d

prop_solve =
    forAll (genSparseMatrix (10,5) :: Gen (SparseMatrix Integer)) $ \m ->
    forAll (genSparseVector  10     :: Gen (SparseVector Integer)) $ \b ->
        let x = solveLinear m b
        in collect (m,b,x) $
           case x of
                Nothing -> True
                Just x' -> m ×· x' == b
    --forAll (genSparseVector  5    :: Gen (SparseVector Integer)) $ \x ->
    --    let b = m ×· x
    --        x' = solveLinear m b
    --    in collect (m,x,b) $
    --        Just x == solveLinear m b
