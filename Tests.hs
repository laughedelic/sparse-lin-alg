module Main where

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import Tests.Vector     (testgr_vector)
import Tests.Matrix     (testgr_matrix)
import Tests.Algorithms (testgr_algorithms)

----------------------------------------------------------------------
main = defaultMain 
     [ testgr_vector
     , testgr_matrix
     , testgr_algorithms
     ]
