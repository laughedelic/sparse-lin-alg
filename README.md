# Synopsis #

This is a small Haskell library providing effective representation for sparse matrices and some linear algebra operations.

# Installation #

1. Dowload it:

    $ git clone git://github.com/laughedelic/sparse-matrix.git

2. Configure
   - There are QuickCheck test, so if additional depensies are ok for you, configure with tests:

    $ cabal configure --enable-tests

   - Else, the only build-depensy is `base`, so just configure

    $ cabal configure

3. Build package:

    $ cabal build

4. Build Haddock documentation (optional):

    $ cabal haddock

5. Run tests, if you configured package with `--enable-tests` option

    $ ./dist/build/tests/tests

  You can use some options with this test-script (for example run tests concurrently), see option `-h`.

6. Install package:

    $ cabal install
