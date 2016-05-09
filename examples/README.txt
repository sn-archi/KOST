This directory contains some examples of how to use KOST:

----------------------
    accuracy
----------------------
Tests round-trip accuracy for kostStateVector2Elements and kostElements2StateVector, for a variety of different state vectors.

----------------------
    gnuplot
----------------------
This contains a small test application that uses GNUplot for displaying the shape of an orbit. It probably works best on UNIX-like systems, like Linux.

----------------------
    KOST_Orbiter
----------------------
This directory contains a copy of the KOST source code, set up for use in Orbiter. It should be identical to ../KOST, except for the following:
* All .c files are renamed to .cpp. This allows KOST to use Orbiter's data types, which are defined in C++-style header files (header files which also contain classes etc.)
* kost_settings.h is modified.

----------------------
    KOSTplanet
----------------------
An implementation of a planet module for Orbiter, based on KOST. The project is set up to compile when the entire KOST directory is placed in the Samples directory, so that this directory becomes Samples/KOST/examples/KOSTplanet.


