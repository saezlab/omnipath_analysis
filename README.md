# The OmniPath2 paper analysis suite

This repository contains a Python module and an R package.
The purpose of both is to process and analyse the database contents of
OmniPath, export tables and create figures.
As the data always comes from the `pypath` Python module, the `omnipath2`
Python module in this repository is responsible for the extraction of the
data from all main database domains of OmniPath (network, enzyme-substrate,
complexes, annotations, inter-cellular communication).

The `omnipath2.database` module builds all databases according to the
parameters provided and saves them into `pickle` files. This way to load
all databases takes only 2-5 minutes.

Various modules extract information from the databases and generate tables
or create figures. For almost all plotting we use the `matplotlib` powered
base class in `omnipath2.plot`.

Some of the tables are supplementary tables, but we compile a number of
tables to provide input for the R package.

The `omnipath2` R package in the `R` directory reads these tables, further
processes the data and creates further figures.
