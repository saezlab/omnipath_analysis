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

## Workflow

The `main` module contains a description of the workflow and a class to call
all its elements.

## Database

The class `database.Database` has the methods to compile the data in order
to make sure we use the databases always the same way. Also it saves
each database to pickles and serves them to all the other modules. It
is instantiated in `__init__.py` hence you can do for example to get the
omnipath network:
    
```
import omnipath2
op = omnipath2.data.get_db('omnipath')
```

The pickle dumps I use you can find here (to be extracted to the `pickles`
directory):

http://denes.omnipathdb.org/54b510889336eb2591d8beff/omnipath2_pickles.tar.gz
(116M)

## Networks

At the moment I created 2 PPI networks, one is called `omnipath` and
corresponds more or less to the data in the web service, it contains
all the `extra` datasets; and the other one called `curated` and
contains only the literature curated data. In addition we have the
TF-target, miRNA-mRNA and TF-miRNA networks.

## Extracting data from the databases

A number of methods for combining the network with the intercell
annotations are already in the `pypath.annot.CustomAnnotation` class,
hence we don't need to implement these in `omnipath2`. Also for each
database domains lots of methods to query various numbers and subsets
are in the relevant class in `pypath`, so we only need to call these.
Just 2 examples:
`pypath.intercell.IntercellAnnotation.degree_inter_class_network_inhibitory`
will provide the degree distribution of inhibitory edges between 2
intercell classes, or
`pypath.main.PyPath.interactions_mutual_by_resource` will return the
mutual interactions grouped by resource, and so on.

## Plotting in R

The `r_preprocess` module  exports tables for the R plotting scripts
-- By default figures and tables are exported to the `figures` and
`tables` directories into a subdirectory with the current date e.g.
`tables/20191028/intercell_classes_20191028.tsv`, also each file name
by default contains the date. In the `R` directory there is an R package
called `omnipath2`; this works from the tables exported by the Python module.

## Supplementary tables

The `supptables` module compiles and exports the supplementary tables S2-S6
with lots of numbers about each database domain.

## Parameters

The file names are in the `omnipath2.settings` module and might
contain variable fields as for certain figures we have more different
versions, e.g. compiled from different networks

## Color palettes

The `omnipath2.colors` module reads the palettes and provides them
to the other modules.

## Plotting in Python

The `omnipath2.plot` module contains a base class `PlotBase` which we use
for most of the plots.
The `omnipath2.intercell_plots`, `omnipath2.annotation_plots`,
`omnipath2.network_plots`, `omnipath2.complexes_plots` modules have the
plots created by Nico and wrapped into classes inheriting from
`omnipath2.plot.PlotBase`.

## Directory structure

The directories can be customized in the `settings` module.
One directory, by default `pickles`, contains the pickle dumps of the
databases.
The `figures` and `tables` directories contain the figures and tables
generated by the module. Both the Python and the R methods use the
same directories. Optionally the figures and tables can be collected into
timestamped subdirectories which helps to keep track of their development,
check or remove the old ones or find the most recent ones. Also optionally
a timestamp can be added to each file name. The timestamps by default are
an 8 digit representation of the date, so every day a new directory will
be started (except if you keep running the same session over midnight, then
still the old directory will be used).
The `files.json` file keeps track of the most recent versions of all tables
and figures generated both for the Python and R part as well as the earlier
versions.

## Resource usage

With all databases loaded Python requires maybe around 5G memory.
The inter-cellular network data frame requires little more than 1G.
`pandas` operations can result 1-3G peaks on top of that baseline.
Overall it's good to have 8G RAM available.

## How to run

You can run the whole workflow by calling `main.py`.

```
python main.py
```

You can run selected parts or tasks only:

```
from omnipath2 import main

# 2 parts selected (each a series of tasks)
workflow = main.Main(parts = ['intercell_plots', 'network_plots'])
workflow.main()

# only selected tasks:
workflow = main.Main(steps = main.workflow['intercell_plots'][0])
workflow.main()
```

The R workflow is called by the `r_plotting` part of the Python one, but
also can be run separately:

```
Rscript omnipath2_workflow.r
```
