#!/bin/bash

./diff_timeline ne8  /gpfs/alpinetds/stf006/proj-shared/imn/HOMME_ACME/preqx.cpu/jw-ne8-nlev72-qsize40-thread1-nodes1-tasks84/movies/asp_baroclinic2.nc \
                     /gpfs/alpinetds/stf006/proj-shared/imn/HOMME_ACME/preqx.openacc/jw-ne8-nlev72-qsize40-thread1-nodes1-tasks6/movies/asp_baroclinic2.nc diff

cat diff_*.dat