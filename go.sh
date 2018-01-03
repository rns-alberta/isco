#!/bin/bash

eos="APR"
e_pow="e14"
ee="8"
emax="22"

make isco

isco -f "../eos/eos$eos" -e "$ee$e_pow" -l "$emax$e_pow" -n 20 -s 300

times
