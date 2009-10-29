#!/bin/bash

[ -h "genericloss.cpp" ] || ln -s ../bmrm-2.1/loss/genericloss.cpp .
[ -h "genericloss.hpp" ] || ln -s ../bmrm-2.1/loss/genericloss.hpp .
[ -h "genericdata.cpp" ] || ln -s ../bmrm-2.1/data/genericdata.cpp .
[ -h "genericdata.hpp" ] || ln -s ../bmrm-2.1/data/genericdata.hpp .

cd ../bmrm-2.1/linear-bmrm/

make && cp linear-bmrm-predict ../../sourcefiles/ && cp linear-bmrm-train ../../sourcefiles/
