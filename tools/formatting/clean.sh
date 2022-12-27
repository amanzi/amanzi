# to be run from AMANZI_SRC_DIR
cd src/physics/ats
python ../../../tools/formatting/clean.py --ats .
cd ../../
python ../tools/formatting/clean.py --exclude=physics .
. ../tools/formatting/clang-format.sh
