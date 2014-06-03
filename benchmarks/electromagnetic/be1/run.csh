#/bin/csh -v
#mkdir -p $g4version
taskset -c 0  $G4BIN/Linux-g++/be1 $BENCHMARKS/electromagnetic/be1/cms10gev.mac >& em1_standard.out
taskset -c 0  $G4BIN/Linux-g++/be1 $BENCHMARKS/electromagnetic/be1/cms10gev_emv.mac >& em1_emv.out
taskset -c 0  $G4BIN/Linux-g++/be1 $BENCHMARKS/electromagnetic/be1/cms10gev_emx.mac >& em1_emx.out
#
