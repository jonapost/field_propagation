#/bin/csh -v
#mkdir -p $g4version
taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar.mac >& em2_standard.out
taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar_emv.mac >& em2_emv.out
taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar_emx.mac >& em2_emx.out

taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar_20um.mac >& em3_standard.out
taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar_emv_20um.mac >& em3_emv.out
taskset -c 0 $G4BIN/Linux-g++/be2 $BENCHMARKS/electromagnetic/be2/atlasbar_emx_20um.mac >& em3_emx.out
#
