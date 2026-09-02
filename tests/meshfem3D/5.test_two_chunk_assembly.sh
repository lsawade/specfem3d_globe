#!/bin/bash

testdir=`pwd`
var=test_two_chunk_assembly

echo >> $testdir/results.log
echo "test: $var" >> $testdir/results.log
echo >> $testdir/results.log

mkdir -p bin

make -f $var.makefile $var >> $testdir/results.log 2>&1
if [ ! -e ./bin/$var ]; then
  echo "compilation of $var failed, please check..." >> $testdir/results.log
  exit 1
fi

error_log=`mktemp "${TMPDIR:-/tmp}/specfem-$var.XXXXXX"`
mpirun -np 4 ./bin/$var >> $testdir/results.log 2>$error_log
if [[ $? -ne 0 ]]; then
  echo "test failed"; cat $error_log; exit 1
fi
if [[ -s $error_log ]]; then
  echo "returned ERROR output:" >> $testdir/results.log
  cat $error_log >> $testdir/results.log
  exit 1
fi

echo "successfully tested: `date`" >> $testdir/results.log
