#!/bin/bash
# A simple script to iterate over all wdl files in the pipelines directory and run 'miniwdl check' on them.
# Do not exit on the first failure.

mainExitCode=0
for file in $(find pipelines/ -name '*.wdl' -type f -print)
do
  miniwdl check $file
  exitCode=$?
  if [ $exitCode != 0 ]; then
    mainExitCode=$exitCode
  fi
done
exit $mainExitCode