# Running a pipeline

This document describes basic steps needed to run a pipeline.

## Building a release.

Since a pipeline wdl may use other workflows, tasks, and structs that are not stored in 
the same file, we need to build a dependency zip for submission to cromwell along
with the main pipeline wdl.  To package the wdl, dependency zip and options file into a single
release, there is a script in the scripts directory named 'build_pipeline_release'.

Note that when this script builds the dependency zip file it *only* includes external files that are referenced in the main pipeline wdl.

An example invocation is:

`./scripts/build_pipeline_release.sh -w pipelines/short-read-alignment-vector/farm5/ShortReadAlignment.wdl -v 0.1.0 -o releases`

This builds a directory named 'releases' which contains the pipeline wdl, dependency zip and options file

To execute this wdl (using cromshell):

`cromshell submit -w releases/ShortReadAlignment/ShortReadAlignment_0.1.0.wdl pipelines/short-read-alignment-vector/farm5/input_files/small/AV0079-C.json releases/ShortReadAlignment/ShortReadAlignment_0.1.0.options.json releases/ShortReadAlignment/ShortReadAlignment_0.1.0.zip`

NB:  If you are using cromshell > v0.3.9, there is a bug where cromshell tries to invoke womtool to validate the WDL but doesn't pass in the dependencies zip.  See https://github.com/broadinstitute/cromshell/issues/139.  To workaround this, ensure that womtool is not in your PATH, or use cromshell v0.3.9.
