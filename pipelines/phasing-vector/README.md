# ShortReadAlignment (vector) pipeline

This folder contains the pipeline wdl and workflow inputs for the Vector phasing pipelines as specified in https://github.com/malariagen/pipelines/blob/797a6b1c152d7e262404575c15a6c36553af5c41/docs/specs/phasing-vector.md

There are gcp and farm5-specific directories which contain versions specific to gcp (Google Cloud Platform) and farm5 (Sanger LSF) backed instances of cromwell.  For the wdl (pipelines) the ONLY differences are in the import statements.  For the example inputs the paths in those are specific to the platform in use.

There are three pipelines implemented here.  ReadBackedPhasing.wdl, StatisticalPhasing.wdl and Phasing.wdl.  The latter is just a wrapper around the first two.
