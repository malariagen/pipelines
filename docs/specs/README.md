# Pipeline specifications

This folder contains pipeline specification documents. Each document
provides a written specification of a pipeline, intended as a basis
for implementation.


## Specification versioning policy

Each specification document will include a three-part version
number. The version number will be incremented as follows:

* Major version bump - For any changes that will qualitatively change
  outputs.

* Minor version bump - For changes to pipeline that should not
  qualitatively change outputs, but may in practice produce slightly
  different outputs. Also for changes that add new steps producing
  additional outputs, without qualitatively changing any previous
  outputs.

* Micro version bump - Add or clarify information in spec, but no
  change to actual pipeline.


## Specification documents

### Short read alignment (vector)

* [Latest version](short-read-alignment-vector.md)
