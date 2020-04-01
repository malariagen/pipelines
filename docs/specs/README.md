# Pipeline specifications

This folder contains pipeline specification documents. Each document
provides a written specification of a pipeline, intended as a basis
for implementation.

## Versioning policy

Each specification will include a three-part version number. The
version number will be incremented as follows:

* Micro version bump - Add or clarify information in spec, but no
  change to actual pipeline.

* Minor version bump - For changes to pipeline that should not
  qualitatively change outputs, but may in practice produce slightly
  different outputs. Also for changes that add new steps producing
  additional outputs, without qualitatively changing any previous
  outputs.

* Major version bump - For any changes that will qualitatively change
  outputs.

## Specifications

### Short read alignment (vector)

* [Latest version](short-read-alignment-vector.md)
