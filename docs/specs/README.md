# Pipeline specifications

This folder contains pipeline specification documents. Each document
provides a written specification of a pipeline, intended as a basis
for implementation.


## Specification documents

### Mosquito short read alignment pipeline

* [Version 1.2.1](https://github.com/malariagen/pipelines/blob/c1531ef6120021106cb0159150a297a5d8473e07/docs/specs/short-read-alignment-vector.md)
* [Version 1.1.1](https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md)

### Mosquito SNP genotyping pipeline

* [Version 1.4.1](https://github.com/malariagen/pipelines/blob/a910bde50e824623a709d86465658d9db44e1ccd/docs/specs/snp-genotyping-vector.md)
* [Version 1.4.0](https://github.com/malariagen/pipelines/blob/c1531ef6120021106cb0159150a297a5d8473e07/docs/specs/snp-genotyping-vector.md)

### Mosquito phasing pipeline

* [Version 0.0.0](https://github.com/malariagen/pipelines/blob/9364efde80c6b2745290bba120b0f0a07a98def2/docs/specs/phasing-vector.md)


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
