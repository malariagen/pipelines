#!/usr/bin/env bash

detect_missing_samples() {
  local manifest=$1
  shift
  echo Diffing ${manifest}. Bucket vs sample count: $(s3cmd ls s3://$manifest | grep -Ev '^$|provenance' | sed 's:.*s3\://\(.*\)/\([^.]*\)\..*:\2:g' | sort -u | wc -l) vs $(tail -n +2 $manifest.tsv | grep -v '^$' | awk '{print $1}' | sort -u | wc -l)
  diff <(cat $manifest.tsv | tail -n +2 | grep -v '^$' | awk '{print $1}' | sort -u) <(s3cmd ls s3://$manifest| grep -Ev '^$|provenance' | sed 's:.*s3\://\(.*\)/\([^.]*\)\..*:\2:g' | sort -u)
}

for bucket in $@; do
  detect_missing_samples "${bucket}"	
  echo
done
