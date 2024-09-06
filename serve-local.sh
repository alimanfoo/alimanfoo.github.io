#!/bin/bash

set -xeuo pipefail

# https://github.com/guard/listen/blob/master/README.md#adapter-warnings
export LISTEN_GEM_ADAPTER_WARN_BEHAVIOR=silent

bundle exec jekyll serve
