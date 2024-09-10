#!/bin/bash

set -xeuo pipefail

# https://github.com/guard/listen/blob/master/README.md#adapter-warnings
export LISTEN_GEM_ADAPTER_WARN_BEHAVIOR=silent

# Assume running on a Vertex AI Workbench VM which runs jupyter-server-proxy.
bundle exec jekyll serve --baseurl /proxy/absolute/4000

