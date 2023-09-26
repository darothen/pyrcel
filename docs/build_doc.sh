#!/usr/bin/env bash
# Quick script to build the documentation

# Generate .rst files from the notebooks containing pre-rendered examples
cd examples
make

# Go back and generate the stub automethod templates in the API reference
cd ../
make gen

# Now make the html documents
make html
