#!/bin/bash

preload_size=16G

krakenuniq \
  --paired \
  --preload-size ${preload_size} \
  --threads ${threads} \
  --db ${db} \
  --report ${report} \
  --output off \
  ${R1} ${R2}
