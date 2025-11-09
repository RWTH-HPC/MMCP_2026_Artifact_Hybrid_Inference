#!/bin/bash

if [ $# -eq 1 ]; then
  filename="$1"
else
  filename=coupling_workshop_overview
fi

pandoc ${filename}.rst -o ${filename}.pdf
