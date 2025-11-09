#!/bin/bash

# Check arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 <varname>" 1>&2
  exit 2
fi

# Save arguments
varname=$1

# Check correct folder
if [ "`basename $(pwd)`" != "src" ]; then
  echo "error: must execute the script from 'src/' directory of ZFS." 1>&2
  exit 1
fi

# Check correct compiler/build type
compiler=$(cd ..; make what 2>&1 | grep "Compiler:" | awk '{print $2}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g")
build_type=$(cd ..; make what 2>&1 | grep "Build type:" | awk '{print $3}' | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g")
if [ "$compiler" != "GNU" ] || [ "$build_type" != "debug" ]; then
  echo "error: compiler/build type must be 'GNU/debug'." 1>&2
  exit 1
fi

# Set list of make targets
targets="zfscartesiangrid_hilbert.o zfscartesiangrid_inst_avg.o zfscartesiangrid_inst_fv.o zfscartesiangrid_inst_lb.o zfscartesiangrid_inst_dg.o zfsgridgenpar.o"

# Make clean
make clean &>/dev/null

# Iterate over make targets
for target in $targets; do
  # Info to user
  echo "Next up: $target"

  # Reset iteration number
  iteration=1

  while :; do
    # Run make
    printf "  Running make (iteration: $iteration)... "
    errorlog=$(make $target 2>&1)
    returncode=$?
    echo "done"

    if [ $returncode -eq 0 ]; then
      echo "  Make finished suceessfully, moving on..."
      break
    fi

    # Set additional lines
    ((additional_lines=iteration-1))

    # Get shadowing occurrences
    occurrences=$(echo "$errorlog" | grep "declaration of .${varname}. shadows a member of" | awk -F: '{print $1 " " $2}' | sort | uniq)
    no_occ=$(echo "$occurrences" | wc -l)
    if [ -z "$occurrences" ]; then
      no_occ=0
    fi
    echo "  Found $no_occ place(s) where '$varname' was shadowed."

    if [ "$no_occ" != "0" ]; then
      # Replace shadowing occurrences
      printf "  Replacing shadowed variables."
      echo "$occurrences" | while read occurrence; do
        file=$(echo $occurrence | awk '{print $1}')
        line=$(echo $occurrence | awk '{print $2}')
        from=$line
        ((to=line+additional_lines))
        sed -i "$from,$to s/\([^.>]\)\<${varname}\>\([^(]\)/\1${varname}_\2/g" $file
        printf "."
      done
      echo
    fi

    # Get other occurrences
    occurrences=$(echo "$errorlog" | grep "^`dirname $(pwd)`" | grep "error\|note" | awk -F: '{print $1 " " $2}' | sort | uniq)
    no_occ=$(echo "$occurrences" | wc -l)
    if [ -z "$occurrences" ]; then
      no_occ=0
    fi
    echo "  Found $no_occ place(s) where '$varname' was otherwise used."

    # Replace occurrences
    if [ "$no_occ" != "0" ]; then
      printf "  Replacing other occurrences of '${varname}'."
      echo "$occurrences" | while read occurrence; do
        file=$(echo $occurrence | awk '{print $1}')
        line=$(echo $occurrence | awk '{print $2}')
        from=$line
        ((to=line+additional_lines))
        sed -i "$from,$to s/\([^.>]\)\<${varname}\>\([^(]\)/\1${varname}_\2/g" $file
        printf "."
      done
      echo
    fi

    # Increment counter
    ((iteration++))
  done
  echo
done

echo 'Finished!'
