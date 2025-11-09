#!/bin/bash

# Checks for unused testcase properties (not used by any rank)
# 1. Run testcases with zfs compiled with ZFS_WRITE_ACCESS_PROPERTIES_FILE defined in zfsconfig.h
# 2. Call this script from e.g. testcases/

baseDir=$PWD

# Loop over all testcases found in current directory
for tc in `find . -maxdepth 3 -mindepth 0 -type d`; do
  cd $tc
  isTc=$(ls -l | grep properties.*toml | wc -l)
  if [ $isTc -eq 0 ]; then
    cd $baseDir
    continue
  fi

  # number of property access files in the testcase
  noPropFiles=$(ls -l | grep access_properties_domain0 | wc -l)
  echo "Testcase: $tc, noAccessPropertyFiles=$noPropFiles"

  if [ $noPropFiles -eq 0 ]; then
    echo ""
    cd $baseDir
    continue
  fi

  # loop over all property access files (for grid generation, run, restart, ...)
  for i in `seq 0 $(($noPropFiles-1))`; do
    logfile=unused_props_${i}.log

    # number of ranks
    noRanks=$(ls -l access_properties_domain*_$i | wc -l)

    # find properties that are unused on all ranks (Accesses 0 appears noRanks times)
    grep -h "Accesses 0 " access_properties_domain*_$i | sort | uniq -c | sort -n -k 1 | grep ^[[:space:]]*$noRanks > $logfile

    noUnusedProps1=$(cat $logfile | wc -l)

    # Remove testcase properties not read by zfs but used by canary
    sed -i '/noDomains/d' $logfile
    sed -i '/noGridGenProcs/d' $logfile
    sed -i '/omp_num_threads/d' $logfile
    # Properties that are not used at a restart
    sed -i '/gridInputFileName/d' $logfile

    noUnusedProps2=$(cat $logfile | wc -l)

    echo "Access properties #$i: $noUnusedProps2 unused properties (with excluded properties, e.g. noDomains: $noUnusedProps1)"


    # Prepare script to delete unused properties
    # TODO does not work for multiblock properties with .default/.X values
    # TODO does only delete line with the property name in it
    awk '{print $3}' $logfile > unused_props.tmp
    mapfile -t plist < ./unused_props.tmp
    pfile='$1'
    deletescript="remove_unused_props_$i.sh"
    > $deletescript
    echo "# pass property*.toml file as argument" > $deletescript

    for i in ${plist[@]}
    do
      echo "echo 'deleting property " $i "'" >> $deletescript
      echo "sed -i '/\<"$i"\>/d' $pfile" >> $deletescript
      # echo "sed -i '/"$i" /d' $pfile" >> $deletescript
      # echo "sed -i '/"$i"=/d' $pfile" >> $deletescript
    done

    # this removes double blank lines
    echo "awk 'NF > 0 {blank=0} NF == 0 {blank++} blank < 2' $pfile > tmpfile.tmp" >> $deletescript
    echo "cat tmpfile.tmp > $pfile" >> $deletescript

    chmod u+x $deletescript

    rm unused_props.tmp

  done
  echo ""
  cd $baseDir
done
