#!/bin/sh

#author: Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>

NC="\033[0m"
RED="\033[1;31m"
GREEN="\033[32m"

# Default values
cmd_fix=""
cmd_j="-j 12"

# Parse options
while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo -e "This scripts runs clang-tidy-diff in parallel and provides a summary of the checked files.\n"
      echo -e "Available options:\n"
      echo -e "  -h, --help      Print this help"
      echo -e "  -j, --num-jobs  Number of jobs used (Default: 12)"
      echo -e "  --fix           Apply fixes (if available)"
      exit 0 
      ;;
    -j|--num-jobs)
      cmd_j="-j $2"
      shift # past argument
      shift # past value 
      ;;
    --fix)
      cmd_fix="-fix"
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

# Check if compile commands are available
res=`ls | grep "compile_commands.json"`
if [ -z $res ]
then
  echo -e "Missing compile_commands.json in current directory. Generate via 'configure.py 1 2 --cc' before running clang-tidy."
  exit 1 
fi

# Check if clang tidy config is available
res=`ls -a | grep ".clang-tidy"`
if [ -z $res ]
then
  echo -e "Missing .clang-tidy in current directory. Please use script in root directory of project."
  exit 1 
fi

# Execure clang-tidy-diff in parallel
forkPoint=`git merge-base --fork-point remotes/origin/master`
out=`git diff $forkPoint -U0 | /pds/opt/llvm/share/clang/clang-tidy-diff.py -clang-tidy-binary=/pds/opt/llvm/bin/clang-tidy -p1 $cmd_j $cmd_fix 2>/dev/null`

# Check for compilation errors
if [[ "$out" == *"error:"* ]]
then
  echo -e "$out"
  echo -e "Compilation error found! Fix the compilation before running clang-tidy."
  exit 1 
fi

# Parse output
exitCode=0
summary="--clangFormat: summary------------------\n"
for f in $(git diff $forkPoint --name-only | grep -e '\.h' -e '\.cpp')
do
  summary="$summary$f:"
  res=`echo "$out" | grep $f`
  if [ -z "$res" ]
  then
    summary="$summary ${GREEN}PASSED${NC}\n"
  else
    summary="$summary ${RED}FAILED${NC}\n"
    exitCode=1
  fi
done
summary="$summary----------------------------------------"

# Print output
if [ ! -z "$out" ]
then
  echo -e "$out"
fi

# Summary of suggest --fix
if [ -z $cmd_fix ]
then
  echo -e "$summary"
else
  echo "Please rerun this script ($0) without --fix to check if all issues have been resolved."
  echo "Note: clang-format should be checked after automatic fixes have been applied."
  exit 0;
fi

if [ $exitCode == 1 ]
then
  echo "Please rerun this script ($0) with --fix or fix by hand."
fi

exit $exitCode
