#!/bin/sh

################################################################################
# Settings
################################################################################

# Old name, new name, and new name for G cells
OLDNAME="m_variables"
NEWNAME="a_pvariable"
NEWNAMEG="${NEWNAME}G"

# Number of variables in [] to capture
NO_PATTERNS=2

# Set whether this variable should have grid() scope (1 - yes, 0 - no)
GRID_SCOPE=0


################################################################################
# Prepare file
################################################################################

# File to be changed
FILE=$1
BASE="$(basename "${FILE}")"

# Temporary newline replacement (must be present in DEPTH0 pattern below!)
NEWLINE="@@"

# Return early if search string (old name) is not even found
grep -q "$OLDNAME" $FILE || exit 0

# Make sure that names are only matched until end of the OLDNAME search string
OLDNAME="$OLDNAME\b"

# Read in file and change newlines to $NEWLINE
CONTENT=$(cat $FILE | sed ':a; N; $! b a; s/\n/'"$NEWLINE"'/g')


################################################################################
# Basic regex elements
################################################################################

# Save tab for convenience
TAB="	"

# Patterns
DEPTH0="[0-9a-zA-Z>*.,+%\/_() -:@${TAB}]*"
DEPTH1="\[\(${DEPTH0}\)*\]"
DEPTH2="\[\(${DEPTH0}\|${DEPTH1}\)*\]"
DEPTH3="\[\(${DEPTH0}\|${DEPTH2}\)*\]"
PATTERN="\[\(\(${DEPTH0}\|${DEPTH3}\)*\)\]"

# If variable should be a grid variable, set scope accordingly
if [ $GRID_SCOPE -eq 1 ]; then
  GRID="grid()."
else
  GRID=""
fi

# For function without explicit cell collector
# Pattern Reference is \2 for cellId
# Pattern Reference is \(2+(depth+1))+1 for position 1 == \7
# Pattern Reference is \(7+(depth+1))+1 for position 2 == \12
# ...

# For function with explicit cell collector
# Pattern Reference is \3 for cellId
# Pattern Reference is \(3+(depth+1))+1 for position 1 == \8
# Pattern Reference is \(8+(depth+1))+1 for position 2 == \13
# ...

# Variable to search for, uses as much patterns as there are positions to catch
# for example:
# SUFFIX="\.m_parentId"
# SUFFIX="\.m_neighborIds${PATTERN}${PATTERN}"

SUFFIX="[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*"
for i in `seq 1 $NO_PATTERNS`; do
  SUFFIX="${SUFFIX}${PATTERN}"
done

# Prefix
# Group 1 catches m_cells->a cells m_pCells
PREF1_1="\<m_cells->a\>"
PREF1_2="\<cells\>"
PREF1="\(${PREF1_1}\|${PREF1_2}\)"

# Group 2 catches input_cells->a cpu_cells->a
PREF2_1="\<input_cells\>"
PREF2_2="\<cpu_cells\>"
PREF2="\(\(\<${PREF2_1}\|\<${PREF2_2}\|\<${PREF2_3}\)->a\>\)"

# Group 3 catches cellsInput.a
PREF3_1="\<cellsInput\>"
PREF3="\(\(\<${PREF3_1}\).a\>\)"

# Group 4 catches gCells
PREF4_1="\<gCells\>"
PREF4_2="\<m_gCells->a\>"
PREF4="\(${PREF4_1}\|${PREF4_2}\)"


################################################################################
# Helper functions
################################################################################

# Take a list of filetypes and return true iff file belongs to all types
ftype() {
  for t in $*; do
    echo "$FILE" | grep -q "$t" || return 1
  done

  return 0
}

# Take a string, a prefix/suffix pair and a replacement pattern and return
# processed string
regex() {
  STRING="$1"
  PREFIX="$2"
  SUFFIX="$3"
  REPLACE="$4"

  STRING=$(echo "$STRING" | sed -e ":loop" -e "s/${PREFIX}${SUFFIX}/${REPLACE}/g" -e "t loop")
  echo "$STRING"
}

################################################################################
# Remove m_variables direct access pointer in zfsfvpartcont.cpp and zfsfvparticlebase.cpp
################################################################################

PREFIX=""
SUFFIX="m_variables = block->m_variables;"
REPLACE=""
if ftype zfsfvpart; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

PREFIX=""
SUFFIX="ZFSFvParticleBase::m_variables = m_variables;"
REPLACE=""
if ftype zfsfvpart; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

PREFIX=""
SUFFIX="ZFSFloat\*\* ZFSFvParticleBase::m_variables;"
REPLACE=""
if ftype zfsfvpart; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

################################################################################
# Patterns: m_block->m_variables[][]
################################################################################

PREFIX="m_block->"
SUFFIX="[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*"
for i in `seq 1 $NO_PATTERNS`; do
  SUFFIX="${SUFFIX}${PATTERN}"
done
REPLACE="m_block->${NEWNAME}(\1, \6)"
if ftype fv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

################################################################################
# Patterns: m_variables[][] in zfsfvpartcont*
################################################################################

PREFIX=""
SUFFIX="[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*"
for i in `seq 1 $NO_PATTERNS`; do
  SUFFIX="${SUFFIX}${PATTERN}"
done
REPLACE="blockPtr->${NEWNAME}(\1, \6)"
if ! ftype zfsfvblock.h && ftype fvpart; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

################################################################################
# Patterns: m_variables[][] in fv*
################################################################################

PREFIX="\([- +<=*/${TAB})(]\)"
SUFFIX="[ ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*"
for i in `seq 1 $NO_PATTERNS`; do
  SUFFIX="${SUFFIX}${PATTERN}"
done
REPLACE="\1${NEWNAME}(\2, \7)"
if ! ftype zfsfvblock.h && ftype fv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

################################################################################
# Patterns: &m_variables[][] in fv*
################################################################################

PREFIX="&"
SUFFIX="${OLDNAME}[${NEWLINE} ${TAB}]*"
for i in `seq 1 $NO_PATTERNS`; do
  SUFFIX="${SUFFIX}${PATTERN}"
done
REPLACE="\&${NEWNAME}(\1, \6)"
if ! ftype zfsfvblock.h && ftype fv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

PREFIX="&"
SUFFIX="${OLDNAME}[${NEWLINE} ${TAB}]*"
SUFFIX="${SUFFIX}${PATTERN}"

REPLACE="\&${NEWNAME}(\1, 0)"
if ! ftype zfsfvblock.h && ftype fv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

PREFIX="\([- +<=*/${TAB})(]\)"
SUFFIX="${OLDNAME}[${NEWLINE} ${TAB}]*"

REPLACE="\1pvariableptr()"
if ! ftype zfsfvblock*.h && ! ftype zfsfvsurface* && ftype zfsfv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

################################################################################
# Cleanup: write contents back to file and change back newlines
################################################################################

# Write content back to file and replace $NEWLINE with newlines
echo "$CONTENT" | sed "s/${NEWLINE}/\n/g" > $FILE
