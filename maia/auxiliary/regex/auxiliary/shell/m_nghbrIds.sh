#!/bin/sh

################################################################################
# Settings
################################################################################

# Old name, new name, and new name for G cells
OLDNAME="m_nghbrIds"
NEWNAME="a_neighborId"
NEWNAMEG="a_neighborGId"

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

SUFFIX="\.[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*"
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

  STRING=$(echo "$STRING" | sed -e ":loop" -e "s/${PREFIX}${PATTERN}${SUFFIX}/${REPLACE}/g" -e "t loop")
  echo "$STRING"
}


################################################################################
# Patterns: m_cells, cells
################################################################################

# Specific replacements for files with m_block
if [ $GRID_SCOPE -eq 1 ]; then
  REPLACE="grid().${NEWNAME}(\2, \7)"
else
  REPLACE="m_block->${NEWNAME}(\2, \7)"
fi
PREFIX="m_block->${PREF1}"
if ftype bnd || ftype interface; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi
PREFIX="${PREF1}"
if ftype bnd || ftype interface; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

# Specific replacements for files with blockPtr
if [ $GRID_SCOPE -eq 1 ]; then
  REPLACE="grid().${NEWNAME}(\2, \7)"
else
  REPLACE="blockPtr->${NEWNAME}(\2, \7)"
fi
PREFIX="blockPtr->${PREF1}"
if ftype fvpartcont || ftype fvparticlebase; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi
PREFIX="${PREF1}"
if ftype fvpartcont || ftype fvparticlebase; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

# Specific replacements for files with m_blockPtr
PREFIX="m_blockPtr->${PREF1}"
REPLACE="m_blockPtr->${NEWNAME}(\2, \7)"
if ftype zfslbmpartcont || ftype zfslbmparticle || ftype zfspartcont || ftype zfsparticle; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi
PREFIX="${PREF1}"
REPLACE="m_blockPtr->${NEWNAME}(\2, \7)"
if ftype zfslbmpartcont || ftype zfslbmparticle || ftype zfspartcont || ftype zfsparticle; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

# Specific replacements for files with grid-> pointer
PREFIX="grid->${PREF1}"
REPLACE="grid->${NEWNAME}(\2, \7)"
if ftype postprocessingblock; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

# Replacements for all other files
PREFIX="${PREF1}"
if ftype cartesiangrid || ftype avgblock; then
  REPLACE="${NEWNAME}(\2, \7)"
else
  REPLACE="${GRID}${NEWNAME}(\2, \7)"
fi
CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")


################################################################################
# Patterns: input_cells, cpu_cells
################################################################################
PREFIX="${PREF2}"
REPLACE="${NEWNAME}(\2, \3, \8)"
if ftype cartesiangrid; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi

SUFFIX2="\.[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*${PATTERN}"
if ftype cartesiangrid; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX2" "$REPLACE")
fi

################################################################################
# Patterns: cellsInput
################################################################################
PREFIX="${PREF3}"
REPLACE="${GRID}${NEWNAME}(\2, \3, \8)"
if ftype fvblock; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi


################################################################################
# Patterns: gCells, m_gCells
################################################################################

PREFIX="${PREF4}"
SUFFIX="\.[${NEWLINE} ${TAB}]*${OLDNAME}[${NEWLINE} ${TAB}]*${PATTERN}"
REPLACE="${NEWNAMEG}(\2, \7)"
if ftype fv; then
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi


################################################################################
# Patterns: pointer access
################################################################################

# Since this is rather seldom, only perform the following regexes if old name
# still exists at all
if echo "$CONTENT" | grep -q "$OLDNAME"; then
  SUFFIX="\.${OLDNAME}"

  # Specific replacements for files with grid-> pointer
  PREFIX="grid->${PREF1}"
  REPLACE="\&grid->${NEWNAME}(\2, 0)"
  if ftype postprocessingblock; then
    CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
  fi

  # Specific replacements for files with block ptr
  PREFIX="m_blockPtr->${PREF1}"
  REPLACE="\&m_blockPtr->${NEWNAME}(\2, 0)"
  if ftype zfslbmpartcont || ftype zfslbmparticle || ftype zfspartcont || ftype zfsparticle; then
    CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
  fi

  # For cells/m_cells
  PREFIX="${PREF1}"

  # Specific replacements for files with m_block
  REPLACE="\&m_block->${NEWNAME}(\2, 0)"
  if ftype bnd || ftype interface; then
    CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
  fi

  # Specific replacements for files with block ptr
  REPLACE="\&blockPtr->${NEWNAME}(\2, 0)"
  if ftype fvpartcont || ftype fvparticlebase; then
    CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
  fi

  # Replacements for all other files
  REPLACE="\&${NEWNAME}(\2, 0)"
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")

  # Replacements for gCells
  PREFIX="${PREF4}"
  REPLACE="\&${NEWNAMEG}(\2, 0)"
  CONTENT=$(regex "$CONTENT" "$PREFIX" "$SUFFIX" "$REPLACE")
fi


################################################################################
# Cleanup: write contents back to file and change back newlines
################################################################################

# Write content back to file and replace $NEWLINE with newlines
echo "$CONTENT" | sed "s/${NEWLINE}/\n/g" > $FILE
