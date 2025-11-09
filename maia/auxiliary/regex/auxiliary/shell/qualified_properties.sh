#!/bin/bash

# File to be changed
FILE=$1

# Exit early if file does not use properties
grep -q "b_properties" $FILE || exit 0

# Temporary newline replacement (must be present in patterns below where newlines are allowed!)
NEWLINE="@@"

# Save tab for convenience
TAB="	"

# Save whitespace pattern for convenience
WS="[ $TAB@]"

# Read in file and change newlines to $NEWLINE
CONTENT=$(cat $FILE | sed ':a; N; $! b a; s/\n/'"$NEWLINE"'/g')

# Define regex function
regex() {
  # Save arguments
  STRING="$1"
  NUM="$2"
  REPL="$3"

  # Call sed
  STRING=$(echo "$STRING" | sed \
      -e "s/b_properties${WS}*\[${WS}*${NUM}${WS}*\]/b_properties[Cell::$REPL]/g")
  echo "$STRING"
}

# Apply regex
echo "Switching to qualified property names in $FILE..."
CONTENT=$(regex "$CONTENT" 0 IsDummy)
CONTENT=$(regex "$CONTENT" 1 HasChildren)
CONTENT=$(regex "$CONTENT" 2 IsBndry)
CONTENT=$(regex "$CONTENT" 3 IsGhost)
CONTENT=$(regex "$CONTENT" 4 IsInterface)
CONTENT=$(regex "$CONTENT" 5 IsPeriodic)
CONTENT=$(regex "$CONTENT" 6 IsCutOff)
CONTENT=$(regex "$CONTENT" 7 IsNotGradient)
CONTENT=$(regex "$CONTENT" 8 IsValidMGC)
CONTENT=$(regex "$CONTENT" 9 IsExchange)
CONTENT=$(regex "$CONTENT" 10 IsFlux)
CONTENT=$(regex "$CONTENT" 11 IsActive)
CONTENT=$(regex "$CONTENT" 12 MayBeCoarsened)
CONTENT=$(regex "$CONTENT" 13 IsOnCurrentMGLevel)
CONTENT=$(regex "$CONTENT" 14 IsInSpongeLayer)
CONTENT=$(regex "$CONTENT" 15 IsHalo)
CONTENT=$(regex "$CONTENT" 16 IsWindow)
CONTENT=$(regex "$CONTENT" 17 IsWindowOrHalo)
CONTENT=$(regex "$CONTENT" 18 IsWindowToDelete)
CONTENT=$(regex "$CONTENT" 19 IsHaloToDelete)
CONTENT=$(regex "$CONTENT" 20 AtStructuredRegion)
CONTENT=$(regex "$CONTENT" 21 IsPeriodicWithRot)
CONTENT=$(regex "$CONTENT" 22 IsParentGap)
CONTENT=$(regex "$CONTENT" 23 IsMinCell)
CONTENT=$(regex "$CONTENT" 24 IsSplitCell)
CONTENT=$(regex "$CONTENT" 25 IsSplitChild)
CONTENT=$(regex "$CONTENT" 26 HasSplitFace)
CONTENT=$(regex "$CONTENT" 27 IsLsInactive)
CONTENT=$(regex "$CONTENT" 28 WasLsInactive)
CONTENT=$(regex "$CONTENT" 29 IsTempLinked)

# Write content back to file and replace $NEWLINE with newlines
echo "$CONTENT" | sed "s/${NEWLINE}/\n/g" > $FILE
