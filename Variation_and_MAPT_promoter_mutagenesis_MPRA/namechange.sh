#!/bin/bash
# Replace all occurences of becky with oligoMPRA, bri with internal
# Create a backup first (filename.oldnames)

FILE="$1"
TEMP="temp_${FILE}"
NEWFILE="new_${FILE}"

# Make a backup first
cp "$FILE" "$TEMP"
#replace names
sed -i '' -e ' ''/becky/ s/_scrambled//g' "$TEMP"
sed -i '' -e ' ''/briNegative/ s/_scrambled//g' "$TEMP"
sed -i '' -e 's/becky/oligoMPRA/g' -e 's/bri/internal/g' "$TEMP"
sed -i '' -e 's/Negative([^_])/Negative_scramble\1/g' "$TEMP"


cat "$TEMP" > "$NEWFILE"

rm "$TEMP"

echo "Replacements complete in $FILE"
echo "A new file has been created as $NEWFILE"
echo "Temporary files removed"











