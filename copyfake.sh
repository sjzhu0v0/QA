# path of the script
SCRIPT_DIR=$(realpath "$(dirname "$(readlink -f "$0")")")
echo "Script directory: $SCRIPT_DIR"

# Print all the directories in the script directory
echo -n "Directories in script directory: "
for dir in "$SCRIPT_DIR"/*/; do
    if [ -d "$dir" ]; then
        echo -n "$(basename "$dir"), "
    fi
done
echo

# create soft link to the directories
for dir in "$SCRIPT_DIR"/*/; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        echo "Creating symlink for $dir_name"
        ln -s "$SCRIPT_DIR/$dir_name"
    fi
done
