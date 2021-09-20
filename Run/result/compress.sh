# tar and compress every folder

# 7z tar and compress
find . -maxdepth 1 -mindepth 1 -type d -exec 7zr a -mx=5 -mmt=on -sdel {}.7z {} \;

# remove empty folders
find . -maxdepth 1 -mindepth 1 -type d -empty -exec rmdir {}  \;
