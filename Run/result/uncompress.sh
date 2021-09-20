#! /bin/bash

for filename in ./*.7z; do
	7zr x $filename
done
