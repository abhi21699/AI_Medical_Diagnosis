#!/bin/bash

if [ "$#" -ne 1]; then
    echo "Usage: git_commit <file_name>"
    exit 1
fi

file="$1"

git add "$file"

git commit -m "commit $file"

git push origin main