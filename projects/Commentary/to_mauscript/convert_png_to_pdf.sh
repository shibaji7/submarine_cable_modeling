#!/usr/bin/env bash
set -u

shopt -s nullglob

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir" || exit 1

png_files=(*.png)

if ((${#png_files[@]} == 0)); then
    echo "No PNG files found in $script_dir"
    exit 0
fi

failed=0

for png_file in "${png_files[@]}"; do
    pdf_file="${png_file%.png}.pdf"
    echo "Converting $png_file -> $pdf_file"

    if ! convert "$png_file" "$pdf_file"; then
        echo "Failed: $png_file" >&2
        failed=1
    fi
done

exit "$failed"
