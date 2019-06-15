#! /bin/bash

# Compiles process tex files from docs/processes/src to pdf

# Usage:
#   compile_process_pdfs.sh path_to_tex_file

set -eu

src_dir=$(dirname $1)
file_base=$(basename $1 .tex)

cd $src_dir

# Compile pdf from tex with proper citation linking
pdflatex ${file_base}.tex
bibtex ${file_base}.aux
pdflatex ${file_base}.tex
pdflatex ${file_base}.tex

# Move compiled pdf to parent directory
mv ${file_base}.pdf ../

# Clean intermediate files
rm *.log *.aux *.bbl *.blg
