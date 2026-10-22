#!/usr/bin/env bash
# Copyright (C) The IMR developers. SPDX-License-Identifier: GPL-3.0-or-later
# Check source Markdown links and compile its math; no website is built.
set -euo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_root"
for dependency in git pandoc jq pdflatex; do
    command -v "$dependency" >/dev/null || {
        printf 'Missing documentation tool: %s\n' "$dependency" >&2
        exit 1
    }
done
doc_tmp=$(mktemp -d)
trap 'rm -rf "$doc_tmp"' EXIT

if (( $# )); then
    doc_files=("$@")
else
    mapfile -d '' -t doc_files < <(git ls-files --cached --others --exclude-standard -z -- '*.md')
fi
if (( ${#doc_files[@]} == 0 )); then
    printf 'No Markdown files to check.\n' >&2
    exit 1
fi

cat > "$doc_tmp/equations.tex" <<'TEX'
\documentclass{article}
\usepackage{amsmath,amssymb}
\pagestyle{empty}
\begin{document}
TEX
math_count=0
link_count=0
for source in "${doc_files[@]}"; do
    pandoc --from=gfm+tex_math_dollars --to=json "$source" > "$doc_tmp/source.json"
    # An unmatched dollar would otherwise become ordinary prose and escape
    # compilation. Literal dollars in this documentation belong in code.
    if jq -e '.. | objects | select(.t? == "Str") | select(.c | contains("$"))' \
        "$doc_tmp/source.json" >/dev/null; then
        printf '%s: unparsed dollar sign; check math delimiters or put a literal dollar in code.\n' "$source" >&2
        exit 1
    fi
    while IFS= read -r link; do
        if [[ -z "$link" || "$link" == \#* || "$link" == //* || "$link" =~ ^[a-zA-Z][a-zA-Z0-9+.-]*: ]]; then
            continue
        fi
        target=${link%%[?#]*}
        # Decode percent escapes without evaluating link text as shell code.
        decoded=''
        while [[ "$target" =~ ^([^%]*)%([[:xdigit:]]{2})(.*)$ ]]; do
            printf -v character '%b' "\\x${BASH_REMATCH[2]}"
            decoded+="${BASH_REMATCH[1]}$character"
            target=${BASH_REMATCH[3]}
        done
        target="$decoded$target"
        if [[ "$target" == /* ]]; then
            destination="$repo_root$target"
        else
            destination="$(dirname -- "$source")/$target"
        fi
        if [[ ! -e "$destination" ]]; then
            printf '%s: missing relative link or image: %s\n' "$source" "$link" >&2
            exit 1
        fi
        link_count=$((link_count+1))
    done < <(jq -r '.. | objects | select(.t? == "Link" or .t? == "Image") | .c[-1][0]' "$doc_tmp/source.json")

    count=$(jq '[.. | objects | select(.t? == "Math")] | length' "$doc_tmp/source.json")
    math_count=$((math_count+count))
    jq -r --arg source "$source" '
        [.. | objects | select(.t? == "Math")] | to_entries[] |
        "\\typeout{IMR math: \($source) expression \(.key+1)}\n" +
        (if .value.c[0].t == "DisplayMath" then "\\[\n" else "\\(\n" end) +
        (.value.c[1] | gsub("^[\\r\\n\\t ]+|[\\r\\n\\t ]+$"; "") |
            gsub("\n[\t\r ]*\n+"; "\n")) +
        (if .value.c[0].t == "DisplayMath" then "\n\\]\n" else "\n\\)\n\n" end)
    ' "$doc_tmp/source.json" >> "$doc_tmp/equations.tex"
done
if (( math_count == 0 )); then
    printf 'No mathematical expressions found; refusing an empty equation check.\n' >&2
    exit 1
fi
printf '\\end{document}\n' >> "$doc_tmp/equations.tex"
if ! pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error \
    -file-line-error -output-directory="$doc_tmp" "$doc_tmp/equations.tex" \
    > "$doc_tmp/compile.log" 2>&1; then
    tail -60 "$doc_tmp/compile.log" >&2
    exit 1
fi
printf 'PASS: %s Markdown files, %s relative destinations, %s LaTeX expressions.\n' \
    "${#doc_files[@]}" "$link_count" "$math_count"
