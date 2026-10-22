#!/usr/bin/env bash
# Copyright (C) The IMR developers. SPDX-License-Identifier: GPL-3.0-or-later
# Prove malformed input fails, including math that Markdown leaves as prose.
set -euo pipefail
checker="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/check.sh"
fixture_tmp=$(mktemp -d)
trap 'rm -rf "$fixture_tmp"' EXIT
touch "$fixture_tmp/linked file.md"
cat > "$fixture_tmp/valid.md" <<'MD'
Inline $x^2$ and [a relative link](linked%20file.md#section).

$$
\frac{1}{2}+\sqrt{x}
$$

`$unclosed` is code, as is the deliberately invalid example below:

```matlab
bad = '$\imrUndefinedCommand{z}$';
[a missing example link](missing.md)
```
MD
bash "$checker" "$fixture_tmp/valid.md"

expect_failure() {
    local fixture=$1 expected=$2
    if bash "$checker" "$fixture" > "$fixture_tmp/result.log" 2>&1; then
        printf 'Checker incorrectly accepted %s\n' "$fixture" >&2
        exit 1
    fi
    if ! grep -F "$expected" "$fixture_tmp/result.log" >/dev/null; then
        cat "$fixture_tmp/result.log" >&2
        printf 'Checker failed for an unexpected reason: %s\n' "$fixture" >&2
        exit 1
    fi
}
printf '$\\imrUndefinedCommand{z}$\n' > "$fixture_tmp/undefined.md"
expect_failure "$fixture_tmp/undefined.md" 'Undefined control sequence'
printf '$\\frac{1}{2$\n' > "$fixture_tmp/braces.md"
expect_failure "$fixture_tmp/braces.md" 'unparsed dollar sign'
printf 'Valid $x$ followed by $unclosed\n' > "$fixture_tmp/delimiters.md"
expect_failure "$fixture_tmp/delimiters.md" 'unparsed dollar sign'
printf 'Valid $x$ with [missing target](missing.md).\n' > "$fixture_tmp/link.md"
expect_failure "$fixture_tmp/link.md" 'missing relative link'
printf 'No equations here.\n' > "$fixture_tmp/empty.md"
expect_failure "$fixture_tmp/empty.md" 'No mathematical expressions'
printf 'PASS: code exclusion, encoded links, invalid TeX, unmatched delimiters and empty math checks.\n'
