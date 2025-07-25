"""Attempts to format short if, elseif, else, and for statements:

If the developer places braces, short statements are not allowed on the same line, e.g.:

```
    if (true) {
      do_something();
    }
```

If the developer does NOT place braces, short statements are preferred
on the same line, under the condition that the line is less than a given number of characters, e.g.:

```
    if (true) do_something();
```

Note, this assumes that clang-format has already been run with the options:

```
AllowShortLoopsOnASingleLine: true
AllowShortIfStatementsOnASingleLine: AllIfsAndElse
```

This will force both of the above to 1-liners.  This script will then
insert newlines in the first case.
"""

#
# Code written by chatgpt, tested by Ethan Coon
#

import re
import sys

NUM_CHARACTERS = 100

import re


import re

def reformatCppCode(lines):
    """Reformats lines according to the following rule:

    If the developer places braces, short statements are not allowed
    on the same line, e.g.:

    ```
    if (true) {
      do_something();
    }
    ```

    If the developer does NOT place braces, short statements are
    preferred on the same line, under the condition that the line is
    less than a given number of characters, e.g.:
    
    ```
    if (true) do_something();
    ```

    Note, this assumes that clang-format has already been run with the
    options:

    ```
    AllowShortLoopsOnASingleLine: true
    AllowShortIfStatementsOnASingleLine: AllIfsAndElse
    ```

    This will force both of the above to 1-liners.  This script will then
    insert newlines in the first case.

    """
    output_lines = []

    for line in lines:
        # Preserve the original line for indentation
        indent_match = re.match(r'^(\s*)', line)
        base_indent = indent_match.group(1)

        # Strip trailing whitespace only (not leading!)
        line = line.rstrip()

        # Match: if/else if/for (...) { statement; }
        m = re.match(r'^(\s*(if|else if|for)\s*\(.*?\))\s*{\s*(.*?)\s*};?\s*$', line)
        if m:
            output_lines.append(f"{m.group(1)} {{")
            output_lines.append(f"{base_indent}  {m.group(3).strip()}")
            output_lines.append(f"{base_indent}}}")
            continue

        # Match: else { statement; }
        m = re.match(r'^(\s*else)\s*{\s*(.*?)\s*};?\s*$', line)
        if m:
            output_lines.append(f"{m.group(1)} {{")
            output_lines.append(f"{base_indent}  {m.group(2).strip()}")
            output_lines.append(f"{base_indent}}}")
            continue

        # Match: if/else if/for (...) statement; (no braces)
        m = re.match(r'^(\s*(if|else if|for)\s*\(.*?\))\s*(?!{)(.+);$', line)
        if m:
            control = m.group(1)
            stmt = m.group(3).strip()
            full_line = f"{control} {stmt};"
            if len(full_line) <= 100:
                output_lines.append(full_line)
            else:
                output_lines.append(control)
                output_lines.append(f"{base_indent}  {stmt};")
            continue

        # Match: else statement; (no braces)
        m = re.match(r'^(\s*else)\s+(.+);$', line)
        if m:
            else_kw = m.group(1)
            stmt = m.group(2).strip()
            full_line = f"{else_kw} {stmt};"
            if len(full_line) <= 100:
                output_lines.append(full_line)
            else:
                output_lines.append(else_kw)
                output_lines.append(f"{base_indent}  {stmt};")
            continue

        # Pass through everything else unchanged
        output_lines.append(line)

    return output_lines


def main():
    if len(sys.argv) != 2:
        print("Usage: python reformat_cpp.py <file.cpp>.  Acts in-place on file.cpp")
        return

    input_file = sys.argv[1]

    with open(input_file, 'r') as f:
        lines = f.readlines()

    formatted_lines = reformatCppCode(lines)

    with open(input_file, 'w') as f:
        f.write('\n'.join(formatted_lines) + '\n')
        

if __name__ == "__main__":
    main()
