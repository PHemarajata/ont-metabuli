#!/usr/bin/env python3
"""
Generate the self-contained Terra WDL from the template + the *tested* Python
scripts in ../bin.

Why generate instead of hand-copying: Terra works best with a single WDL file
that needs no auxiliary staging, so the analysis scripts are embedded in the
task command blocks. Generating them from bin/*.py guarantees the embedded
copies are byte-identical to the ones validated by the Nextflow pipeline, and
makes updating them a one-command operation.

Usage:
    python3 wdl/build_wdl.py            # writes wdl/ont_metabuli.wdl
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, 'bin')
TEMPLATE = os.path.join(HERE, 'ont_metabuli.wdl.tmpl')
OUTPUT = os.path.join(HERE, 'ont_metabuli.wdl')

MARKER = re.compile(r'^([ \t]*)@@INSERT:([A-Za-z0-9_.]+)@@[ \t]*$')


def main():
    with open(TEMPLATE) as fh:
        lines = fh.readlines()

    out = []
    inserted = []
    for line in lines:
        m = MARKER.match(line.rstrip('\n'))
        if not m:
            out.append(line)
            continue
        script = m.group(2)
        path = os.path.join(BIN, script)
        if not os.path.exists(path):
            sys.exit(f"ERROR: {path} not found (cannot embed {script})")
        with open(path) as sf:
            body = sf.read()
        # Emit at column 0: the surrounding heredoc is quoted (<<'PYEOF'), so
        # content is preserved literally and Python's own indentation is intact.
        if not body.endswith('\n'):
            body += '\n'
        out.append(body)
        inserted.append(script)

    text = ''.join(out)

    # Safety: WDL interpolates ~{...} inside command blocks. If an embedded
    # script ever contains that sequence it would be silently substituted.
    if '~{' in ''.join(open(os.path.join(BIN, s)).read() for s in inserted):
        sys.exit("ERROR: an embedded script contains '~{' which WDL would "
                 "interpolate. Refactor the script before embedding.")

    with open(OUTPUT, 'w') as fh:
        fh.write(text)

    print(f"wrote {OUTPUT}")
    print(f"embedded: {', '.join(inserted)}")
    print(f"{len(text.splitlines())} lines")


if __name__ == '__main__':
    main()
