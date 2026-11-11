# -*- coding: utf-8 -*-
"""One-time path adapter for a fresh copy of this repository.

The analysis scripts were developed with absolute ROOT constants pointing at
the original workstation. Run this ONCE after extracting the package:

    python code/make_portable.py

It rewrites only those string literals to the location of this copy. It does
not touch any numeric input, result table or figure. Running it twice is safe.
"""
from pathlib import Path

CODE = Path(__file__).resolve().parent
REPO = CODE.parent
OLD_ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
OLD_PROJ = 'F:/博士/国自然科学基金/2025神自然'
TOKEN = '@@A1_REPO_ROOT_TOKEN@@'

new_root = REPO.as_posix()
# In a portable copy the project root IS the extracted repository, so every
# generated document stays inside this copy instead of leaking to its parent.
new_proj = REPO.as_posix()
(REPO / '_markdown源').mkdir(exist_ok=True)

changed = []
for p in sorted(CODE.glob('*.py')):
    if p.name == 'make_portable.py':
        continue
    s = p.read_text(encoding='utf-8')
    if OLD_ROOT not in s and OLD_PROJ not in s:
        continue
    # Longest match first via a sentinel, so the shorter project-root literal
    # cannot corrupt a repository-root path that happens to contain it.
    t = s.replace(OLD_ROOT, TOKEN)
    t = t.replace(OLD_PROJ, new_proj)
    t = t.replace(TOKEN, new_root)
    # a1_postreview_patch.py builds ROOT as PROJ + '/A1_启动'. In a portable
    # copy the repository root IS the project root, so that suffix must go.
    t = t.replace(" + '/A1_启动'", "")
    if t != s:
        p.write_text(t, encoding='utf-8')
        changed.append(p.name)

print('repository root :', new_root)
print('project root    :', new_proj)
print('adapted scripts :', len(changed))
for n in changed:
    print('  -', n)
print('Done. You can now run: bash code/run_all.sh')
