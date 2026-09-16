#!/usr/bin/env python3
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
failures = []
warnings = []
passes = []


def pass_(msg): passes.append(msg)
def warn(msg): warnings.append(msg)
def fail(msg): failures.append(msg)

# Required reviewer-facing files added by this preparation pass.
for rel in [
    'README.md', 'REVIEWER_CHECKLIST.md', 'docs/INSTALLATION.md',
    'docs/RUN_ON_YOUR_DATA.md', 'demo/run_demo.py', 'demo/README.md'
]:
    if (ROOT / rel).exists(): pass_(f'present: {rel}')
    else: fail(f'missing: {rel}')

# Known critical source dependency.
if (ROOT / '2_python_script/gtbulksearch_py3.py').exists() or (ROOT / 'gtbulksearch_py3.py').exists():
    pass_('gtbulksearch_py3.py is present')
else:
    fail('gtbulksearch_py3.py is absent but is called by 2_Somatic_mutation_calling/04_Variant_call.sh')

# The distributed extract_tags.py is a provenance stub, not the upstream implementation.
extract_stub = ROOT / '2_python_script' / 'extract_tags.py'
if extract_stub.exists() and 'Please refer to the following link' in extract_stub.read_text():
    warn('2_python_script/extract_tags.py is a provenance stub; fetch/install the pinned NanoSeq implementation before running preprocessing')

# Project license.
if any((ROOT / n).exists() for n in ['LICENSE', 'LICENSE.txt', 'LICENSE.md', 'COPYING']):
    pass_('project license file is present')
else:
    fail('project license file is absent')

# Static shell syntax, including the historical extensionless bulk script.
shell_paths = list(ROOT.rglob('*.sh'))
extensionless = ROOT / '1_Bulk_WGS' / '02_Bulk_variant_calling'
if extensionless.exists(): shell_paths.append(extensionless)
for path in sorted(set(shell_paths)):
    if '.ipynb_checkpoints' in path.parts: continue
    p = subprocess.run(['bash', '-n', str(path)], capture_output=True, text=True)
    if p.returncode == 0: pass_(f'bash syntax: {path.relative_to(ROOT)}')
    else: fail(f'bash syntax: {path.relative_to(ROOT)}: {p.stderr.strip()}')

# Python compilation.
for path in sorted(ROOT.rglob('*.py')):
    if '.ipynb_checkpoints' in path.parts: continue
    p = subprocess.run([sys.executable, '-m', 'py_compile', str(path)], capture_output=True, text=True)
    if p.returncode == 0: pass_(f'python compile: {path.relative_to(ROOT)}')
    else: fail(f'python compile: {path.relative_to(ROOT)}: {p.stderr.strip()}')

# Historical placeholders.
pat = re.compile(r'xxxxxx|xxxxxxxx|\bAYY\b|\bAXX\b|\bAXXX\b|\bBYY\b|\bAAA\b')
placeholder_files = []
for path in list((ROOT/'1_Bulk_WGS').glob('*')) + list((ROOT/'2_Somatic_mutation_calling').glob('*')) + list((ROOT/'3_Phylogeny').glob('*')):
    if path.is_file():
        try: txt = path.read_text()
        except UnicodeDecodeError: continue
        if pat.search(txt): placeholder_files.append(str(path.relative_to(ROOT)))
if placeholder_files:
    warn('study-specific placeholders remain in historical templates: ' + ', '.join(placeholder_files))

# Demo.
p = subprocess.run([sys.executable, str(ROOT/'demo/run_demo.py')], cwd=ROOT, capture_output=True, text=True)
if p.returncode == 0:
    pass_('reviewer demo passes')
else:
    fail('reviewer demo fails: ' + p.stdout + p.stderr)

print('=== PASS ===')
for x in passes: print('PASS:', x)
print('\n=== WARN ===')
for x in warnings: print('WARN:', x)
print('\n=== FAIL / AUTHOR ACTION REQUIRED ===')
for x in failures: print('FAIL:', x)
print(f'\nSummary: {len(passes)} pass, {len(warnings)} warn, {len(failures)} fail')
sys.exit(1 if failures else 0)
