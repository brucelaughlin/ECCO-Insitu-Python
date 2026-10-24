"""The Sphinx documentation builds with no warnings (docs/sphinx), and every public function and
class in the package has a docstring."""
import ast
import importlib.util
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / 'src' / 'ecco_insitu'


def _public_without_docstring():
    missing = []
    for p in sorted(SRC.rglob('*.py')):
        tree = ast.parse(p.read_text())
        if not ast.get_docstring(tree):
            missing.append(f'{p.relative_to(ROOT)}: module')
        # what Sphinx documents: module-level functions and classes, and the methods of those classes
        # (functions nested inside functions are internal and not documented)
        nodes = list(tree.body)
        nodes += [m for c in tree.body if isinstance(c, ast.ClassDef) and not c.name.startswith('_') for m in c.body]
        for node in nodes:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
                    and not node.name.startswith('_') and not ast.get_docstring(node):
                missing.append(f'{p.relative_to(ROOT)}:{node.lineno} {node.name}')
    return missing


def test_every_public_function_has_a_docstring():
    missing = _public_without_docstring()
    assert not missing, f'{len(missing)} without docstrings:\n  ' + '\n  '.join(missing)


@pytest.mark.slow
@pytest.mark.skipif(importlib.util.find_spec('sphinx') is None or importlib.util.find_spec('myst_parser') is None,
                    reason='sphinx / myst-parser not installed (pip install -e ".[docs]")')
def test_docs_build_without_warnings(tmp_path):
    r = subprocess.run([sys.executable, '-m', 'sphinx', '-b', 'html', '-q', '-W', '--keep-going', '-E',
                        '-d', str(tmp_path / 'doctrees'), str(ROOT / 'docs' / 'sphinx'), str(tmp_path / 'html')],
                       capture_output=True, text=True)
    assert r.returncode == 0, 'Sphinx warnings or errors:\n' + r.stderr[-6000:]
