"""
Check that two versions of Python files differ only in docstrings and comments.

Each file is parsed, every docstring (module, class, function) is removed from the syntax tree, and
the trees are compared. Comments and formatting never reach the tree. So "same" means the code
that runs is identical.

Usage:
    python tests/tools/same_code.py <git-ref> [files...]    # compare the working tree with <git-ref>
"""
import ast
import subprocess
import sys


def _strip(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], 'value', None), ast.Constant) \
                    and isinstance(body[0].value.value, str):
                node.body = body[1:] or [ast.Pass()]
    return ast.dump(tree, include_attributes=False)


def same(old_src, new_src):
    """True if the two sources are the same code apart from docstrings and comments."""
    return _strip(ast.parse(old_src)) == _strip(ast.parse(new_src))


def main(argv):
    ref, files = argv[0], argv[1:]
    if not files:
        files = subprocess.run(['git', 'diff', '--name-only', ref, '--', '*.py'], capture_output=True,
                               text=True, check=True).stdout.split()
    bad = []
    for f in files:
        try:
            old = subprocess.run(['git', 'show', f'{ref}:{f}'], capture_output=True, text=True, check=True).stdout
        except subprocess.CalledProcessError:
            print(f'new file (not compared): {f}')
            continue
        if not same(old, open(f).read()):
            bad.append(f)
    for f in bad:
        print(f'CODE CHANGED: {f}')
    print(f'{len(files)} files compared, {len(bad)} with code changes')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
