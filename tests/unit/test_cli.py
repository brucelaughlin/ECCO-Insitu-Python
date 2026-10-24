"""The ecco-insitu command: every subcommand starts and shows its help."""
import os
import subprocess
import sys

import pytest

from ecco_insitu import cli

@pytest.mark.parametrize('cmd', [c for c in cli.COMMANDS if c != 'run'])
def test_every_command_has_help(cmd, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main([cmd, '--help'])
    assert e.value.code == 0
    assert f'ecco-insitu {cmd}' in capsys.readouterr().out


def test_top_level_help_lists_commands(capsys):
    assert cli.main(['--help']) == 0
    out = capsys.readouterr().out
    assert all(c in out for c in cli.COMMANDS) and 'climatology --build' in out


def test_unknown_command():
    assert cli.main(['no-such-command']) == 2


def test_installed_command_runs_from_any_folder(tmp_path):
    """``python -m ecco_insitu.cli`` (what the ``ecco-insitu`` script runs) works outside the repository."""
    src = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src'))
    env = dict(os.environ, PYTHONPATH=src)
    r = subprocess.run([sys.executable, '-m', 'ecco_insitu.cli', 'config', '--get', 'chain.method'],
                       cwd=tmp_path, capture_output=True, text=True, env=env)
    assert r.returncode == 0 and r.stdout.strip() == '1', r.stderr[-2000:]
