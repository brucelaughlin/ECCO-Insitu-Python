"""The ecco-insitu command and the old script names: every command starts and shows its help."""
import os
import subprocess
import sys

import pytest

from ecco_insitu import cli

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OLD_SCRIPTS = ['NCEI.py', 'step11_spatial_scaling.py', 'step12_profiles_compact.py', 'step13_pack_ecco_arrays.py',
               'optimize_profile_files.py', 'validate_inputs.py', 'setup_climatology.py', 'config_loader.py',
               'build_woa23_climatology.py', 'prebake_woa23_climatology.py', 'csv_to_nc.py',
               'diagnostics/check_run_files.py', 'diagnostics/check_ecco_arrays.py', 'diagnostics/check_run_qc.py',
               'diagnostics/report_profile_counts.py', 'diagnostics/report_costs.py',
               'diagnostics/report_spatial_scaling.py']


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


@pytest.mark.parametrize('script', OLD_SCRIPTS)
def test_old_script_names_still_work(script, tmp_path):
    """Run from another folder, so nothing depends on the current directory."""
    r = subprocess.run([sys.executable, os.path.join(REPO, script), '--help'], cwd=tmp_path,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
