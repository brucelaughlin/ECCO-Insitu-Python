"""config_loader: merging order, validation (unknown keys, wrong types), path resolution, saving."""
import pytest
import yaml

import ecco_insitu.config as cl


def _load(tmp_path, *layers, overrides=()):
    files = []
    for i, d in enumerate(layers):
        f = tmp_path / f'layer{i}.yaml'
        f.write_text(yaml.safe_dump(d))
        files.append(f)
    return cl.load(files, overrides, use_local=False, use_env=False)


def test_defaults_alone_load(tmp_path):
    cfg = _load(tmp_path)
    assert cfg['chain']['new_floor'] == {'prof_T': 0, 'prof_S': 0.005}
    assert cfg['chain']['convert_insitu_T_to_potential'] == 'auto'
    assert cfg['high_latitude_sigma']['south_constant'] == {'prof_T': 10.0, 'prof_S': 5.0}


def test_later_layers_win(tmp_path):
    cfg = _load(tmp_path, {'chain': {'new_floor': {'prof_S': 0.01}}},
                {'chain': {'new_floor': {'prof_S': 0.02}}}, overrides=['chain.new_floor.prof_S=0.03'])
    assert cfg['chain']['new_floor'] == {'prof_T': 0, 'prof_S': 0.03}     # prof_T kept from defaults


def test_unknown_key_is_an_error(tmp_path):
    with pytest.raises(cl.ConfigError, match='unknown setting "chain.new_flor"'):
        _load(tmp_path, {'chain': {'new_flor': {'prof_S': 1}}})


@pytest.mark.parametrize('bad', [{'chain': {'apply_gamma_factor': 3}},          # bool expected
                                 {'chain': {'distance_tolerance': 'far'}},      # number expected
                                 {'chain': {'variables': 'prof_T'}},            # list expected
                                 {'chain': {'new_floor': 0.1}}])                # mapping expected
def test_wrong_type_is_an_error(tmp_path, bad):
    with pytest.raises(cl.ConfigError):
        _load(tmp_path, bad)


def test_convert_T_accepts_auto_or_bool(tmp_path):
    assert _load(tmp_path, {'chain': {'convert_insitu_T_to_potential': False}})['chain'][
        'convert_insitu_T_to_potential'] is False


def test_relative_paths_use_data_root(tmp_path):
    cfg = _load(tmp_path, {'data_root': str(tmp_path), 'paths': {'grid_dir': 'grid_llc90'}})
    assert cfg['paths']['grid_dir'] == str(tmp_path / 'grid_llc90')


def test_relative_path_without_data_root_is_an_error(tmp_path):
    with pytest.raises(cl.ConfigError, match='data_root is not set'):
        _load(tmp_path, {'paths': {'grid_dir': 'grid_llc90'}})


def test_missing_required_paths_are_listed(tmp_path):
    cfg = _load(tmp_path, {'paths': {'grid_dir': str(tmp_path / 'nope')}})
    with pytest.raises(cl.ConfigError) as e:
        cl.check_required_paths(cfg)
    assert 'paths.grid_dir: does not exist' in str(e.value) and 'paths.climatology_file: not set' in str(e.value)


def test_dump_round_trips(tmp_path):
    cfg = _load(tmp_path, overrides=['chain.method=1'])
    out = tmp_path / 'run' / 'run_config.yaml'
    cl.dump(cfg, out)
    back = yaml.safe_load(out.read_text())
    assert back['chain'] == cfg['chain'] and '_sources' not in back
    assert out.read_text().startswith('# Settings used for this run')


def test_example_local_file_is_valid(tmp_path):
    """config/local.example.yaml must only use real keys."""
    cl.load([cl.CONFIG_DIR / 'local.example.yaml'], use_local=False, use_env=False)


def test_defaults_ship_with_the_package():
    assert cl.DEFAULTS_FILE.parent.name == 'ecco_insitu' and cl.DEFAULTS_FILE.is_file()
