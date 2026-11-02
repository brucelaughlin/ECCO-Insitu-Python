"""run_provenance: input fingerprints detect changed inputs; the record names the commit and saves a patch."""
import json
import os

import ecco_insitu.run_provenance as pv


def test_fingerprint_changes_with_inputs(tmp_path):
    d = tmp_path / 'inputs'
    d.mkdir()
    (d / 'a.nc').write_bytes(b'1234')
    before = pv.fingerprint(d)
    assert before['n_files'] == 1 and before['bytes'] == 4
    assert pv.fingerprint(d) == before                       # no change, same fingerprint
    (d / 'b.nc').write_bytes(b'5')
    assert pv.fingerprint(d)['listing_sha256'] != before['listing_sha256']
    os.remove(d / 'b.nc')
    os.utime(d / 'a.nc', (1, 1))                               # same name and size, other time
    assert pv.fingerprint(d)['listing_sha256'] != before['listing_sha256']


def test_write_and_compare(tmp_path):
    inp = tmp_path / 'in'
    inp.mkdir()
    (inp / 'x.nc').write_bytes(b'abc')
    run = tmp_path / 'run'
    rec = pv.write(run, {'paths': {'climatology_file': str(inp / 'x.nc')}}, inp)
    saved = json.loads((run / 'documentation' / 'provenance.json').read_text())
    assert saved['code']['commit'] == rec['code']['commit']
    assert set(saved['inputs']) == {'input_dir', 'paths.climatology_file'}
    if rec['code'].get('uncommitted_changes'):
        assert (run / 'documentation' / 'code.patch').exists()
    assert not [d for d in pv.compare(run) if not d.startswith('code:')]
    (inp / 'x.nc').write_bytes(b'abcd')
    assert any('input_dir' in d for d in pv.compare(run))
