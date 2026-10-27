"""Step 06 code 12 (chain.xbt_max_depth, KI-21): XBT values deeper than the limit get zero weight; other
profiles and shallower values are untouched; None turns the test off."""
import numpy as np
import xarray as xr

from ecco_insitu.steps import step06


def _ds():
    z = np.array([10., 1000., 2000., 3000.])
    T = np.array([[20., 5., 3., 2.5], [20., 5., 3., 2.5]])
    return xr.Dataset({
        'prof_T': (('iPROF', 'iDEPTH'), T),
        'prof_Tweight': (('iPROF', 'iDEPTH'), np.ones_like(T)),
        'prof_Tclim': (('iPROF', 'iDEPTH'), T.copy()),
        'prof_Tuncertainty_source': (('iPROF',), np.zeros(2)),
        'prof_depth': (('iDEPTH',), z),
        'prof_descr': (('iPROF',), np.array([b'WOD-XBT//4128244:JSVY 1994', b'WOD-CTD//3395198:Au1191_'])),
        'prof_lat': (('iPROF',), np.array([-40., -40.])),
        'prof_lon': (('iPROF',), np.array([150., 150.])),
        'prof_YYYYMMDD': (('iPROF',), np.array([19940101., 19940101.])),
        'prof_HHMMSS': (('iPROF',), np.array([120000., 120000.])),
    })


def _run(limit):
    return step06.main(_ds(), {'prof_T'}, True, 60, True, xbt_max_depth=limit)


def test_xbt_values_below_limit_zeroed():
    out = _run(2900)
    assert list(out.prof_Tweight.values[0]) == [1, 1, 1, 0]          # XBT: 3000 m zeroed
    assert list(out.prof_Tweight.values[1]) == [1, 1, 1, 1]          # CTD: untouched
    assert int(out.prof_Tweight_code.values[0, 3]) & 2 ** 11


def test_no_limit():
    out = _run(None)
    assert (out.prof_Tweight.values == 1).all()
    assert not (out.prof_Tweight_code.values.astype(int) & 2 ** 11).any()
