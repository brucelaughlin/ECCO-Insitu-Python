"""Every value of a sample of profiles in the detail files vs the yearly file (and the model files)."""
import sys, glob, numpy as np, netCDF4
from pathlib import Path
import ecco_insitu.diagnostics.bin_viewer.export as e
from ecco_insitu.diagnostics.profile_analyses.compute import EstimIndex
from ecco_insitu.diagnostics import common
OUT, SRC = Path(sys.argv[1]), sys.argv[2]
RUN = common.yearly_dir('/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_20261027_142247/yearly')
est = EstimIndex('/Users/brucel/ecco/yip/ecco_v4r6_model_data/model_equivalents_20261028_142247/netCDF', RUN)
files = [f for f in sorted(RUN.rglob('*.nc')) if common.classify(f, RUN)[0] == SRC]
tree, _ = e._geodesic_tree()
fi, ri, lon, lat = e._ordered_rows(SRC, files, tree)
rng = np.random.default_rng(1); sample = np.sort(rng.choice(len(fi), min(400, len(fi)), replace=False))
chunks = {}; worst = {k: 0.0 for k in ('obs', 'clim', 'model', 'sig')}; nval = 0; bad = []
for j in sample:
    c = j // e.DETAIL_CHUNK
    if c not in chunks: chunks[c] = e.decode_detail((OUT / 'detail' / SRC / f'c{c:04d}.bin').read_bytes())
    D = chunks[c]; q = j - c * e.DETAIL_CHUNK; f = files[fi[j]]; r = ri[j]
    with netCDF4.Dataset(f) as ds:
        ds.set_auto_mask(False); n = ds.dimensions['iPROF'].size
        assert D['date'][q] == ds['prof_YYYYMMDD'][r], ('date', j)
        assert D['descr'][q] == b''.join(ds['prof_descr'][r]).decode(errors='replace').strip(), ('descr', j)
        m = est.get(f, n, ds.dimensions['iDEPTH'].size)
        for V in ('T', 'S'):
            if f'prof_{V}' not in ds.variables: continue
            X = D[V]; a = int(X['cnt'][:q].sum()); b = a + int(X['cnt'][q])
            v = ds[f'prof_{V}'][r].astype(float); have = np.isfinite(v) & (np.abs(v) < 900)
            if not np.array_equal(np.nonzero(have)[0], X['lev'][a:b]): bad.append(('levels', j, V)); continue
            ref = {'obs': v[have], 'clim': ds[f'prof_{V}clim'][r].astype(float)[have],
                   'sig': ds[f'prof_{V}uncertainty'][r].astype(float)[have] if f'prof_{V}uncertainty' in ds.variables else None,
                   'model': np.asarray(m[V][r], float)[have] if m else None}
            code = ds[f'prof_{V}weight_code'][r][have].astype(int)
            if not np.array_equal(code, X['code'][a:b]): bad.append(('code', j, V))
            for k, bit in (('obs', None), ('clim', 1), ('model', 2), ('sig', 4)):
                if ref[k] is None: continue
                ok = np.isfinite(ref[k]) & (np.abs(ref[k]) < 900) & ((ref[k] != 0) if k == 'model' else True)
                got = X[k][a:b]
                if bit is not None and not np.array_equal(((X['mask'][a:b] & bit) > 0), ok): bad.append(('mask', j, V, k))
                if ok.any(): worst[k] = max(worst[k], float(np.max(np.abs(got[ok] - ref[k][ok]))))
            nval += int(have.sum())
print(f'{SRC}: {len(sample)} profiles, {nval} values checked; problems: {bad[:5] if bad else "none"}')
print('  largest |stored - file|:', {k: f'{v:.6f}' for k, v in worst.items()}, '(rounding: 0.0005; sigma 0.00005)')
