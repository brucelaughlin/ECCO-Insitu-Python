"""Rebuild PFL vs_model T mean term at 2400 m per bin from the profile details, vs the saved analysis."""
import numpy as np, netCDF4, gzip, json
from pathlib import Path
import ecco_insitu.diagnostics.bin_viewer.export as e
D = Path('/Users/brucel/ecco/yip/bin_viewer_data/run_20261027_142247')
b = gzip.decompress((D/'PFL/platforms.bin').read_bytes()); n = int(np.frombuffer(b,'<u4',1)[0])
cb = e._unplanes(b[4:4+2*n], n, '<u2')
K = 78
sums = {}
for f in sorted((D/'detail/PFL').glob('c*.bin')):
    c = int(f.stem[1:]); X = e.decode_detail(f.read_bytes()); T = X['T']; start = np.r_[0, np.cumsum(T['cnt'])]
    sel = np.nonzero(T['lev'] == K)[0]
    if not len(sel): continue
    prof = np.searchsorted(start, sel, side='right') - 1
    ok = ((T['mask'][sel] & 2) > 0) & ((T['mask'][sel] & 4) > 0) & (T['code'][sel] == 0) & (T['sig'][sel] > 0)
    for i, q in zip(sel[ok], prof[ok]):
        j = c * e.DETAIL_CHUNK + q; bb = int(cb[j]); m = T['obs'][i] - T['model'][i]; s2 = T['sig'][i] ** 2
        a = sums.setdefault(bb, [0.0, 0.0, 0]); a[0] += m; a[1] += s2; a[2] += 1
ds = netCDF4.Dataset('/Users/brucel/ecco/yip/insitu_reports/profile_analyses_20261028_142247/PFL/analysis_PFL_20261027_142247_vs_model.nc'); ds.set_auto_mask(False)
ref = ds['meanterm__T__true__by_bin_depth'][:, K]
mine = {bb: (lambda mean, s, n: (1 / max(s / n, s / 9)) * mean ** 2)(a[0] / a[2], a[1] / a[2], a[2]) for bb, a in sums.items()}
common = [bb for bb in mine if 1 <= bb <= len(ref) and np.isfinite(ref[bb - 1])]
rel = np.array([abs(mine[bb] - ref[bb - 1]) / max(abs(ref[bb - 1]), 1e-9) for bb in common])
print(f'bins: details {len(mine)}, analysis {int(np.isfinite(ref).sum())}, both {len(common)}')
print(f'relative difference per bin: median {np.median(rel):.2e}, 99th pct {np.percentile(rel, 99):.2e}, max {rel.max():.2e}')
print(f'mean over bins: details {np.mean([mine[bb] for bb in common]):.3f}, analysis {np.nanmean(ref):.3f}')
bad = sorted(common, key=lambda bb: -abs(mine[bb] - ref[bb-1]) / max(abs(ref[bb-1]), 1e-9))[:6]
cnt = ds['bias_own_count_T'][:, K]; bias = ds['bias_own_T'][:, K]
for bb in bad:
    a = sums[bb]; print(f'bin {bb}: details n={a[2]} mean misfit {a[0]/a[2]:+.4f} -> {mine[bb]:.4g}; analysis n={int(cnt[bb-1])} bias {bias[bb-1]:+.4f} -> {ref[bb-1]:.4g}')
absd = np.array([abs(mine[bb] - ref[bb-1]) for bb in common]); tot = sum(ref[bb-1] for bb in common)
big = [bb for bb in common if ref[bb-1] > 0.1]
relb = np.array([abs(mine[bb] - ref[bb-1]) / ref[bb-1] for bb in big])
print(f'absolute difference per bin: max {absd.max():.4f}; summed {absd.sum():.3f} = {100*absd.sum()/tot:.3f} % of the summed term')
print(f'bins with term > 0.1 ({len(big)}): relative difference max {relb.max():.2e}, 99th pct {np.percentile(relb,99):.2e}')
print('counts equal in every bin:', all(sums[bb][2] == int(cnt[bb-1]) for bb in common))
