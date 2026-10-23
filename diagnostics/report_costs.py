"""
Cost vs climatology report: recompute (V - Vclim)^2 * Vweight for T and S from
the final chain files, store it per profile, and write an HTML + PDF report.

    python diagnostics/report_costs.py                         # newest finished step-11 run (~1 min + report)
    python diagnostics/report_costs.py --processed <run_dir>
    python diagnostics/report_costs.py --report_only [<cost_dir>]   # re-make report from existing costs

Output: <diagnostics.cost_report_dir>/cost_<run_ts>__<YYYYMMDD_HHMMSS>/ (from the chain settings)
See diagnostics/lib/costs_compute.py for definitions (fill rules, final vs
true unscaled cost) and the stored files.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
import costs_compute  # noqa: E402
import costs_report   # noqa: E402

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='T/S cost vs climatology report.')
    ap.add_argument('cost_dir', nargs='?', default=None, help='With --report_only: existing cost dir')
    ap.add_argument('--processed', default=None, help='Run dir (default: newest finished step-11 run)')
    ap.add_argument('--out', default=None)
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--sources', nargs='*', default=None, help='Only these sources (quick tests)')
    ap.add_argument('--report_only', action='store_true', help='Skip recomputing; re-make the report')
    ap.add_argument('--log', default=None,
                    help='Chain log for the decimation section (default: z_logs/ncei_log_<run ts>.log)')
    a = ap.parse_args()

    if a.report_only:
        cost_dir = Path(a.cost_dir) if a.cost_dir else max(costs_compute.REPORT_PARENT.glob('cost_*'),
                                                            key=lambda p: p.name)
    else:
        cost_dir = costs_compute.build(a.processed, a.out, a.n_workers, a.sources)
    h = costs_report.write_report(cost_dir, a.log)
    print(f'HTML: {h}\nPDF:  {costs_report.to_pdf(h)}')
