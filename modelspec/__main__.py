"""python -m modelspec build | check | impact | serve

build   parse input_data/ and the source tree -> modelspec/out/model_spec.json + ./4DWCM-GUI.html (needs pandas,
        openpyxl, biopython, libsbml: run it in the 4DWCM container)
check   validate a perturbation YAML against the spec and print the edits and static impact (standard library only)
impact  static impact of knocking out one or more genes
serve   serve the model browser on localhost; its Save button validates and writes perturbations/<name>.yaml
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HEAD = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'out')
GUI = os.path.join(HEAD, '4DWCM-GUI.html')      # the model browser, at the repo root


def _load_spec(path):
    if not os.path.exists(path):
        sys.exit('%s not found: run "python -m modelspec build" (in the container) first' % path)
    return json.load(open(path))


def cmd_build(a):
    from . import gui
    if a.gui_only:                                   # re-inline the page from the last spec (no container needed)
        spec = _load_spec(os.path.join(a.out, 'model_spec.json'))
        hp = gui.write_html(spec, GUI)
        gui.write_artifact(hp, os.path.join(a.out, '4DWCM-GUI.artifact.html'))
        print('wrote %s from the existing spec (+ modelspec/out/4DWCM-GUI.artifact.html for claude.ai)' % hp)
        return
    from .load import build_spec
    from .impact import single_knockouts
    spec = build_spec(a.head, with_code=not a.no_code)
    spec['single_knockouts'] = single_knockouts(spec)
    os.makedirs(a.out, exist_ok=True)
    jp = os.path.join(a.out, 'model_spec.json')
    with open(jp, 'w') as f:
        json.dump(spec, f, separators=(',', ':'))
    hp = gui.write_html(spec, GUI)
    gui.write_artifact(hp, os.path.join(a.out, '4DWCM-GUI.artifact.html'))
    n_proc_species = sum(len(p.get('species_roles', {})) for p in spec['processes'].values())
    print('genes %d  species %d  reactions %d (%d active)  processes %d (%d species roles)' % (
        len(spec['genes']), len(spec['species']), len(spec['reactions']),
        sum(r['active'] for r in spec['reactions'].values()), len(spec['processes']), n_proc_species))
    rd = spec.get('rdme')
    if rd:
        print('RDME: %d reactions in %d families (%d per-gene types); CME: %d transcription + %d tRNA-charging reactions; '
              '%d rate constants, %d diffusion profiles' % (rd['n_reactions'], len(rd['families']), sum(f['per_gene'] for f in rd['families'].values()),
              rd['n_cme_transcription'], rd['n_cme_charging'], len(rd['rate_constants']), len(rd['diffusion_profiles'])))
    else:
        print('RDME: not recorded (build in the container)')
    if 'code' in spec:
        print('code refs: %d ids referenced, %d anchors resolved' % (len(spec['code']['refs']), len(spec['code']['anchors'])))
        unused = {}
        for r in spec['inputs']:
            if not r['used']:
                unused.setdefault(r['file'], []).append(r['sheet'])
        print('inputs never read by the code:')
        for f, sheets in unused.items():
            whole = all(r['file'] != f or not r['used'] for r in spec['inputs'])
            print('  %s%s' % (f, '' if whole else ': ' + ', '.join(s for s in sheets if s)))
    print('wrote %s (%.1f MB)\n      %s (%.1f MB)' % (jp, os.path.getsize(jp) / 1e6, hp, os.path.getsize(hp) / 1e6))


def _print_resolution(spec, res):
    print('edits:')
    for e in res['edits']:
        print('  %-40s %s -> %s   [%s]%s' % (e['target'], e['from'], e['to'], e['where'], ('  ' + e['note']) if e['note'] else ''))
    imp = res['impact']
    for k in ('processes', 'reactions', 'metabolites'):
        if imp[k]:
            print('%s affected:' % k)
            for i, v in sorted(imp[k].items(), key=lambda kv: (kv[1]['status'] != 'blocked', kv[0])):
                print('  %-8s %-24s %s' % (v['status'], i, '; '.join(v['why'])))


def cmd_check(a):
    from . import perturbation as P
    spec = _load_spec(a.spec)
    pert = P.load(a.file)
    err, warn = P.validate(spec, pert)
    for w in warn:
        print('warning: ' + w)
    if err:
        for e in err:
            print('error: ' + e)
        sys.exit(1)
    print('%s: valid' % a.file)
    _print_resolution(spec, P.resolve(spec, pert))


def cmd_impact(a):
    from . import perturbation as P
    spec = _load_spec(a.spec)
    pert = P.empty('impact')
    pert['knockouts'] = [{'gene': g if g.startswith('JCVI') else 'JCVISYN3A_' + g.zfill(4), 'mode': a.mode} for g in a.genes]
    err, _ = P.validate(spec, pert)
    if err:
        sys.exit('\n'.join(err))
    _print_resolution(spec, P.resolve(spec, pert))


def cmd_serve(a):
    from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
    from . import perturbation as P
    spec = _load_spec(a.spec)
    if not os.path.exists(GUI):
        sys.exit('%s not found: run "python -m modelspec build" (or build --gui-only) first' % GUI)
    html = open(GUI, 'rb').read()
    os.makedirs(a.dir, exist_ok=True)

    class H(BaseHTTPRequestHandler):
        def _send(self, code, body, ctype='application/json'):
            b = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header('Content-Type', ctype)
            self.send_header('Content-Length', str(len(b)))
            self.end_headers()
            self.wfile.write(b)

        def do_GET(self):
            if self.path in ('/', '/index.html'):
                return self._send(200, html, 'text/html; charset=utf-8')
            if self.path == '/api/list':
                return self._send(200, {'dir': a.dir, 'files': sorted(f for f in os.listdir(a.dir) if f.endswith('.yaml'))})
            if self.path.startswith('/api/file/'):
                name = os.path.basename(self.path[len('/api/file/'):])
                p = os.path.join(a.dir, name)
                if not os.path.isfile(p):
                    return self._send(404, {'error': 'no such file'})
                return self._send(200, {'name': name, 'perturbation': P.load(p)})
            self._send(404, {'error': 'not found'})

        def do_POST(self):
            n = int(self.headers.get('Content-Length', 0))
            try:
                body = json.loads(self.rfile.read(n))
            except ValueError:
                return self._send(400, {'error': 'bad json'})
            pert = body.get('perturbation')
            err, warn = P.validate(spec, pert)
            if self.path == '/api/validate' or err:
                res = None if err else P.resolve(spec, pert)
                return self._send(200 if not err else 422, {'errors': err, 'warnings': warn, 'resolution': res})
            if self.path == '/api/save':
                p = os.path.join(a.dir, pert['name'] + '.yaml')
                if os.path.exists(p) and not body.get('overwrite'):
                    return self._send(409, {'errors': ['%s exists' % p], 'exists': True})
                with open(p, 'w') as f:
                    f.write(P.dump(pert))
                return self._send(200, {'saved': os.path.relpath(p, HEAD), 'warnings': warn})
            self._send(404, {'error': 'not found'})

        def log_message(self, fmt, *args):
            sys.stderr.write('[modelspec] ' + fmt % args + '\n')

    srv = ThreadingHTTPServer((a.host, a.port), H)
    print('model browser on http://%s:%d  (saving to %s)' % (a.host, a.port, a.dir))
    if a.host == '127.0.0.1':
        print('from your laptop: ssh -L %d:localhost:%d <this node>, then open http://localhost:%d' % (a.port, a.port, a.port))
    srv.serve_forever()


def main():
    ap = argparse.ArgumentParser(prog='python -m modelspec', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    spec_default = os.path.join(OUT, 'model_spec.json')
    b = sub.add_parser('build', help='parse inputs and source, write spec + browser')
    b.add_argument('--head', default=HEAD, help='4DWCM repo root (default: %(default)s)')
    b.add_argument('--out', default=OUT)
    b.add_argument('--no-code', action='store_true', help='skip the source-code index')
    b.add_argument('--gui-only', action='store_true', help='only rebuild 4DWCM-GUI.html from the existing model_spec.json')
    b.set_defaults(fn=cmd_build)
    c = sub.add_parser('check', help='validate a perturbation file')
    c.add_argument('file')
    c.add_argument('--spec', default=spec_default)
    c.set_defaults(fn=cmd_check)
    i = sub.add_parser('impact', help='static knockout impact')
    i.add_argument('genes', nargs='+', help='locus tags or numbers (JCVISYN3A_0415 or 415)')
    i.add_argument('--mode', default='full', choices=['full', 'expression_only', 'initial_only'])
    i.add_argument('--spec', default=spec_default)
    i.set_defaults(fn=cmd_impact)
    s = sub.add_parser('serve', help='serve the browser with save-to-repo')
    s.add_argument('--port', type=int, default=8765)
    s.add_argument('--host', default='127.0.0.1')
    s.add_argument('--dir', default=os.path.join(HEAD, 'perturbations'))
    s.add_argument('--spec', default=spec_default)
    s.set_defaults(fn=cmd_serve)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
