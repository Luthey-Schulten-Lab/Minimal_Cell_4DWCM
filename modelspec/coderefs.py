"""Index the simulator's source for everything the model browser shows in developer mode.

- literal references: every string literal equal to a species or reaction id ('P_0001', 'M_atp_c', 'RNAP'), with file,
  line, enclosing function and whether the line is commented out
- pattern references: literals that build ids for every gene ('P_' + locusNum, 'G_' + locusNum + '_C1')
- input reads: every sheet_name= and input_data/ path, which also tells which spreadsheet sheets the code never reads
- process anchors: 'path::function' from roles.PROCESSES resolved to line numbers; a missing function fails the build
"""

import ast
import glob
import io
import os
import re
import tokenize

SOURCE_GLOBS = ['Whole_Cell_Minimal_Cell.py', 'Restart_Whole_Cell_Minimal_Cell.py', 'cythonCompiledFunctions.pyx',
                'processes/*.py', 'processes/*.pyx', 'modules/*.py', 'utility/*.py', 'restart/*.py']

PREFIXES = {'P_': 'protein', 'C_P_': 'membrane-protein precursor', 'R_': 'RNA', 'G_': 'gene copy', 'RB_': 'ribosome-bound mRNA',
            'RP_': 'RNAP on gene', 'S_': 'SecY-bound protein', 'D_': 'degradosome-bound mRNA', 'DT_': 'degradation tracker',
            'PM_': 'proteins-made counter', 'RPM_': 'transcripts-made counter', 'DM_': 'mRNA-degraded counter', 'M_': 'metabolite'}
PREFIX_KINDS = {'protein': ['P_', 'C_P_', 'S_', 'PM_'], 'mRNA': ['R_', 'G_', 'RB_', 'RP_', 'D_', 'DT_', 'RPM_', 'DM_'],
                'tRNA': ['R_', 'G_', 'RP_', 'RPM_'], 'rRNA': ['R_', 'G_', 'RP_', 'RPM_'],
                'metabolite': ['M_'], 'medium': ['M_']}

QUOTED = re.compile(r"""(?P<q>['"])(?P<s>[^'"\n]{1,120})(?P=q)""")
SHEET = re.compile(r"""sheet_name\s*=\s*['"]([^'"]+)['"]""")
INPUT = re.compile(r"""input_data/([A-Za-z0-9_.\-]+)""")


class _Spans(ast.NodeVisitor):
    def __init__(self):
        self.spans, self.stack = [], []

    def _fn(self, node):
        self.stack.append(node.name)
        end = getattr(node, 'end_lineno', None) or max(getattr(n, 'lineno', node.lineno) for n in ast.walk(node))  # py3.7
        self.spans.append(('.'.join(self.stack), node.lineno, end))
        self.generic_visit(node)
        self.stack.pop()

    visit_FunctionDef = visit_AsyncFunctionDef = visit_ClassDef = _fn


def _function_spans(path, text):
    if path.endswith('.py'):
        try:
            v = _Spans()
            v.visit(ast.parse(text))
            return v.spans
        except SyntaxError:
            pass
    spans, lines = [], text.splitlines()          # .pyx or unparsable: top-level defs by indentation
    starts = [(i + 1, m.group(2)) for i, l in enumerate(lines) for m in [re.match(r'^(c?p?def|class)\s+(?:\w+\s+)*?(\w+)\s*\(', l)] if m]
    for k, (ln, name) in enumerate(starts):
        spans.append((name, ln, (starts[k + 1][0] - 1) if k + 1 < len(starts) else len(lines)))
    return spans


def _enclosing(spans, line):
    best = None
    for name, a, b in spans:
        if a <= line <= b and (best is None or a >= best[1]):
            best = (name, a, b)
    return best[0] if best else '<module>'


def _literals(path, text):
    """(line, literal, commented) for every quoted string, from tokens for .py (comments scanned too), regex otherwise."""
    out = []
    if path.endswith('.py'):
        try:
            for tok in tokenize.generate_tokens(io.StringIO(text).readline):
                if tok.type == tokenize.STRING:
                    try:
                        s = ast.literal_eval(tok.string)
                    except Exception:
                        continue
                    if isinstance(s, str):
                        out.append((tok.start[0], s, False))
                elif tok.type == tokenize.COMMENT:
                    for m in QUOTED.finditer(tok.string):
                        out.append((tok.start[0], m.group('s'), True))
            return out
        except (tokenize.TokenError, IndentationError, SyntaxError):
            out = []
    for i, line in enumerate(text.splitlines(), 1):
        commented = line.lstrip().startswith('#')
        for m in QUOTED.finditer(line):
            out.append((i, m.group('s'), commented))
    return out


def scan(head, ids):
    files = []
    for g in SOURCE_GLOBS:
        files += sorted(glob.glob(os.path.join(head, g)))
    refs, patterns, sheets, inputs, spans_by_file, lines_by_file = {}, {}, {}, {}, {}, {}
    for path in files:
        rel = os.path.relpath(path, head)
        text = open(path, encoding='utf-8', errors='replace').read()
        lines = text.splitlines()
        spans = _function_spans(path, text)
        spans_by_file[rel], lines_by_file[rel] = spans, lines

        def ref(line, commented):
            return {'f': rel, 'l': line, 'fn': _enclosing(spans, line), 't': lines[line - 1].strip()[:170], 'c': commented}
        seen = set()
        for line, s, commented in _literals(path, text):
            key = (line, s)
            if key in seen:
                continue
            seen.add(key)
            if s in ids:
                refs.setdefault(s, []).append(ref(line, commented))
            elif s in PREFIXES:
                patterns.setdefault(s, []).append(ref(line, commented))
        for i, l in enumerate(lines, 1):
            commented = l.lstrip().startswith('#')
            for m in SHEET.finditer(l):
                sheets.setdefault(m.group(1), []).append(ref(i, commented))
            for m in INPUT.finditer(l):
                inputs.setdefault(m.group(1), []).append(ref(i, commented))
    return {'refs': refs, 'patterns': patterns, 'sheets': sheets, 'inputs': inputs}, spans_by_file, lines_by_file


def resolve_anchor(anchor, spans_by_file):
    rel, fn = anchor.split('::')
    if rel not in spans_by_file:
        raise RuntimeError('roles.py anchor %s: file not found' % anchor)
    hits = [(n, a, b) for n, a, b in spans_by_file[rel] if n == fn or n.endswith('.' + fn)]
    if not hits:
        raise RuntimeError('roles.py anchor %s: no function %s in %s (renamed or removed?)' % (anchor, fn, rel))
    n, a, b = hits[0]
    return {'f': rel, 'fn': n, 'l': a, 'end': b}


def input_table(head, code):
    """Every sheet of every workbook, with where the code reads it."""
    import pandas as pd
    rows = []
    for path in sorted(glob.glob(os.path.join(head, 'input_data', '*'))):
        name = os.path.basename(path)
        if os.path.isdir(path):
            continue
        file_reads = code['inputs'].get(name, [])
        if name.endswith('.xlsx'):
            for sheet in pd.ExcelFile(path).sheet_names:
                reads = [r for r in code['sheets'].get(sheet, []) if not r['c']]
                n = len(pd.read_excel(path, sheet_name=sheet))
                rows.append({'file': name, 'sheet': sheet, 'rows': n, 'reads': code['sheets'].get(sheet, []),
                             'used': bool(reads) and any(not r['c'] for r in file_reads)})
        else:
            rows.append({'file': name, 'sheet': None, 'reads': file_reads, 'used': any(not r['c'] for r in file_reads),
                         'size': os.path.getsize(path)})
    return rows


def attach(spec, head):
    ids = set(spec['species']) | set(spec['reactions'])
    code, spans, lines = scan(head, ids)
    anchors = {}
    for p in spec['processes'].values():
        p['anchor_refs'] = []
        for a in p['anchors']:
            anchors[a] = resolve_anchor(a, spans)
            p['anchor_refs'].append(anchors[a])
            info = anchors[a]
            named = p.setdefault('code_species', [])
            for sid, rs in code['refs'].items():            # species named inside the anchor functions (dev view only)
                if sid in spec['species'] and sid not in named and any(
                        r['f'] == info['f'] and info['l'] <= r['l'] <= info['end'] and not r['c'] for r in rs):
                    named.append(sid)
    code['anchors'] = anchors
    code['prefixes'] = PREFIXES
    code['prefix_kinds'] = PREFIX_KINDS
    spec['code'] = code
    spec['inputs'] = input_table(head, code)
    from .load import link_roles
    link_roles(spec)
