"""Inline template.html + style.css + app.js + KaTeX + the spec into one self-contained page (no network needed)."""

import base64
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def katex_inline():
    """vendor/katex (KaTeX 0.16.11, MIT): the stylesheet with its woff2 fonts as data URIs, and the script."""
    d = os.path.join(HERE, 'vendor', 'katex')
    css = open(os.path.join(d, 'katex.min.css'), encoding='utf-8').read()

    def font(m):
        with open(os.path.join(d, m.group(1)), 'rb') as f:
            return 'url(data:font/woff2;base64,%s) format("woff2")' % base64.b64encode(f.read()).decode()
    # keep only the woff2 source of each @font-face (woff/ttf are not vendored)
    css = re.sub(r'url\((fonts/[^)]+\.woff2)\) format\("woff2"\)(,url\([^)]+\) format\("[a-z]+"\))*', font, css)
    js = open(os.path.join(d, 'katex.min.js'), encoding='utf-8').read()
    return '<style>\n%s\n</style>\n<script>\n%s\n</script>' % (css, js.replace('</script', '<\\/script'))


def write_html(spec, path):
    def read(name):
        return open(os.path.join(HERE, name), encoding='utf-8').read()
    data = json.dumps(spec, separators=(',', ':')).replace('</', '<\\/')
    html = (read('template.html').replace('<!--KATEX-->', katex_inline()).replace('/*STYLE*/', read('style.css'))
            .replace('/*APP*/', read('app.js')).replace('/*SPEC*/', data))
    with open(path, 'w', encoding='utf-8') as f:
        f.write(html)
    return path


def write_artifact(html_path, path):
    """The same page for a claude.ai artifact: the viewer wraps it in its own doctype/head/body with charset and viewport,
    so drop the document wrapper and the meta tags; the <title> stays first."""
    html = open(html_path, encoding='utf-8').read()
    for tag in ('<!DOCTYPE html>', '<html lang="en">', '<head>', '</head>', '</html>', '<meta charset="utf-8">',
                '<meta name="viewport" content="width=device-width, initial-scale=1">'):
        html = html.replace(tag, '', 1)
    html = html.replace('<body>', '', 1)
    html = html[::-1].replace('>ydob/<', '', 1)[::-1]        # the last </body>
    with open(path, 'w', encoding='utf-8') as f:
        f.write(html.lstrip())
    return path
