#!/usr/bin/env python3
"""report/report.md -> report/report.html -> report/report.pdf (headless Chromium via Playwright).

  python3 scripts/build_report.py [output.pdf]
"""
import os
import subprocess
import sys

import markdown

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(ROOT, "report")

CSS = """
@page { size: A4; margin: 16mm 14mm 16mm 14mm; }
body { font-family: "WenQuanYi Zen Hei", "DejaVu Sans", sans-serif; font-size: 10.2pt; line-height: 1.55;
       color: #111; background: #fff; }
h1 { font-size: 19pt; margin: 0 0 4pt; }
h2 { font-size: 14.5pt; border-bottom: 1.5px solid #2a78d6; padding-bottom: 2pt; margin-top: 18pt;
     break-after: avoid; }
h3 { font-size: 11.5pt; margin-top: 12pt; break-after: avoid; color: #1d4f8f; }
h4 { font-size: 10.5pt; margin-top: 10pt; break-after: avoid; }
p, li { text-align: justify; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 8.8pt; background: #f2f1ee; padding: 0 2px;
       border-radius: 2px; }
pre { background: #f6f5f2; border: 1px solid #e4e3df; padding: 6pt 8pt; font-size: 8pt; line-height: 1.35;
      overflow-x: hidden; white-space: pre-wrap; break-inside: avoid; }
pre code { background: none; padding: 0; font-size: 8pt; }
table { border-collapse: collapse; margin: 6pt 0; font-size: 8.8pt; break-inside: avoid; }
th, td { border: 1px solid #d6d5d0; padding: 2pt 5pt; }
th { background: #eef3fa; }
img { max-width: 100%; display: block; margin: 6pt auto; break-inside: avoid; }
img.shot { max-width: 100%; }
.title { text-align: center; margin: 30pt 0 18pt; }
.title .t1 { font-size: 20pt; font-weight: bold; }
.title .t2 { font-size: 14pt; margin-top: 6pt; }
.caption { text-align: center; font-size: 8.5pt; color: #52514e; margin-top: -2pt; margin-bottom: 8pt; }
blockquote { border-left: 3px solid #2a78d6; margin: 6pt 0; padding: 2pt 10pt; color: #333; background: #f7f9fc; }
.pagebreak { break-after: page; }
"""


def main():
    out_pdf = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPORT, "report.pdf")
    with open(os.path.join(REPORT, "report.md"), encoding="utf-8") as f:
        md = f.read()
    body = markdown.markdown(md, extensions=["tables", "fenced_code", "toc", "attr_list", "md_in_html"])
    html = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<title>Advanced C++ Homework 2</title><style>{CSS}</style></head><body>{body}</body></html>"""
    html_path = os.path.join(REPORT, "report.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    npm_root = subprocess.check_output(["npm", "root", "-g"], text=True).strip()
    js = f"""
const {{ chromium }} = require('{npm_root}/playwright');
(async () => {{
  const b = await chromium.launch();
  const p = await b.newPage();
  await p.goto('file://{html_path}', {{ waitUntil: 'networkidle' }});
  await p.pdf({{ path: '{out_pdf}', format: 'A4', printBackground: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="font-size:8px;width:100%;text-align:center;color:#777"><span class="pageNumber"></span> / <span class="totalPages"></span></div>',
    margin: {{ top: '14mm', bottom: '16mm', left: '14mm', right: '14mm' }} }});
  await b.close();
}})();
"""
    subprocess.run(["node", "-e", js], check=True)
    print("wrote", html_path, "and", out_pdf)


if __name__ == "__main__":
    main()
