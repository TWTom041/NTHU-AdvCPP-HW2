// Render excerpts of the captured console logs (results/*.txt) as terminal
// screenshots for the report:  node scripts/screenshots.js
//
// report/shots.json: [{ "name": "q3", "file": "q3_fastmath.txt", "from": "$ build/q3_strict", "lines": 20 }, ...]
// "from" is the first line to include (substring match), "lines" how many lines.
const fs = require('fs');
const path = require('path');
const { chromium } = require(path.join(require('child_process').execSync('npm root -g').toString().trim(), 'playwright'));

const root = path.resolve(__dirname, '..');
const outDir = path.join(root, 'report', 'shots');
fs.mkdirSync(outDir, { recursive: true });
const shots = JSON.parse(fs.readFileSync(path.join(root, 'report', 'shots.json'), 'utf8'));

const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

function excerpt(shot) {
  const lines = fs.readFileSync(path.join(root, 'results', shot.file), 'utf8').split('\n');
  let start = 0;
  if (shot.from) {
    start = lines.findIndex((l) => l.includes(shot.from));
    if (start < 0) throw new Error(`${shot.name}: "${shot.from}" not found in ${shot.file}`);
  }
  const body = lines.slice(start, start + (shot.lines || 40));
  while (body.length && body[body.length - 1].trim() === '') body.pop();
  return body;
}

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 2 });
  for (const shot of shots) {
    const body = excerpt(shot)
      .map((l) => (l.startsWith('$ ') ? `<span class="p">user@hw2:~/NTHU-AdvCPP-HW2</span>${esc(l)}` : esc(l)))
      .join('\n');
    await page.setContent(`<!doctype html><html><head><meta charset="utf-8"><style>
      body { margin: 0; background: #fff; }
      .win { display: inline-block; margin: 8px; border-radius: 8px; overflow: hidden;
             box-shadow: 0 2px 10px rgba(0,0,0,.25); background: #1e1f22; }
      .bar { height: 26px; background: #2d2f33; display: flex; align-items: center; padding-left: 10px; gap: 7px;
             color: #aaa; font: 12px sans-serif; }
      .dot { width: 11px; height: 11px; border-radius: 50%; }
      pre { margin: 0; padding: 10px 14px; color: #e6e6e6; font: 12.5px/1.35 "DejaVu Sans Mono", monospace; }
      .p { color: #7ec16e; }
    </style></head><body><div class="win"><div class="bar">
      <span class="dot" style="background:#ff5f57"></span><span class="dot" style="background:#febc2e"></span>
      <span class="dot" style="background:#28c840"></span><span style="margin-left:10px">bash — ${esc(shot.file)}</span>
    </div><pre>${body}</pre></div></body></html>`);
    await page.locator('.win').screenshot({ path: path.join(outDir, shot.name + '.png') });
    console.log('shot', shot.name);
  }
  await browser.close();
})();
