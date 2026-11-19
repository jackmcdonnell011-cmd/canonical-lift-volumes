"""Build volume_explorer.html: volumes of the canonical lift complements on the modular surface and the
once-punctured torus, plotted against word length or hyperbolic length.

usage:  python make_explorer.py modular_geodesics_L25.csv.gz volume_explorer.html torus11.json torus11_twists.json
"""
import json, sys
import numpy as np, pandas as pd

csv, out = sys.argv[1], sys.argv[2]
FULL_UP_TO, SAMPLE_PER_LENGTH = 16, 3000

df = pd.read_csv(csv, usecols=["word", "hyperbolic_length", "volume"]).rename(columns={"hyperbolic_length": "geodesic_length"})
df["L"] = df.word.str.len()

def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))

keep = set(df.loc[df.L <= FULL_UP_TO, "word"])
rng = np.random.default_rng(0)
for L, grp in df[df.L > FULL_UP_TO].groupby("L"):
    keep.update(grp.word.values[rng.choice(len(grp), SAMPLE_PER_LENGTH, replace=False)])
for _, grp in df.groupby("L"):
    keep.add(grp.loc[grp.volume.idxmin(), "word"])
    keep.add(grp.loc[grp.volume.idxmax(), "word"])
families = {
    "X^kY": [canon("X" * k + "Y") for k in range(1, 25)],
    "X^kY^2": [canon("X" * k + "YY") for k in range(1, 24)],
    "X^kY^3": [canon("X" * k + "YYY") for k in range(1, 23)],
}
for ws in families.values():
    keep.update(ws)

sub = df[df.word.isin(keep)].sort_values(["L", "word"])
counts = df.groupby("L").size()
surface = {
    "label": "Modular surface",
    "maxL": int(df.L.max()),
    "fullUpTo": FULL_UP_TO,
    "counts": {int(k): int(v) for k, v in counts.items()},
    "w": sub.word.tolist(),
    "len": [round(x, 3) for x in sub.geodesic_length],
    "vol": [round(x, 3) for x in sub.volume],
    "families": {k: [w for w in ws if w in keep] for k, ws in families.items()},
}
data = {"mod": surface}
# the once-punctured torus: JSON files {"1,1": {label, maxL, fullUpTo, counts, w, len, vol, ...}} merged per
# surface, so the second file adds the key "twists" ({name: {k, w, code, wl, len, vol}})
for extra in sys.argv[3:]:
    for key, val in json.load(open(extra, encoding="utf-8")).items():
        data.setdefault(key, {}).update(val)

html = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Canonical Lift Volumes</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Lato:wght@400;700&family=Playfair+Display:wght@600&display=swap" rel="stylesheet">
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
<style>
  :root { color-scheme: light; --title:#4e4e4e; --text:#5b5b5b; --light:#858585; --line:#dadada; --grid:#dddddd;
          --wash:#f2f2f2; --hover:#f7f7f7; --blue:#3a5dae; --dot:#4C6A9C; }
  * { box-sizing: border-box; }
  html, body { margin:0; background:#f7f7f7; color:var(--text); font:14px/1.4 Lato, "Helvetica Neue", Arial, sans-serif; }
  .page { padding:32px 16px 48px; }
  .grapher { max-width:1000px; margin:0 auto; background:#fff; padding:20px 20px 14px; display:flex; flex-direction:column;
             box-shadow:0 0 0 1px rgba(0,0,0,.04), 0 2px 2px rgba(0,0,0,.08); }
  .grapher:fullscreen { max-width:none; height:100vh; overflow:auto; }
  .g-head h1 { font:600 24px/1.2 "Playfair Display", Georgia, serif; color:var(--title); margin:0; letter-spacing:-.01em; }
  .g-head p { margin:6px 0 0; font-size:15px; line-height:1.45; color:var(--text); max-width:90ch; }

  .g-controls { display:flex; flex-wrap:wrap; gap:8px; align-items:center; justify-content:space-between; margin:16px 0 10px; }
  .tabs { display:inline-flex; background:var(--wash); border-radius:4px; padding:2px; }
  .tabs button { display:inline-flex; align-items:center; gap:6px; height:28px; padding:0 12px; border:0; border-radius:3px;
                 background:none; color:var(--text); font:13px Lato, sans-serif; cursor:pointer; }
  .tabs button[aria-selected=true] { background:#fff; color:var(--title); font-weight:700; box-shadow:0 0 0 1px rgba(0,0,0,.06), 0 1px 2px rgba(0,0,0,.1); }
  .right { display:flex; flex-wrap:wrap; gap:8px; align-items:center; }
  .btn, .pick { display:inline-flex; align-items:center; gap:6px; height:32px; padding:0 10px; border:1px solid var(--line);
                border-radius:4px; background:#fff; color:var(--text); font:13px Lato, sans-serif; cursor:pointer; white-space:nowrap; }
  .btn:hover, .pick:hover { background:var(--hover); }
  .pick select, .pick input { font:inherit; color:var(--title); font-weight:700; border:0; background:none; padding:0; outline:none; }
  .pick select { cursor:pointer; }
  .pick input { width:40px; }
  svg.i { width:14px; height:14px; flex:none; fill:none; stroke:currentColor; stroke-width:1.6; stroke-linecap:round; stroke-linejoin:round; }
  details.settings { position:relative; }
  details.settings summary { list-style:none; }
  details.settings summary::-webkit-details-marker { display:none; }
  .menu { position:absolute; right:0; top:38px; z-index:20; width:260px; background:#fff; border-radius:4px; padding:14px 16px 10px;
          box-shadow:0 0 0 1px rgba(0,0,0,.06), 0 4px 16px rgba(0,0,0,.14); }
  .menu h4 { margin:0 0 10px; font-size:13px; color:var(--title); }
  .menu .sec { margin:10px 0 6px; font-size:11px; font-weight:700; letter-spacing:.06em; text-transform:uppercase; color:var(--light); }
  .menu label { display:flex; align-items:center; gap:8px; padding:4px 0; font-size:13px; cursor:pointer; }
  .menu input { accent-color:var(--blue); margin:0; }

  .legend { display:flex; flex-wrap:wrap; gap:4px 16px; font-size:13px; color:var(--text); margin:0 0 2px; min-height:18px; }
  .legend i { display:inline-block; width:9px; height:9px; border-radius:50%; background:var(--dot); opacity:.75; margin-right:6px; vertical-align:-1px; }
  #plot { width:100%; height:480px; }
  .grapher:fullscreen #plot { height:calc(100vh - 290px); }
  .empty { height:480px; display:flex; flex-direction:column; align-items:center; justify-content:center; text-align:center; color:var(--light); gap:6px; }
  .empty b { font-size:16px; color:var(--title); }

  .tbl { max-height:480px; overflow:auto; border-top:1px solid var(--line); }
  .grapher:fullscreen .tbl { max-height:calc(100vh - 290px); }
  table { width:100%; border-collapse:collapse; font-size:13px; font-variant-numeric:tabular-nums; }
  th { position:sticky; top:0; background:#fff; text-align:left; font-weight:700; color:var(--title); padding:8px 10px;
       border-bottom:1px solid #b9b9b9; cursor:pointer; white-space:nowrap; user-select:none; }
  th .arr { color:var(--light); font-size:10px; margin-left:4px; }
  td { padding:6px 10px; border-bottom:1px solid #ececec; color:var(--text); }
  th.num, td.num { text-align:right; }
  td.w { font-family:"SFMono-Regular", Consolas, monospace; font-size:12px; color:var(--title); }
  tbody tr:hover td { background:var(--hover); }
  .tcap { font-size:12px; color:var(--light); margin:8px 0 0; }

  .timeline { display:flex; align-items:center; gap:10px; margin:12px 0 4px; }
  .play { width:32px; height:32px; border-radius:50%; border:1px solid var(--line); background:#fff; color:var(--title);
          display:inline-flex; align-items:center; justify-content:center; cursor:pointer; flex:none; }
  .play:hover { background:var(--hover); }
  .play svg { width:12px; height:12px; fill:currentColor; }
  .tl-lab { font-size:13px; color:var(--title); border:1px solid var(--line); border-radius:3px; padding:3px 7px; white-space:nowrap;
            font-variant-numeric:tabular-nums; }
  input[type=range] { flex:1; min-width:80px; accent-color:var(--blue); height:4px; }

  .g-foot { display:flex; flex-wrap:wrap; gap:10px 24px; justify-content:space-between; align-items:flex-end;
            border-top:1px solid #e7e7e7; margin-top:10px; padding-top:10px; }
  .g-foot .txt { flex:1 1 420px; font-size:13px; color:var(--text); }
  .g-foot p { margin:0 0 4px; }
  .g-foot .small { font-size:12px; color:var(--light); }
  .actions { display:flex; gap:8px; }
  button:focus-visible, .pick:focus-within, summary:focus-visible, input:focus-visible { outline:2px solid var(--blue); outline-offset:1px; }
  [hidden] { display:none !important; }
  @media (max-width:600px) { .page { padding:0; } .grapher { padding:14px 12px; } .g-head h1 { font-size:20px; }
                             .btn .lbl { display:none; } #plot, .empty { height:380px; } }
</style>
</head>
<body>
<div class="page">
<figure class="grapher" id="grapher" style="margin:0 auto">
  <header class="g-head">
    <h1 id="title"></h1>
    <p id="subtitle"></p>
  </header>

  <div class="g-controls">
    <div class="tabs" role="tablist">
      <button role="tab" id="tab-chart" aria-selected="true"><svg class="i" viewBox="0 0 16 16"><path d="M1.5 14.5h13M3 11l3.5-4 3 2.5L14 3.5"/></svg>Chart</button>
      <button role="tab" id="tab-table" aria-selected="false"><svg class="i" viewBox="0 0 16 16"><rect x="1.5" y="2.5" width="13" height="11" rx="1"/><path d="M1.5 6.5h13M1.5 10h13M6 6.5v7"/></svg>Table</button>
    </div>
    <div class="right">
      <label class="pick">Surface
        <select id="surface">
          <option value="mod">Modular surface</option>
          <option value="1,1">Once-punctured torus</option>
        </select>
      </label>
      <details class="settings" id="settings">
        <summary class="btn"><svg class="i" viewBox="0 0 16 16"><path d="M2 4h12M2 8h12M2 12h12"/><circle cx="5" cy="4" r="1.6" fill="#fff"/><circle cx="11" cy="8" r="1.6" fill="#fff"/><circle cx="7" cy="12" r="1.6" fill="#fff"/></svg><span class="lbl">Settings</span></summary>
        <div class="menu">
          <h4>Chart settings</h4>
          <div class="sec">X axis</div>
          <label><input type="radio" name="x" value="word" checked> Word length</label>
          <label><input type="radio" name="x" value="len"> Hyperbolic length</label>
          <div class="sec">Show</div>
          <label><input type="checkbox" id="fam" checked> Twist families</label>
        </div>
      </details>
    </div>
  </div>

  <div id="body">
    <div class="legend" id="legend"></div>
    <div id="plot"></div>
    <div id="table" hidden></div>
    <div class="empty" id="empty" hidden></div>
  </div>

  <div class="timeline" id="timeline">
    <button class="play" id="play" aria-label="Play"><svg viewBox="0 0 12 12" id="playicon"><path d="M2.5 1.2v9.6L10.5 6z"/></svg></button>
    <span class="tl-lab" id="Lmin">2</span>
    <input type="range" id="L" min="2" step="1" aria-label="Maximum word length">
    <span class="tl-lab" id="Lval"></span>
  </div>

  <footer class="g-foot">
    <div class="txt">
      <p><b>Data source:</b> M. Banarse, T. Cremaschi and J. McDonnell, <i>Volumes of canonical knot complements in tangent bundles</i> (2026); volumes computed with SnapPy<span id="verified"></span>.</p>
      <p class="small" id="note"></p>
    </div>
    <div class="actions">
      <button class="btn" id="dl"><svg class="i" viewBox="0 0 16 16"><path d="M8 2v8M4.5 6.5 8 10l3.5-3.5M2 13.5h12"/></svg><span class="lbl">Download</span></button>
      <button class="btn" id="share"><svg class="i" viewBox="0 0 16 16"><path d="M6.5 9.5 9.5 6.5M7 4.5l1.5-1.5a2.5 2.5 0 0 1 3.5 3.5L10.5 8M9 11.5 7.5 13a2.5 2.5 0 0 1-3.5-3.5L5.5 8"/></svg><span class="lbl" id="sharelbl">Copy link</span></button>
      <button class="btn" id="fs" aria-label="Enter full-screen"><svg class="i" viewBox="0 0 16 16"><path d="M2 6V2h4M10 2h4v4M14 10v4h-4M6 14H2v-4"/></svg></button>
    </div>
  </footer>
</figure>
</div>
<script>
const DATA = __DATA__;
const FAM = ["#C05917", "#6D3E91", "#58AC8C", "#BC8E5A"];
const state = { s: "mod", L: 12, x: "word", tab: "chart", sort: { key: "vol", dir: 1 } };
const $ = id => document.getElementById(id);
const S = () => DATA[state.s];
const WL = (s, i) => s.wl ? s.wl[i] : s.w[i].length;
const nameOf = s => s.label.replace(/\s*\(.*\)\s*$/, "");
const viaOf = s => (s.label.match(/\((.*)\)\s*$/) || [null, ""])[1];
const pretty = n => n.replace("^k", "ᵏ").replace("^2", "²").replace("^3", "³");
const available = () => Object.keys(DATA).filter(k => k !== "mod").map(k => `(${k})`).join(", ");

function chooseSurface() {
  state.s = $("surface").value;
  setup(); render();
}

function setup() {
  const s = S();
  $("timeline").hidden = !s;
  if (!s) return;
  $("L").max = s.maxL; state.L = Math.max(2, Math.min(state.L, s.maxL)); $("L").value = state.L;
}

function render() {
  const s = S(), chart = state.tab === "chart";
  $("empty").hidden = !!s;
  $("legend").hidden = $("plot").hidden = !s || !chart;
  $("table").hidden = !s || chart;
  if (!s) {
    const [g, n] = state.s.split(",");
    $("title").textContent = `Volume vs. word length: genus ${g} with ${n} puncture${n === "1" ? "" : "s"}`;
    $("subtitle").textContent = "";
    $("empty").innerHTML = `<b>No data available</b><span>No curves have been computed on this surface yet.</span><span>Surfaces with data: ${available()}.</span>`;
    $("note").textContent = ""; $("verified").textContent = "";
    Plotly.purge("plot");
    return;
  }
  const L = state.L, byWord = state.x === "word", mod = state.s === "mod";
  const idx = [];
  for (let i = 0; i < s.w.length; i++) if (WL(s, i) <= L) idx.push(i);
  let total = 0;
  for (let k = 2; k <= L; k++) total += s.counts[k] || 0;
  let lo = idx[0];
  idx.forEach(i => { if (s.vol[i] < s.vol[lo]) lo = i; });

  $("Lval").textContent = `|w| ≤ ${L}`;
  $("title").textContent = `Volume vs. ${byWord ? "word length" : "hyperbolic length"}: ${mod ? "the modular surface" : nameOf(s).toLowerCase()}`;
  const via = viaOf(s);
  $("subtitle").textContent = mod
    ? `Hyperbolic volume of the complement of the trefoil and the canonical lift of each closed geodesic, for the ${total.toLocaleString()} geodesics with code word length up to ${L}. The smallest is ${s.vol[lo]} (${s.w[lo]}).`
    : `Hyperbolic volume of the canonical lift complement of each closed curve that fills the surface, for the ${total.toLocaleString()} curves with word length up to ${L} in a free basis of π₁(S). The smallest is ${s.vol[lo]} (${s.w[lo]}).${via ? ` The surface is uniformised ${via}.` : ""}`;
  $("verified").textContent = mod ? "" : ", with hyperbolicity checked by an exact test of whether each curve fills the surface and every volume verified with interval arithmetic in SageMath";
  const fam = $("fam").checked && (Object.keys(s.families || {}).length || Object.keys(s.twists || {}).length);
  $("note").textContent = (L > s.fullUpTo ? `Note: every code word up to length ${s.fullUpTo} is shown; above that a random sample, together with the smallest and largest volume at each length. `
    : mod ? "" : "Note: each curve is counted once up to orientation; curves that do not fill the surface have non-hyperbolic complements and are not shown. ")
    + (fam && s.twists && state.tab === "chart" ? "Lines are twist families u tᵏ, computed separately for larger k and drawn in full." : "");
  $("legend").innerHTML = `<span><i></i>${mod ? "Closed geodesic" : "Closed curve on S"}</span>` +
    (byWord ? `<span style="color:#858585">Larger dots are several curves with the same word length and volume</span>` : "");
  if (chart) drawChart(s, idx, byWord, L); else drawTable(s, idx);
}

function drawChart(s, idx, byWord, L) {
  const row = (head, w, code, wl, len, vol) => `<b>${head}</b><br>${code ? `Code word: ${code}<br>` : ""}Word length: ${wl}<br>Hyperbolic length: ${len}<br>Volume: <b>${vol}</b>`;
  // against word length, curves with the same length and volume are one dot, sized by how many there are
  const groups = new Map();
  idx.forEach(i => {
    const key = byWord ? `${WL(s, i)}|${s.vol[i]}` : `${i}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(i);
  });
  const G = [...groups.values()], base = idx.length < 300 ? 7 : 4;
  const many = g => `<b>${g.length} curves</b><br>Word length: ${WL(s, g[0])}<br>Volume: <b>${s.vol[g[0]]}</b><br>` +
    g.slice(0, 8).map(i => s.w[i] + (s.code ? ` (${s.code[i]})` : "")).join("<br>") + (g.length > 8 ? `<br>and ${g.length - 8} more` : "");
  const traces = [{
    type: "scattergl", mode: "markers", x: G.map(g => byWord ? WL(s, g[0]) : s.len[g[0]]), y: G.map(g => s.vol[g[0]]),
    text: G.map(g => g.length > 1 ? many(g) : row(s.w[g[0]], s.w[g[0]], s.code && s.code[g[0]], WL(s, g[0]), s.len[g[0]], s.vol[g[0]])),
    hovertemplate: "%{text}<extra></extra>",
    marker: { size: G.map(g => Math.min(base + 2.2 * Math.sqrt(g.length - 1), 22)), color: "#4C6A9C", opacity: 0.6,
              line: { color: "#fff", width: 0.5 } },
  }];
  const notes = [];
  if ($("fam").checked) {
    const lines = [];
    const pos = new Map(s.w.map((w, i) => [w, i]));
    Object.entries(s.families || {}).forEach(([name, ws]) => {
      const ii = ws.map(w => pos.get(w)).filter(i => i !== undefined && WL(s, i) <= L);
      if (ii.length) lines.push({ name: pretty(name), x: ii.map(i => byWord ? WL(s, i) : s.len[i]), y: ii.map(i => s.vol[i]),
                                  text: ii.map(i => row(`${pretty(name)}: ${s.w[i]}`, s.w[i], null, WL(s, i), s.len[i], s.vol[i])) });
    });
    // twist families were computed separately, so each is drawn in full, also past L and the computed range
    Object.entries(s.twists || {}).forEach(([name, f]) => {
      const l = { name, x: [], y: [], text: [] };
      f.w.forEach((w, i) => {
        if (l.x.length && f.k[i] !== f.k[i - 1] + 1) { l.x.push(null); l.y.push(null); l.text.push(""); }   // a member is missing
        l.x.push(byWord ? f.wl[i] : f.len[i]); l.y.push(f.vol[i]);
        l.text.push(row(`${name}, k = ${f.k[i]}`, w, `${f.code[i]}<br>Surface word: ${w}`, f.wl[i], f.len[i], f.vol[i]));
      });
      lines.push(l);
    });
    lines.forEach((l, j) => {
      const c = FAM[j % FAM.length];
      traces.push({ type: "scatter", mode: "lines+markers", x: l.x, y: l.y, text: l.text, hovertemplate: "%{text}<extra></extra>",
                    line: { color: c, width: 1.5 }, marker: { size: 5, color: c } });
      const last = l.x.length - 1;
      notes.push({ xref: "paper", x: 1, y: l.y[last], text: l.name, showarrow: false, xanchor: "left", xshift: 8,
                   font: { family: "Lato, sans-serif", size: 13, color: c } });
    });
  }
  const axis = (title, extra) => Object.assign({ title: { text: title, font: { size: 12, color: "#858585" } }, showgrid: true,
    gridcolor: "#dddddd", griddash: "3px,2px", zeroline: false, ticks: "", tickfont: { size: 12, color: "#858585" } }, extra);
  Plotly.react("plot", traces, {
    paper_bgcolor: "#fff", plot_bgcolor: "#fff", showlegend: false, hovermode: "closest", annotations: notes,
    font: { family: "Lato, sans-serif", size: 13, color: "#5b5b5b" },
    margin: { l: 44, r: notes.length ? 96 : 12, t: 6, b: 42 },
    hoverlabel: { bgcolor: "#fff", bordercolor: "#d0d0d0", align: "left", font: { family: "Lato, sans-serif", size: 13, color: "#4e4e4e" } },
    xaxis: axis(byWord ? (s.wlName || "word length |w|") : "hyperbolic length",
                Object.assign({ showline: true, linecolor: "#a1a1a1" }, byWord ? { dtick: Math.max(...traces.flatMap(t => t.x.filter(v => v !== null))) > 16 ? 2 : 1 } : {})),
    yaxis: axis("volume", { showline: false }),
  }, { responsive: true, displayModeBar: false });
}

function drawTable(s, idx) {
  const mod = state.s === "mod";
  const cols = [["w", mod ? "Code word" : "Surface word", 0], ...(s.code ? [["code", "Code word", 0]] : []),
                ["wl", "Word length", 1], ["len", "Hyperbolic length", 1], ["vol", "Volume", 1]];
  const val = (k, i) => k === "wl" ? WL(s, i) : s[k][i];
  const { key, dir } = state.sort;
  const rows = idx.slice().sort((a, b) => { const x = val(key, a), y = val(key, b); return (x < y ? -1 : x > y ? 1 : 0) * dir; });
  const shown = rows.slice(0, 1000);
  $("table").innerHTML = `<div class="tbl"><table><thead><tr>${cols.map(([k, t, n]) =>
      `<th data-k="${k}" class="${n ? "num" : ""}" aria-sort="${k === key ? (dir > 0 ? "ascending" : "descending") : "none"}">${t}<span class="arr">${k === key ? (dir > 0 ? "▲" : "▼") : ""}</span></th>`).join("")}</tr></thead>` +
    `<tbody>${shown.map(i => `<tr>${cols.map(([k, , n]) => `<td class="${n ? "num" : "w"}">${val(k, i)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>` +
    `<p class="tcap">${rows.length > shown.length ? `Showing the first ${shown.length.toLocaleString()} of ${rows.length.toLocaleString()} rows; click a column to sort.` : `${rows.length.toLocaleString()} rows; click a column to sort.`}</p>`;
  $("table").querySelectorAll("th").forEach(th => th.onclick = () => {
    const k = th.dataset.k;
    state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : 1 };
    render();
  });
}

let timer = null;
function stop() { clearInterval(timer); timer = null; $("playicon").innerHTML = '<path d="M2.5 1.2v9.6L10.5 6z"/>'; $("play").setAttribute("aria-label", "Play"); }
$("play").onclick = () => {
  if (timer) return stop();
  const max = +$("L").max;
  if (state.L >= max) { state.L = 2; }
  $("playicon").innerHTML = '<path d="M2.5 1.5h2.5v9H2.5zM7 1.5h2.5v9H7z"/>'; $("play").setAttribute("aria-label", "Pause");
  render();
  timer = setInterval(() => { if (state.L >= max) return stop(); state.L++; $("L").value = state.L; render(); }, 700);
};
$("L").addEventListener("input", () => { stop(); state.L = +$("L").value; render(); });
document.querySelectorAll('input[name=x]').forEach(r => r.onchange = () => { state.x = r.value; render(); });
$("fam").onchange = render;
$("surface").onchange = chooseSurface;
const tab = t => { state.tab = t; $("tab-chart").setAttribute("aria-selected", t === "chart"); $("tab-table").setAttribute("aria-selected", t === "table"); render(); };
$("tab-chart").onclick = () => tab("chart");
$("tab-table").onclick = () => tab("table");
document.addEventListener("click", e => { if (!$("settings").contains(e.target)) $("settings").open = false; });
$("dl").onclick = () => {
  const s = S();
  if (!s) return;
  const name = `volumes_${state.s.replace(",", "_")}_L${state.L}`;
  if (state.tab === "chart")
    return Plotly.downloadImage("plot", { format: "png", width: 1200, height: 700, scale: 2, filename: name });
  const head = ["word", ...(s.code ? ["code_word"] : []), "word_length", "hyperbolic_length", "volume"];
  const lines = [head.join(",")];
  for (let i = 0; i < s.w.length; i++)
    if (WL(s, i) <= state.L) lines.push([s.w[i], ...(s.code ? [s.code[i]] : []), WL(s, i), s.len[i], s.vol[i]].join(","));
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/csv" }));
  a.download = name + ".csv"; a.click(); URL.revokeObjectURL(a.href);
};
$("share").onclick = () => {
  const h = new URLSearchParams();
  if (state.s !== "mod") h.set("s", "torus");
  h.set("L", state.L); if (state.x !== "word") h.set("x", state.x); if (state.tab !== "chart") h.set("tab", state.tab);
  const url = location.href.split("#")[0] + "#" + h.toString();
  history.replaceState(null, "", url);
  const done = () => { $("sharelbl").textContent = "Link copied"; setTimeout(() => $("sharelbl").textContent = "Copy link", 1600); };
  navigator.clipboard ? navigator.clipboard.writeText(url).then(done, done) : done();
};
$("fs").onclick = () => document.fullscreenElement ? document.exitFullscreen() : $("grapher").requestFullscreen();
document.addEventListener("fullscreenchange", () => S() && state.tab === "chart" && Plotly.Plots.resize("plot"));

const hash = new URLSearchParams(location.hash.slice(1));
if (hash.get("s") === "torus") $("surface").value = "1,1";
if (hash.has("L")) state.L = +hash.get("L");
if (hash.get("x") === "len") { state.x = "len"; document.querySelector('input[name=x][value=len]').checked = true; }
if (hash.get("tab") === "table") tab("table");
chooseSurface();
</script>
</body>
</html>
"""
html = html.replace("__DATA__", json.dumps(data, separators=(",", ":")))
open(out, "w", encoding="utf-8").write(html)
print(len(sub), "points;", round(len(html) / 1e6, 2), "MB")
