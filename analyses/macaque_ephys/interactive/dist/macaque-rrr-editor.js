(() => {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const STORAGE_KEY = "macaque-e4-t-e-rrr-reference-layout-v2";
  const DATA_VERSION = "macaque-e4-margins-20260922";
  const DATASETS = {
    gc_merged: {label: "Macaque E1–E4 · Ca + Pu + NAc", stem: "macaque_e4", classLabel: "E class"}
  };
  const MERGED_GROUPS = ["E1", "E2", "E3", "E4"];
  const UNMERGED_GROUPS = ["GC1", "GC2", "GC3", "GC4", "GC5", "GC6", "GC7", "GC8", "GC9", "GC10", "GC11"];
  let CLASS_GROUPS = [...MERGED_GROUPS];
  const D_GROUPS = ["D1", "D2"];
  const T_POINT_GROUPS = ["D1", "D2"];
  const DEFAULT_COLORS = {
    D1: "#D95F02", D2: "#008F7A",
    E1: "#F8766D", E2: "#7CAE00", E3: "#00BFC4", E4: "#C77CFF",
    GC1: "#4E79A7", GC2: "#F28E2B", GC3: "#59A14F", GC4: "#E15759", GC5: "#B07AA1",
    GC6: "#76B7B2", GC7: "#EDC948", GC8: "#FF9DA7", GC9: "#9C755F", GC10: "#BAB0AC", GC11: "#2F4B7C"
  };
  const E_FEATURE_ABBREVIATIONS = {};
  const PANEL_ORDER = ["T12", "T13", "E12", "E13"];
  const PANEL_META = {
    T12: {space: "T", pair: [1, 2], side: -1, count: "gene"},
    T13: {space: "T", pair: [1, 3], side: -1, count: "gene"},
    E12: {space: "E", pair: [1, 2], side: 1, count: "feature"},
    E13: {space: "E", pair: [1, 3], side: 1, count: "feature"}
  };
  const DEFAULT_PAIRS = Object.fromEntries(PANEL_ORDER.map(id => [id, [...PANEL_META[id].pair]]));
  const DEFAULT_STATE = {
    datasetKey: "gc_merged",
    geneCount: 9, featureCount: 10, fontSize: 11, labelAlign: "auto", boxOpacity: 0.9,
    cellSize: 2, cellOpacity: 0.72, panelLayout: "row", plotRadius: 90,
    panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,
    marginLeft: 0, marginRight: 0, marginTop: 0, marginBottom: 0,
    radialPercentile: 99, ellipseLevel: 0.9, loadingScale: 0.76,
    fontFamily: "Arial", arrowWidth: 1.4, connectorWidth: 0.8, ellipseWidth: 2, ellipseOpacity: 0.055,
    circleWidth: 2, axisWidth: 1, titleSize: 18, axisLabelSize: 13,
    axisColor: "#aeb6c2", axisTextColor: "#111111", showTicks: false,
    tickCount: 5, tickDecimals: 1, tickSize: 9,
    titleT: "Transcriptomic space", titleE: "Electrophysiological space", canvasColor: "#ffffff",
    exportWidth: 7.2, exportDpi: 600,
    showBoxes: true, showConnectors: true, showEllipses: true, showLegend: true, showCircle: true,
    showPoints: true,
    axisDisplay: Object.fromEntries(PANEL_ORDER.map(id => [id, {showX: true, showY: true, labelX: true, labelY: true, tickX: true, tickY: true}])),
    pairs: structuredClone(DEFAULT_PAIRS), positions: {}, aliases: {}, colors: {...DEFAULT_COLORS}
  };
  function mergeLayout(incoming = {}) {
    return {
      ...structuredClone(DEFAULT_STATE),
      ...incoming,
      colors: {...DEFAULT_COLORS, ...(incoming.colors || {})},
      pairs: {...structuredClone(DEFAULT_PAIRS), ...(incoming.pairs || {})},
      axisDisplay: Object.fromEntries(PANEL_ORDER.map(id => [id, {...DEFAULT_STATE.axisDisplay[id], ...(incoming.axisDisplay?.[id] || {})}])),
      positions: incoming.positions || {},
      aliases: incoming.aliases || {}
    };
  }
  function savedDefault(baseLayout = {}) {
    try {
      const saved = JSON.parse(localStorage.getItem(STORAGE_KEY) || "null");
      return mergeLayout(saved || baseLayout);
    } catch {
      return mergeLayout(baseLayout);
    }
  }
  let state = structuredClone(DEFAULT_STATE);
  let builtInLayout = {};
  let cells = [], genes = [], features = [], currentLabels = [], drag = null;
  let datasets = {}, versionManifest = {};
  const $ = id => document.getElementById(id);

  function svgNode(name, attrs = {}, text = "") {
    const n = document.createElementNS(NS, name);
    Object.entries(attrs).forEach(([k, v]) => n.setAttribute(k, v));
    if (text !== "") n.textContent = text;
    return n;
  }

  function cleanCSV(value) {
    return String(value ?? "").replace(/^\uFEFF/, "").replace(/^"|"$/g, "").replace(/""/g, '"');
  }

  function parseCSV(text) {
    const lines = text.trim().split(/\r?\n/);
    const head = lines.shift().split(",").map(cleanCSV);
    return lines.map(line => {
      const cols = line.split(",");
      return Object.fromEntries(head.map((h, i) => {
        const value = cleanCSV(cols[i] ?? "");
        return [h, /^(T_Component|E_Component|Component)\d$/.test(h) ? +value : value];
      }));
    });
  }

  function applyDataset(key, shouldRender = true) {
    const next = datasets[key];
    if (!next) return;
    state.datasetKey = key;
    cells = next.cells;
    genes = next.genes;
    features = next.features;
    CLASS_GROUPS = key === "gc_unmerged"
      ? UNMERGED_GROUPS.filter(group => cells.some(cell => cell.E_class === group))
      : [...MERGED_GROUPS];
    if ($("datasetVersion")) $("datasetVersion").value = key;
    rebuildColorControls();
    if (shouldRender) render();
  }

  function percentile(values, p) {
    const s = [...values].filter(Number.isFinite).sort((a, b) => a - b);
    if (!s.length) return 1;
    const index = (s.length - 1) * p, base = Math.floor(index), rest = index - base;
    return s[base + 1] === undefined ? s[base] : s[base] + rest * (s[base + 1] - s[base]);
  }

  function panelGeometry(id) {
    const i = PANEL_ORDER.indexOf(id), grid = state.panelLayout === "grid";
    const col = grid ? i % 2 : i, row = grid ? Math.floor(i / 2) : 0;
    const x = state.marginLeft + col * (state.panelWidth + state.columnGap);
    const y = state.marginTop + row * (state.panelHeight + state.rowGap);
    return {...PANEL_META[id], id, pair: state.pairs[id], x, y,
      cx: x + state.panelWidth / 2, cy: y + state.panelHeight / 2 + 4};
  }

  function canvasGeometry() {
    const grid = state.panelLayout === "grid";
    return {
      width: state.marginLeft + (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap + state.marginRight,
      height: state.marginTop + (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0) + state.marginBottom
    };
  }

  function scaledScores(panel) {
    const xs = cells.map(c => +c[`${panel.space}_Component${panel.pair[0]}`]);
    const ys = cells.map(c => +c[`${panel.space}_Component${panel.pair[1]}`]);
    const radii = xs.map((x, i) => Math.hypot(x, ys[i]));
    const scale = percentile(radii, state.radialPercentile / 100) || 1;
    return cells.map((c, i) => ({...c, sx: xs[i] / scale, sy: ys[i] / scale}));
  }

  function spread(values, gap = 0.18) {
    if (!values.length) return [];
    const indexed = values.map((v, i) => ({v: Math.max(-0.9, Math.min(0.9, v)), i})).sort((a, b) => a.v - b.v);
    for (let i = 1; i < indexed.length; i++) indexed[i].v = Math.max(indexed[i].v, indexed[i - 1].v + gap);
    let shift = Math.max(...indexed.map(x => x.v)) - 0.92;
    if (shift > 0) indexed.forEach(x => x.v -= shift);
    for (let i = indexed.length - 2; i >= 0; i--) indexed[i].v = Math.min(indexed[i].v, indexed[i + 1].v - gap);
    shift = -0.92 - Math.min(...indexed.map(x => x.v));
    if (shift > 0) indexed.forEach(x => x.v += shift);
    const out = [];
    indexed.forEach(x => { out[x.i] = x.v; });
    return out;
  }

  function panelLoadings(panel) {
    const source = panel.space === "T" ? genes : features;
    const count = state[`${panel.count}Count`];
    const xKey = `Component${panel.pair[0]}`, yKey = `Component${panel.pair[1]}`;
    const top = source.map(d => ({...d, label: panel.space === "T" ? d.Gene : d.Display, vx: +d[xKey], vy: +d[yKey], mag: Math.hypot(+d[xKey], +d[yKey])}))
      .sort((a, b) => b.mag - a.mag).slice(0, count);
    const max = Math.max(1e-12, ...top.map(d => d.mag));
    top.forEach(d => { d.vx = d.vx * state.loadingScale / max; d.vy = d.vy * state.loadingScale / max; });
    const ys = spread(top.map(d => d.vy), count >= 10 ? 0.18 : 0.19);
    const defaultPair = panel.pair.join() === DEFAULT_PAIRS[panel.id].join();
    return top.map((d, i) => {
      const key = `${panel.id}${defaultPair ? "" : `-${panel.pair.join("")}`}:${d.label}`;
      const saved = state.positions[key];
      return {...d, key, lx: saved?.x ?? panel.side * 1.04, ly: saved?.y ?? ys[i]};
    });
  }

  function covarianceEllipse(points, group) {
    const q = points.filter(p => p.D1_D2 === group);
    if (q.length < 3) return [];
    const mx = q.reduce((s, p) => s + p.sx, 0) / q.length;
    const my = q.reduce((s, p) => s + p.sy, 0) / q.length;
    let a = 0, b = 0, c = 0;
    q.forEach(p => { const x = p.sx - mx, y = p.sy - my; a += x * x; b += x * y; c += y * y; });
    a /= q.length - 1; b /= q.length - 1; c /= q.length - 1;
    const tr = a + c, disc = Math.sqrt((a - c) ** 2 + 4 * b * b);
    const l1 = Math.max(0, (tr + disc) / 2), l2 = Math.max(0, (tr - disc) / 2);
    const theta = 0.5 * Math.atan2(2 * b, a - c);
    const k = Math.sqrt(-2 * Math.log(1 - state.ellipseLevel));
    return Array.from({length: 121}, (_, i) => {
      const t = i / 120 * Math.PI * 2;
      const x = k * Math.sqrt(l1) * Math.cos(t), y = k * Math.sqrt(l2) * Math.sin(t);
      return [mx + x * Math.cos(theta) - y * Math.sin(theta), my + x * Math.sin(theta) + y * Math.cos(theta)];
    });
  }

  function closestRectTouch(px, py, x, y, left, top, width, height) {
    const qx = px - x, qy = py - y, right = left + width, bottom = top + height;
    let tx = Math.max(left, Math.min(right, qx)), ty = Math.max(top, Math.min(bottom, qy));
    if (qx >= left && qx <= right && qy >= top && qy <= bottom) {
      const sides = [[Math.abs(qx-left),left,qy],[Math.abs(right-qx),right,qy],[Math.abs(qy-top),qx,top],[Math.abs(bottom-qy),qx,bottom]].sort((u,v)=>u[0]-v[0]);
      [, tx, ty] = sides[0];
    }
    return {x: x + tx, y: y + ty};
  }

  function textPlacement(align, left, width) {
    return align === "left" ? {x: left + 4, anchor: "start"} : align === "center" ? {x: left + width / 2, anchor: "middle"} : {x: left + width - 4, anchor: "end"};
  }

  function safe(value) { return value.replace(/[^A-Za-z0-9_-]/g, "_"); }

  function tickValues(count) {
    return Array.from({length: Math.max(2, count)}, (_, i) => -1 + 2 * i / (Math.max(2, count) - 1));
  }

  function tickLabel(value) {
    const clean = Math.abs(value) < 1e-10 ? 0 : value;
    return clean.toFixed(state.tickDecimals);
  }

  function render() {
    const svg = $("rrrSvg"), canvas = canvasGeometry();
    svg.innerHTML = ""; currentLabels = [];
    svg.setAttribute("viewBox", `0 0 ${canvas.width} ${canvas.height}`);
    svg.setAttribute("width", `${state.exportWidth}in`);
    svg.setAttribute("height", `${(state.exportWidth * canvas.height / canvas.width).toFixed(4)}in`);
    svg.style.fontFamily = state.fontFamily;
    svg.style.minWidth = state.panelLayout === "row" ? "1500px" : "900px";
    svg.style.background = state.canvasColor;
    const defs = svgNode("defs");
    const marker = svgNode("marker", {id: "rrrArrow", viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 5, markerHeight: 5, orient: "auto-start-reverse"});
    marker.append(svgNode("path", {d: "M 0 0 L 10 5 L 0 10 z", fill: "#111111"})); defs.append(marker); svg.append(defs);

    PANEL_ORDER.forEach(id => {
      const panel = panelGeometry(id), R = state.plotRadius, scores = scaledScores(panel);
      const X = v => panel.cx + v * R, Y = v => panel.cy - v * R;
      const g = svgNode("g", {"data-panel": id}); svg.append(g);
      const title = panel.space === "T" ? state.titleT : state.titleE;
      g.append(svgNode("text", {x: panel.cx, y: panel.y + 22, "text-anchor": "middle", "font-size": state.titleSize, fill: "#111111"}, title));
      if (state.showCircle) g.append(svgNode("circle", {cx: panel.cx, cy: panel.cy, r: R, fill: "none", stroke: "#333333", "stroke-width": state.circleWidth}));
      const axis = state.axisDisplay[id];
      if (axis.showX) g.append(svgNode("line", {"data-axis": "x", x1: X(-1.45), y1: panel.cy, x2: X(1.45), y2: panel.cy, stroke: state.axisColor, "stroke-width": state.axisWidth}));
      if (axis.showY) g.append(svgNode("line", {"data-axis": "y", x1: panel.cx, y1: Y(-1.2), x2: panel.cx, y2: Y(1.2), stroke: state.axisColor, "stroke-width": state.axisWidth}));
      if (state.showTicks) {
        tickValues(state.tickCount).forEach(value => {
          if (axis.showX && axis.tickX) {
            g.append(svgNode("line", {"data-axis-tick": "x", x1: X(value), y1: panel.cy - 3, x2: X(value), y2: panel.cy + 3, stroke: state.axisColor, "stroke-width": state.axisWidth}));
            g.append(svgNode("text", {"data-axis-tick-label": "x", x: X(value), y: panel.cy + state.tickSize + 6, "text-anchor": "middle", "font-size": state.tickSize, fill: state.axisTextColor}, tickLabel(value)));
          }
          if (axis.showY && axis.tickY) {
            g.append(svgNode("line", {"data-axis-tick": "y", x1: panel.cx - 3, y1: Y(value), x2: panel.cx + 3, y2: Y(value), stroke: state.axisColor, "stroke-width": state.axisWidth}));
            g.append(svgNode("text", {"data-axis-tick-label": "y", x: panel.cx - 6, y: Y(value) + state.tickSize * 0.34, "text-anchor": "end", "font-size": state.tickSize, fill: state.axisTextColor}, tickLabel(value)));
          }
        });
      }
      if (state.showEllipses) D_GROUPS.forEach(group => {
        const ellipse = covarianceEllipse(scores, group);
        if (ellipse.length) g.append(svgNode("polygon", {points: ellipse.map(([x,y]) => `${X(x)},${Y(y)}`).join(" "), fill: state.colors[group], "fill-opacity": state.ellipseOpacity, stroke: state.colors[group], "stroke-width": state.ellipseWidth}));
      });
      if (state.showPoints) scores.forEach(p => {
        const pointGroup = panel.space === "T" ? p.D1_D2 : p.E_class;
        g.append(svgNode("circle", {cx: X(p.sx), cy: Y(p.sy), r: state.cellSize, fill: state.colors[pointGroup], "fill-opacity": state.cellOpacity}));
      });

      if (state.showLegend) {
        const legendGroups = panel.space === "T" ? T_POINT_GROUPS : CLASS_GROUPS;
        const step = panel.space === "T" ? 38 : Math.max(30, Math.min(42, 210 / legendGroups.length)), total = step * legendGroups.length;
        const start = panel.cx - total / 2 + step / 2, y = panel.y + state.panelHeight - 36;
        legendGroups.forEach((group, index) => {
          const x = start + index * step;
          g.append(svgNode("circle", {cx: x, cy: y, r: 3, fill: state.colors[group]}));
          g.append(svgNode("text", {x: x + 6, y: y + 3, "font-size": Math.max(8, state.axisLabelSize - 2), fill: "#111111"}, group));
        });
      }

      panelLoadings(panel).forEach(d => {
        const defaultLabel = panel.space === "E" ? (E_FEATURE_ABBREVIATIONS[d.label] || d.label) : d.label;
        const display = state.aliases[d.key] ?? defaultLabel;
        currentLabels.push({panel: id, original: d.label, defaultLabel, label: display, key: d.key, x: d.lx, y: d.ly});
        const ex = X(d.vx), ey = Y(d.vy), lx = X(d.lx), ly = Y(d.ly);
        const width = Math.max(34, display.length * state.fontSize * 0.57 + 8), anchor = panel.side > 0 ? 0 : -width;
        const top = -state.fontSize * 0.85, height = state.fontSize * 1.22;
        const touch = closestRectTouch(ex, ey, lx, ly, anchor, top, width, height);
        const align = state.labelAlign === "auto" ? (panel.side > 0 ? "left" : "right") : state.labelAlign;
        const textPos = textPlacement(align, anchor, width);
        g.append(svgNode("line", {x1: panel.cx, y1: panel.cy, x2: ex, y2: ey, stroke: "#111111", "stroke-width": state.arrowWidth, "marker-end": "url(#rrrArrow)"}));
        g.append(svgNode("line", {id: `connector-${safe(d.key)}`, x1: ex, y1: ey, x2: touch.x, y2: touch.y, stroke: "#555555", "stroke-width": state.connectorWidth, "stroke-opacity": state.showConnectors ? 1 : 0}));
        const label = svgNode("g", {class: "drag-label", "data-key": d.key, "data-panel": id, "data-ex": ex, "data-ey": ey, "data-left": anchor, "data-top": top, "data-width": width, "data-height": height, transform: `translate(${lx} ${ly})`});
        label.append(svgNode("rect", {x: anchor, y: top, width, height, rx: 2, fill: state.canvasColor, "fill-opacity": state.showBoxes ? state.boxOpacity : 0, stroke: "#222222", "stroke-width": state.showBoxes ? 0.5 : 0}));
        label.append(svgNode("text", {x: textPos.x, y: 2, "text-anchor": textPos.anchor, "font-size": state.fontSize, fill: "#111111"}, display));
        g.append(label);
      });
      if (axis.labelX) g.append(svgNode("text", {"data-axis-title": "x", x: panel.cx, y: panel.y + state.panelHeight - 16, "text-anchor": "middle", "font-size": state.axisLabelSize, fill: state.axisTextColor}, `Component ${panel.pair[0]}`));
      if (axis.labelY) g.append(svgNode("text", {"data-axis-title": "y", x: panel.x + 18, y: panel.cy, "text-anchor": "middle", "font-size": state.axisLabelSize, fill: state.axisTextColor, transform: `rotate(-90 ${panel.x + 18} ${panel.cy})`}, `Component ${panel.pair[1]}`));
    });
    bindDrag();
    const d1 = cells.filter(cell => cell.D1_D2 === "D1").length;
    const d2 = cells.filter(cell => cell.D1_D2 === "D2").length;
    const manifest = versionManifest[state.datasetKey] || {};
    const agreement = manifest.agreement_n && manifest.agreement_denominator
      ? ` · GC–HC agreement ${manifest.agreement_n}/${manifest.agreement_denominator}` : "";
    const classText = state.datasetKey === "gc_merged" ? "E1–E4" : "raw GC subclusters";
    $("rrrSummary").textContent = `${cells.length} frozen D1/D2 cells · D1 ${d1} · D2 ${d2}${agreement} · ${classText} · rank 3 · ${state.panelLayout === "row" ? "1 x 4" : "2 x 2"} · ${state.exportWidth.toFixed(1)} in wide`;
    $("cohortHeader").textContent = `${DATASETS[state.datasetKey].label} · independently refitted T–E RRR · n = ${cells.length}`;
    $("versionNote").textContent = state.datasetKey === "gc_merged"
      ? "Frozen Macaque E1–E4 analysis: Ca + Pu + NAc HC–GC consensus cells with D1/D2 identity; Hybrid cells are excluded."
      : "Exploratory sensitivity analysis: raw GC K=11 matched to the best Ward.D2 HC K=13 solution; strict D1/D2 cells are shown by raw GC class.";
    $("rrrStatus").textContent = `Loaded ${DATASETS[state.datasetKey].label}. Drag labels to reposition them; loading vectors remain fixed to the frozen Macaque T–E RRR model.`;
  }

  function pointer(svg, event) {
    const p = svg.createSVGPoint(); p.x = event.clientX; p.y = event.clientY;
    return p.matrixTransform(svg.getScreenCTM().inverse());
  }

  function bindDrag() {
    document.querySelectorAll(".drag-label").forEach(label => {
      label.onpointerdown = event => { event.preventDefault(); label.setPointerCapture(event.pointerId); drag = {label, key: label.dataset.key, panel: panelGeometry(label.dataset.panel)}; };
      label.onpointermove = event => {
        if (!drag || drag.label !== label) return;
        const p = pointer($("rrrSvg"), event), x = (p.x - drag.panel.cx) / state.plotRadius, y = (drag.panel.cy - p.y) / state.plotRadius;
        state.positions[drag.key] = {x, y}; label.setAttribute("transform", `translate(${p.x} ${p.y})`);
        const current = currentLabels.find(d => d.key === drag.key); if (current) { current.x = x; current.y = y; }
        const line = $(`connector-${safe(drag.key)}`);
        const touch = closestRectTouch(+label.dataset.ex, +label.dataset.ey, p.x, p.y, +label.dataset.left, +label.dataset.top, +label.dataset.width, +label.dataset.height);
        line.setAttribute("x2", touch.x); line.setAttribute("y2", touch.y);
      };
      label.onpointerup = () => { drag = null; };
      label.ondblclick = event => {
        event.preventDefault(); const current = currentLabels.find(d => d.key === label.dataset.key);
        const next = prompt("Label text", current?.label ?? "");
        if (next !== null) { state.aliases[label.dataset.key] = next.trim() || current.defaultLabel; render(); }
      };
    });
  }

  const numeric = ["geneCount","featureCount","fontSize","boxOpacity","cellSize","cellOpacity","plotRadius","panelWidth","panelHeight","columnGap","rowGap","marginLeft","marginRight","marginTop","marginBottom","radialPercentile","ellipseLevel","loadingScale","arrowWidth","connectorWidth","ellipseWidth","ellipseOpacity","circleWidth","axisWidth","titleSize","axisLabelSize","tickCount","tickDecimals","tickSize","exportWidth","exportDpi"];
  const toggles = ["showBoxes","showConnectors","showEllipses","showLegend","showCircle","showPoints","showTicks"];
  const texts = ["titleT","titleE","canvasColor","axisColor","axisTextColor"];
  const decimals = new Set(["fontSize","boxOpacity","cellSize","cellOpacity","radialPercentile","ellipseLevel","loadingScale","arrowWidth","connectorWidth","ellipseWidth","ellipseOpacity","circleWidth","axisWidth"]);

  function outputValue(id) {
    const out = $(`${id}Value`); if (!out) return;
    const value = +state[id]; out.value = id === "ellipseLevel" ? `${Math.round(value * 100)}%` : decimals.has(id) ? value.toFixed(id === "ellipseOpacity" ? 3 : 2).replace(/0+$/, "").replace(/\.$/, "") : String(value);
  }

  function syncControls() {
    numeric.forEach(id => { if ($(id)) { $(id).value = state[id]; outputValue(id); } });
    toggles.forEach(id => { $(id).checked = state[id]; }); texts.forEach(id => { $(id).value = state[id]; });
    $("panelLayout").value = state.panelLayout; $("labelAlign").value = state.labelAlign; $("fontFamily").value = state.fontFamily;
    if ($("datasetVersion")) $("datasetVersion").value = state.datasetKey;
    Object.entries(state.pairs).forEach(([id, pair]) => { $(`${id}X`).value = pair[0]; $(`${id}Y`).value = pair[1]; });
    Object.entries(state.axisDisplay).forEach(([id, axis]) => {
      $(`${id}ShowX`).checked = axis.showX; $(`${id}ShowY`).checked = axis.showY;
      $(`${id}LabelX`).checked = axis.labelX; $(`${id}LabelY`).checked = axis.labelY;
      $(`${id}TickX`).checked = axis.tickX; $(`${id}TickY`).checked = axis.tickY;
    });
    Object.entries(state.colors).forEach(([id, value]) => { if ($(`color-${id}`)) $(`color-${id}`).value = value; });
  }

  function serializedSVG() {
    const copy = $("rrrSvg").cloneNode(true); copy.setAttribute("xmlns", NS); copy.style.minWidth = "";
    return new XMLSerializer().serializeToString(copy);
  }
  function download(name, blob) { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); }
  function exportSVG() { download(`Macaque_T-E_RRR_${state.datasetKey}_edited.svg`, new Blob([serializedSVG()], {type: "image/svg+xml"})); }
  function exportPNG() {
    const svg = $("rrrSvg"), vb = svg.viewBox.baseVal, width = Math.round(state.exportWidth * state.exportDpi), height = Math.round(width * vb.height / vb.width);
    if (width * height > 8e7) { $("rrrStatus").textContent = "Requested PNG exceeds the safe 80-megapixel browser limit."; return; }
    const img = new Image(), url = URL.createObjectURL(new Blob([serializedSVG()], {type: "image/svg+xml"}));
    img.onload = () => { const canvas = document.createElement("canvas"); canvas.width = width; canvas.height = height; const ctx = canvas.getContext("2d"); ctx.fillStyle = state.canvasColor; ctx.fillRect(0,0,width,height); ctx.drawImage(img,0,0,width,height); URL.revokeObjectURL(url); canvas.toBlob(blob => download(`Macaque_T-E_RRR_${state.datasetKey}_${state.exportDpi}dpi.png`, blob), "image/png"); };
    img.src = url;
  }

  function rebuildColorControls() {
    const palette = $("colorControls");
    if (!palette) return;
    palette.innerHTML = "";
    [...T_POINT_GROUPS, ...CLASS_GROUPS].forEach(id => {
      const label = document.createElement("label"), input = document.createElement("input");
      label.append(document.createTextNode(id));
      input.type = "color";
      input.id = `color-${id}`;
      input.value = state.colors[id] || DEFAULT_COLORS[id] || "#777777";
      input.oninput = event => { state.colors[id] = event.target.value; render(); };
      label.append(input);
      palette.append(label);
    });
  }

  function bindControls() {
    rebuildColorControls();
    numeric.forEach(id => { if (!$(id)) return; $(id).oninput = event => { state[id] = +event.target.value; outputValue(id); render(); }; });
    toggles.forEach(id => { $(id).onchange = event => { state[id] = event.target.checked; render(); }; });
    texts.forEach(id => { $(id).oninput = event => { state[id] = event.target.value; render(); }; });
    $("panelLayout").onchange = event => { state.panelLayout = event.target.value; render(); };
    $("labelAlign").onchange = event => { state.labelAlign = event.target.value; render(); };
    $("fontFamily").onchange = event => { state.fontFamily = event.target.value; render(); };
    $("datasetVersion").onchange = event => applyDataset(event.target.value);
    PANEL_ORDER.forEach(id => [0,1].forEach(axis => { const input = $(`${id}${axis ? "Y" : "X"}`); input.onchange = event => { const pair = state.pairs[id], next = +event.target.value; if (next === pair[1-axis]) { pair[1-axis] = pair[axis]; $(`${id}${axis ? "X" : "Y"}`).value = pair[1-axis]; } pair[axis] = next; render(); }; }));
    PANEL_ORDER.forEach(id => {
      [["ShowX", "showX"], ["ShowY", "showY"], ["LabelX", "labelX"], ["LabelY", "labelY"], ["TickX", "tickX"], ["TickY", "tickY"]].forEach(([suffix, key]) => {
        $(`${id}${suffix}`).onchange = event => { state.axisDisplay[id][key] = event.target.checked; render(); };
      });
    });
    $("saveDefault").onclick = () => {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
      $("rrrStatus").textContent = `Current ${DATASETS[state.datasetKey].label} parameters saved as the browser default.`;
    };
    $("restoreDefault").onclick = () => {
      state = savedDefault(builtInLayout);
      applyDataset(state.datasetKey in datasets ? state.datasetKey : "gc_merged", false);
      syncControls();
      render();
    };
    $("resetLabels").onclick = () => { state.positions = {}; state.aliases = {}; render(); };
    $("resetAll").onclick = () => { state = mergeLayout(builtInLayout); applyDataset(state.datasetKey in datasets ? state.datasetKey : "gc_merged", false); syncControls(); render(); };
    $("exportSvg").onclick = exportSVG; $("exportPng").onclick = exportPNG;
    $("exportCsv").onclick = () => { const rows = [["panel","original_label","display_label","x_normalized","y_normalized"], ...currentLabels.map(d => [d.panel,d.original,d.label,d.x,d.y])]; const csv = rows.map(row => row.map(v => `"${String(v).replaceAll('"','""')}"`).join(",")).join("\n"); download(`Macaque_T-E_RRR_${state.datasetKey}_label_positions.csv`, new Blob([csv], {type: "text/csv"})); };
    $("exportLayout").onclick = () => download(`Macaque_T-E_RRR_${state.datasetKey}_layout.json`, new Blob([JSON.stringify(state, null, 2)], {type: "application/json"}));
    $("importLayout").onchange = async event => { const file = event.target.files[0]; if (!file) return; try { state = mergeLayout(JSON.parse(await file.text())); applyDataset(state.datasetKey in datasets ? state.datasetKey : "gc_merged", false); syncControls(); render(); } catch { $("rrrStatus").textContent = "The selected layout JSON could not be read."; } event.target.value = ""; };
    syncControls();
  }

  try {
    const getText = (path, description) => fetch(`${path}?v=${DATA_VERSION}`).then(response => {
      if (!response.ok) throw new Error(`${description} is unavailable`);
      return response.text();
    });
    Promise.all([
      getText("rrr-data/rrr_scores_gc_merged.csv", "GC-merged cell scores"),
      getText("rrr-data/rrr_gene_loadings_gc_merged.csv", "GC-merged gene loadings"),
      getText("rrr-data/rrr_feature_loadings_gc_merged.csv", "GC-merged E-feature loadings"),
      getText("rrr-data/rrr_scores_gc_unmerged.csv", "GC-unmerged cell scores"),
      getText("rrr-data/rrr_gene_loadings_gc_unmerged.csv", "GC-unmerged gene loadings"),
      getText("rrr-data/rrr_feature_loadings_gc_unmerged.csv", "GC-unmerged E-feature loadings"),
      getText("rrr-data/rrr_gc_versions_manifest.csv", "GC-version manifest"),
      fetch(`rrr-data/Macaque_T-E_RRR_layout_reference.json?v=${DATA_VERSION}`).then(response => {
        if (!response.ok) throw new Error("default RRR layout is unavailable");
        return response.json();
      })
    ]).then(([mergedScores, mergedGenes, mergedFeatures, rawScores, rawGenes, rawFeatures, manifestText, defaultLayout]) => {
      builtInLayout = defaultLayout || {};
      state = savedDefault(builtInLayout);
      datasets = {
        gc_merged: {cells: parseCSV(mergedScores), genes: parseCSV(mergedGenes), features: parseCSV(mergedFeatures)},
        gc_unmerged: {cells: parseCSV(rawScores), genes: parseCSV(rawGenes), features: parseCSV(rawFeatures)}
      };
      versionManifest = Object.fromEntries(parseCSV(manifestText).map(row => [row.version_key, row]));
      applyDataset(state.datasetKey in datasets ? state.datasetKey : "gc_merged", false);
      bindControls();
      render();
    }).catch(error => {
      $("rrrStatus").textContent = `Data could not be loaded: ${error.message}`;
    });
  } catch (error) {
    $("rrrStatus").textContent = `Data could not be loaded: ${error.message}`;
  }
})();
