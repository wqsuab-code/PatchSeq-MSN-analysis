(() => {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const M_GROUPS = ["M1", "M2", "M3", "M4"];
  const D_GROUPS = ["D1", "D2"];
  const T_POINT_GROUPS = ["D1", "D2"];
  const DEFAULT_COLORS = {D1: "#ff379b", D2: "#00eeb3", M1: "#00468B", M2: "#42B540", M3: "#ED0000", M4: "#0099B4"};
  const M_FEATURE_ABBREVIATIONS = {
    "Soma circularity": "Soma_Circ",
    "Soma aspect ratio": "Soma_AR",
    "Max radial distance": "Cell_RadDist_Max",
    "Primary neurite number": "Neurite_N",
    "Mean tortuosity": "BD_AvgTort",
    "Total neurite length": "NeuriteLen_Total",
    "Bifurcation points": "Bifurcation_N",
    "Maximum branch order": "BranchOrder_Max",
    "Minimum trunk angle": "TrunkAng_Min",
    "Maximum trunk angle": "TrunkAng_Max"
  };
  const PANEL_ORDER = ["T12", "T13", "M12", "M13"];
  const PANEL_META = {
    T12: {space: "T", pair: [1, 2], side: -1, count: "gene"},
    T13: {space: "T", pair: [1, 3], side: -1, count: "gene"},
    M12: {space: "M", pair: [1, 2], side: 1, count: "feature"},
    M13: {space: "M", pair: [1, 3], side: 1, count: "feature"}
  };
  const DEFAULT_PAIRS = Object.fromEntries(PANEL_ORDER.map(id => [id, [...PANEL_META[id].pair]]));
  const DEFAULT_STATE = {
    geneCount: 9, featureCount: 10, fontSize: 11, labelAlign: "auto", boxOpacity: 0.9,
    cellSize: 2, cellOpacity: 0.72, panelLayout: "row", plotRadius: 90,
    panelWidth: 440, panelHeight: 370, columnGap: 0, rowGap: 10,
    radialPercentile: 99, ellipseLevel: 0.9, loadingScale: 0.76,
    fontFamily: "Arial", arrowWidth: 1.4, connectorWidth: 0.8, ellipseWidth: 2, ellipseOpacity: 0.055,
    circleWidth: 2, axisWidth: 1, titleSize: 18, axisLabelSize: 13,
    axisColor: "#aeb6c2", axisTextColor: "#111111", showTicks: false,
    tickCount: 5, tickDecimals: 1, tickSize: 9,
    titleT: "Transcriptomic space", titleM: "Morphological space", canvasColor: "#ffffff",
    exportWidth: 7.2, exportDpi: 600,
    showBoxes: true, showConnectors: true, showEllipses: true, showLegend: true, showCircle: true,
    showPoints: true,
    axisDisplay: Object.fromEntries(PANEL_ORDER.map(id => [id, {showX: true, showY: true, labelX: true, labelY: true, tickX: true, tickY: true}])),
    pairs: structuredClone(DEFAULT_PAIRS), positions: {}, aliases: {}, colors: {...DEFAULT_COLORS}
  };
  function stateFromLayout(incoming = {}) {
    return {
      ...structuredClone(DEFAULT_STATE),
      ...incoming,
      colors: {...DEFAULT_COLORS, ...(incoming.colors || {})},
      pairs: {...structuredClone(DEFAULT_PAIRS), ...(incoming.pairs || {})},
      axisDisplay: Object.fromEntries(PANEL_ORDER.map(id => [
        id,
        {...DEFAULT_STATE.axisDisplay[id], ...(incoming.axisDisplay?.[id] || {})}
      ])),
      positions: incoming.positions || {},
      aliases: incoming.aliases || {}
    };
  }
  let state = stateFromLayout(window.MOUSE_T_M_RRR_DEFAULT_LAYOUT || {});
  let cells = [], genes = [], features = [], currentLabels = [], drag = null;
  const $ = id => document.getElementById(id);

  function svgNode(name, attrs = {}, text = "") {
    const n = document.createElementNS(NS, name);
    Object.entries(attrs).forEach(([k, v]) => n.setAttribute(k, v));
    if (text !== "") n.textContent = text;
    return n;
  }

  function parseCSV(text) {
    const lines = text.trim().split(/\r?\n/);
    const head = lines.shift().split(",");
    return lines.map(line => {
      const cols = line.split(",");
      return Object.fromEntries(head.map((h, i) => {
        const value = cols[i] ?? "";
        return [h, /^(T_RRR|M_RRR|RRR)\d$/.test(h) ? +value : value];
      }));
    });
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
    const x = col * (state.panelWidth + state.columnGap);
    const y = row * (state.panelHeight + state.rowGap);
    return {...PANEL_META[id], id, pair: state.pairs[id], x, y,
      cx: x + state.panelWidth / 2, cy: y + state.panelHeight / 2 + 4};
  }

  function canvasGeometry() {
    const grid = state.panelLayout === "grid";
    return {
      width: (grid ? 2 : 4) * state.panelWidth + (grid ? 1 : 3) * state.columnGap,
      height: (grid ? 2 : 1) * state.panelHeight + (grid ? state.rowGap : 0)
    };
  }

  function scaledScores(panel) {
    const xs = cells.map(c => +c[`${panel.space}_RRR${panel.pair[0]}`]);
    const ys = cells.map(c => +c[`${panel.space}_RRR${panel.pair[1]}`]);
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
    const xKey = `RRR${panel.pair[0]}`, yKey = `RRR${panel.pair[1]}`;
    const top = source.map(d => ({...d, label: d.Feature, vx: +d[xKey], vy: +d[yKey], mag: Math.hypot(+d[xKey], +d[yKey])}))
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
    const q = points.filter(p => p.T_identity === group);
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
      const title = panel.space === "T" ? state.titleT : state.titleM;
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
        const pointGroup = panel.space === "T" ? p.T_identity : p.M_class;
        g.append(svgNode("circle", {cx: X(p.sx), cy: Y(p.sy), r: state.cellSize, fill: state.colors[pointGroup], "fill-opacity": state.cellOpacity}));
      });

      if (state.showLegend) {
        const legendGroups = panel.space === "T" ? T_POINT_GROUPS : M_GROUPS;
        const step = panel.space === "T" ? 38 : 42, total = step * legendGroups.length;
        const start = panel.cx - total / 2 + step / 2, y = panel.y + state.panelHeight - 36;
        legendGroups.forEach((group, index) => {
          const x = start + index * step;
          g.append(svgNode("circle", {cx: x, cy: y, r: 3, fill: state.colors[group]}));
          g.append(svgNode("text", {x: x + 6, y: y + 3, "font-size": Math.max(8, state.axisLabelSize - 2), fill: "#111111"}, group));
        });
      }

      panelLoadings(panel).forEach(d => {
        const defaultLabel = panel.space === "M" ? (M_FEATURE_ABBREVIATIONS[d.label] || d.label) : d.label;
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
    $("rrrSummary").textContent = `${cells.length} cells · rank 3 · ${state.panelLayout === "row" ? "1 × 4" : "2 × 2"} · ${state.exportWidth.toFixed(1)} in wide`;
    $("rrrStatus").textContent = "Drag labels to reposition them; loading vectors remain fixed to the fitted Mouse T–M RRR model.";
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

  const numeric = ["geneCount","featureCount","fontSize","boxOpacity","cellSize","cellOpacity","plotRadius","panelWidth","panelHeight","columnGap","rowGap","radialPercentile","ellipseLevel","loadingScale","arrowWidth","connectorWidth","ellipseWidth","ellipseOpacity","circleWidth","axisWidth","titleSize","axisLabelSize","tickCount","tickDecimals","tickSize","exportWidth","exportDpi"];
  const toggles = ["showBoxes","showConnectors","showEllipses","showLegend","showCircle","showPoints","showTicks"];
  const texts = ["titleT","titleM","canvasColor","axisColor","axisTextColor"];
  const decimals = new Set(["fontSize","boxOpacity","cellSize","cellOpacity","radialPercentile","ellipseLevel","loadingScale","arrowWidth","connectorWidth","ellipseWidth","ellipseOpacity","circleWidth","axisWidth"]);

  function outputValue(id) {
    const out = $(`${id}Value`); if (!out) return;
    const value = +state[id]; out.value = id === "ellipseLevel" ? `${Math.round(value * 100)}%` : decimals.has(id) ? value.toFixed(id === "ellipseOpacity" ? 3 : 2).replace(/0+$/, "").replace(/\.$/, "") : String(value);
  }

  function syncControls() {
    numeric.forEach(id => { if ($(id)) { $(id).value = state[id]; outputValue(id); } });
    toggles.forEach(id => { $(id).checked = state[id]; }); texts.forEach(id => { $(id).value = state[id]; });
    $("panelLayout").value = state.panelLayout; $("labelAlign").value = state.labelAlign; $("fontFamily").value = state.fontFamily;
    Object.entries(state.pairs).forEach(([id, pair]) => { $(`${id}X`).value = pair[0]; $(`${id}Y`).value = pair[1]; });
    Object.entries(state.axisDisplay).forEach(([id, axis]) => {
      $(`${id}ShowX`).checked = axis.showX; $(`${id}ShowY`).checked = axis.showY;
      $(`${id}LabelX`).checked = axis.labelX; $(`${id}LabelY`).checked = axis.labelY;
      $(`${id}TickX`).checked = axis.tickX; $(`${id}TickY`).checked = axis.tickY;
    });
    Object.entries(state.colors).forEach(([id, value]) => { $(`color-${id}`).value = value; });
  }

  function serializedSVG() {
    const copy = $("rrrSvg").cloneNode(true); copy.setAttribute("xmlns", NS); copy.style.minWidth = "";
    return new XMLSerializer().serializeToString(copy);
  }
  function download(name, blob) { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); }
  function exportSVG() { download("Mouse_T-M_RRR_edited.svg", new Blob([serializedSVG()], {type: "image/svg+xml"})); }
  function exportPNG() {
    const svg = $("rrrSvg"), vb = svg.viewBox.baseVal, width = Math.round(state.exportWidth * state.exportDpi), height = Math.round(width * vb.height / vb.width);
    if (width * height > 8e7) { $("rrrStatus").textContent = "Requested PNG exceeds the safe 80-megapixel browser limit."; return; }
    const img = new Image(), url = URL.createObjectURL(new Blob([serializedSVG()], {type: "image/svg+xml"}));
    img.onload = () => { const canvas = document.createElement("canvas"); canvas.width = width; canvas.height = height; const ctx = canvas.getContext("2d"); ctx.fillStyle = state.canvasColor; ctx.fillRect(0,0,width,height); ctx.drawImage(img,0,0,width,height); URL.revokeObjectURL(url); canvas.toBlob(blob => download(`Mouse_T-M_RRR_${state.exportDpi}dpi.png`, blob), "image/png"); };
    img.src = url;
  }

  function bindControls() {
    const palette = $("colorControls");
    [...T_POINT_GROUPS, ...M_GROUPS].forEach(id => { const label = document.createElement("label"), input = document.createElement("input"); label.append(document.createTextNode(id)); input.type = "color"; input.id = `color-${id}`; input.value = state.colors[id]; input.oninput = event => { state.colors[id] = event.target.value; render(); }; label.append(input); palette.append(label); });
    numeric.forEach(id => { if (!$(id)) return; $(id).oninput = event => { state[id] = +event.target.value; outputValue(id); render(); }; });
    toggles.forEach(id => { $(id).onchange = event => { state[id] = event.target.checked; render(); }; });
    texts.forEach(id => { $(id).oninput = event => { state[id] = event.target.value; render(); }; });
    $("panelLayout").onchange = event => { state.panelLayout = event.target.value; render(); };
    $("labelAlign").onchange = event => { state.labelAlign = event.target.value; render(); };
    $("fontFamily").onchange = event => { state.fontFamily = event.target.value; render(); };
    PANEL_ORDER.forEach(id => [0,1].forEach(axis => { const input = $(`${id}${axis ? "Y" : "X"}`); input.onchange = event => { const pair = state.pairs[id], next = +event.target.value; if (next === pair[1-axis]) { pair[1-axis] = pair[axis]; $(`${id}${axis ? "X" : "Y"}`).value = pair[1-axis]; } pair[axis] = next; render(); }; }));
    PANEL_ORDER.forEach(id => {
      [["ShowX", "showX"], ["ShowY", "showY"], ["LabelX", "labelX"], ["LabelY", "labelY"], ["TickX", "tickX"], ["TickY", "tickY"]].forEach(([suffix, key]) => {
        $(`${id}${suffix}`).onchange = event => { state.axisDisplay[id][key] = event.target.checked; render(); };
      });
    });
    $("resetLabels").onclick = () => { state.positions = {}; state.aliases = {}; render(); };
    $("resetAll").onclick = () => { state = stateFromLayout(window.MOUSE_T_M_RRR_DEFAULT_LAYOUT || {}); syncControls(); render(); };
    $("exportSvg").onclick = exportSVG; $("exportPng").onclick = exportPNG;
    $("exportCsv").onclick = () => { const rows = [["panel","original_label","display_label","x_normalized","y_normalized"], ...currentLabels.map(d => [d.panel,d.original,d.label,d.x,d.y])]; const csv = rows.map(row => row.map(v => `"${String(v).replaceAll('"','""')}"`).join(",")).join("\n"); download("Mouse_T-M_RRR_label_positions.csv", new Blob([csv], {type: "text/csv"})); };
    $("exportLayout").onclick = () => download("Mouse_T-M_RRR_layout.json", new Blob([JSON.stringify(state, null, 2)], {type: "application/json"}));
    $("importLayout").onchange = async event => { const file = event.target.files[0]; if (!file) return; try { const incoming = JSON.parse(await file.text()); state = stateFromLayout(incoming); syncControls(); render(); } catch { $("rrrStatus").textContent = "The selected layout JSON could not be read."; } event.target.value = ""; };
    syncControls();
  }

  try {
    if (!window.MOUSE_M_RRR_DATA) throw new Error("embedded data are unavailable");
    cells = parseCSV(window.MOUSE_M_RRR_DATA.scores);
    const loadings = parseCSV(window.MOUSE_M_RRR_DATA.loadings);
    genes = loadings.filter(d => d.Domain === "Transcriptomic");
    features = loadings.filter(d => d.Domain === "Morphological");
    bindControls(); render();
  } catch (error) {
    $("rrrStatus").textContent = `Data could not be loaded: ${error.message}`;
  }
})();
