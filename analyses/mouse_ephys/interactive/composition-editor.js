(() => {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const GROUPS = ["E1", "E2", "E3", "E4", "E5"];
  const T_GROUPS = ["D1", "D2", "Ambiguous"];
  const DEFAULT_STORAGE_KEY = "mouse-e-strict-t-composition-default-v1";
  const DEFAULTS = {
    mode: "percent", stackOrder: "D1,D2,Ambiguous", classOrder: "E1,E2,E3,E4,E5",
    canvasWidth: 560, canvasHeight: 460, marginLeft: 72, marginRight: 34, marginTop: 62, marginBottom: 62,
    colorD1: "#ff379b", colorD2: "#00eeb3", colorAmbiguous: "#b8b8b8", barWidth: 0.68, barOpacity: 1,
    borderWidth: 0.6, borderColor: "#ffffff", canvasColor: "#ffffff",
    title: "", titleSize: 14, titleX: 0.5, titleY: 24, yLabel: "Cells (%)",
    fontFamily: "Arial", axisFontSize: 10, yLabelSize: 11, yMin: null, yMax: null,
    tickCount: 5, tickLength: 5, axisWidth: 0.8, axisColor: "#111111",
    showAxes: true, showTicks: true, showTickLabels: true, showXLabels: true, showYLabel: true,
    showAmbiguous: true,
    showSegmentLabels: false, valueFormat: "percent", valueSize: 8, valueColor: "#111111", showTotals: false,
    showLegend: true, legendX: 0.66, legendY: 24, legendSize: 10, legendSwatchWidth: 24, legendGap: 72,
    exportWidth: 2.2, exportDpi: 600
  };
  function savedDefault() {
    try {
      const saved = JSON.parse(localStorage.getItem(DEFAULT_STORAGE_KEY) || "null");
      return saved && typeof saved === "object" ? {...structuredClone(DEFAULTS), ...saved} : structuredClone(DEFAULTS);
    } catch {
      return structuredClone(DEFAULTS);
    }
  }
  let state = savedDefault(), rows = [];
  const $ = id => document.getElementById(id);
  const el = (name, attrs = {}, text = "") => {
    const node = document.createElementNS(NS, name);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    if (text !== "") node.textContent = text;
    return node;
  };
  function parseCSV(text) {
    const clean = value => value.replace(/^\uFEFF/, "").replace(/^"|"$/g, "").replace(/""/g, '"');
    const lines = text.trim().split(/\r?\n/), head = lines.shift().split(",").map(clean);
    return lines.map(line => Object.fromEntries(line.split(",").map((value, i) => [head[i], clean(value)])));
  }
  function summarize(data) {
    return GROUPS.map(group => {
      const subset = data.filter(d => d.E_class === group);
      const counts = Object.fromEntries(T_GROUPS.map(label => [label, subset.filter(d => d.T_identity === label).length]));
      const total = subset.length;
      return {group, ...counts, total, ...Object.fromEntries(T_GROUPS.map(label => [`${label}pct`, counts[label] / total * 100]))};
    });
  }
  function order() {
    const requested = state.classOrder.split(",").map(x => x.trim().toUpperCase()).filter(x => GROUPS.includes(x));
    return [...new Set([...requested, ...GROUPS])].slice(0, GROUPS.length);
  }
  function activeTGroups() { return state.showAmbiguous ? T_GROUPS : ["D1", "D2"]; }
  function displayTotal(row) { return activeTGroups().reduce((sum, group) => sum + row[group], 0); }
  function displayPct(row, group) { const total = displayTotal(row); return total ? row[group] / total * 100 : 0; }
  function domain() {
    const automaticMax = state.mode === "percent" ? 100 : Math.max(...rows.map(displayTotal));
    const min = Number.isFinite(+state.yMin) && state.yMin !== null && state.yMin !== "" ? +state.yMin : 0;
    const max = Number.isFinite(+state.yMax) && state.yMax !== null && state.yMax !== "" ? +state.yMax : automaticMax;
    return max > min ? [min, max] : [0, automaticMax];
  }
  function nice(value) {
    if (state.mode === "percent") return `${Math.round(value)}`;
    return Number.isInteger(value) ? String(value) : value.toFixed(1);
  }
  function segmentText(row, group) {
    const count = row[group], pct = displayPct(row, group);
    if (state.valueFormat === "count") return `${count}`;
    if (state.valueFormat === "both") return `${count} (${pct.toFixed(1)}%)`;
    return `${pct.toFixed(1)}%`;
  }
  function drawLegend(svg) {
    if (!state.showLegend) return;
    const x0 = state.canvasWidth * state.legendX, y = state.legendY, order = activeTGroups();
    order.forEach((group, i) => {
      const x = x0 + i * state.legendGap, color = state[`color${group}`];
      svg.append(el("rect", {x, y: y - state.legendSize, width: state.legendSwatchWidth, height: state.legendSize + 2, fill: color, stroke: "#d0d5dd", "stroke-width": 0.5}));
      svg.append(el("text", {x: x + state.legendSwatchWidth + 5, y, "font-family": state.fontFamily, "font-size": state.legendSize, "font-weight": 600, fill: state.axisColor}, group));
    });
  }
  function render() {
    const svg = $("compositionSvg"), W = state.canvasWidth, H = state.canvasHeight;
    svg.innerHTML = ""; svg.setAttribute("viewBox", `0 0 ${W} ${H}`); svg.setAttribute("width", W); svg.setAttribute("height", H); svg.style.background = state.canvasColor;
    svg.append(el("rect", {x: 0, y: 0, width: W, height: H, fill: state.canvasColor}));
    const left = state.marginLeft, right = W - state.marginRight, top = state.marginTop, bottom = H - state.marginBottom;
    const plotW = Math.max(10, right - left), plotH = Math.max(10, bottom - top), groups = order(), [yMin, yMax] = domain();
    const X = i => left + (i + 0.5) * plotW / groups.length;
    const Y = value => bottom - (value - yMin) / (yMax - yMin) * plotH;
    if (state.title) svg.append(el("text", {x: W * state.titleX, y: state.titleY, "text-anchor": "middle", "font-family": state.fontFamily, "font-size": state.titleSize, "font-weight": 600, fill: state.axisColor}, state.title));
    drawLegend(svg);
    const tickN = Math.max(2, Math.round(state.tickCount));
    if (state.showAxes) {
      svg.append(el("line", {x1: left, y1: top, x2: left, y2: bottom, stroke: state.axisColor, "stroke-width": state.axisWidth}));
      svg.append(el("line", {x1: left, y1: bottom, x2: right, y2: bottom, stroke: state.axisColor, "stroke-width": state.axisWidth}));
    }
    for (let i = 0; i < tickN; i++) {
      const value = yMin + i * (yMax - yMin) / (tickN - 1), y = Y(value);
      if (state.showTicks) svg.append(el("line", {x1: left - state.tickLength, y1: y, x2: left, y2: y, stroke: state.axisColor, "stroke-width": state.axisWidth}));
      if (state.showTickLabels) svg.append(el("text", {x: left - state.tickLength - 4, y: y + state.axisFontSize * 0.34, "text-anchor": "end", "font-family": state.fontFamily, "font-size": state.axisFontSize, fill: state.axisColor}, nice(value)));
    }
    if (state.showYLabel) svg.append(el("text", {x: 16, y: (top + bottom) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bottom) / 2})`, "font-family": state.fontFamily, "font-size": state.yLabelSize, fill: state.axisColor}, state.yLabel));
    const active = activeTGroups(), requestedStack = state.stackOrder.split(",").map(x => x.trim()).filter(x => active.includes(x));
    const stack = [...new Set([...requestedStack, ...active])], slot = plotW / groups.length, width = slot * state.barWidth;
    groups.forEach((group, i) => {
      const row = rows.find(d => d.group === group), x = X(i) - width / 2;
      let cumulative = 0;
      stack.forEach(segment => {
        const value = state.mode === "percent" ? displayPct(row, segment) : row[segment], y0 = Y(cumulative), y1 = Y(cumulative + value), height = Math.abs(y0 - y1);
        svg.append(el("rect", {x, y: Math.min(y0, y1), width, height, fill: state[`color${segment}`], "fill-opacity": state.barOpacity, stroke: state.borderColor, "stroke-width": state.borderWidth}));
        if (state.showSegmentLabels && height >= state.valueSize + 2) svg.append(el("text", {x: X(i), y: (y0 + y1) / 2 + state.valueSize * 0.34, "text-anchor": "middle", "font-family": state.fontFamily, "font-size": state.valueSize, fill: state.valueColor}, segmentText(row, segment)));
        cumulative += value;
      });
      if (state.showXLabels) svg.append(el("text", {x: X(i), y: bottom + state.axisFontSize + 10, "text-anchor": "middle", "font-family": state.fontFamily, "font-size": state.axisFontSize, fill: state.axisColor}, group));
      if (state.showTotals) svg.append(el("text", {x: X(i), y: Y(state.mode === "percent" ? 100 : displayTotal(row)) - 6, "text-anchor": "middle", "font-family": state.fontFamily, "font-size": state.valueSize, fill: state.valueColor}, `n=${displayTotal(row)}`));
    });
    $("compositionSummary").textContent = `${rows.reduce((sum, row) => sum + displayTotal(row), 0)} displayed cells · ${state.showAmbiguous ? "all strict-status categories" : "stable D1/D2 only"} · ${state.mode === "percent" ? "within-class percentage" : "cell counts"} · ${state.exportWidth.toFixed(1)} in wide`;
    $("compositionStatus").textContent = `Stable D1 ${rows.reduce((s,r)=>s+r.D1,0)} · stable D2 ${rows.reduce((s,r)=>s+r.D2,0)} · T-identity ambiguous ${rows.reduce((s,r)=>s+r.Ambiguous,0)} ${state.showAmbiguous ? "(included)" : "(hidden; percentages renormalized)"}`;
    updateTable();
  }
  function updateTable() {
    $("compositionTable").querySelector("tbody").innerHTML = order().map(group => {
      const row = rows.find(d => d.group === group);
      return `<tr><td>${group}</td><td>${row.D1}</td><td>${displayPct(row,"D1").toFixed(2)}</td><td>${row.D2}</td><td>${displayPct(row,"D2").toFixed(2)}</td><td>${row.Ambiguous}</td><td>${state.showAmbiguous ? displayPct(row,"Ambiguous").toFixed(2) : "hidden"}</td><td>${displayTotal(row)}</td></tr>`;
    }).join("");
  }
  const numeric = ["canvasWidth","canvasHeight","marginLeft","marginRight","marginTop","marginBottom","barWidth","barOpacity","borderWidth","titleSize","titleX","titleY","axisFontSize","yLabelSize","tickCount","tickLength","axisWidth","valueSize","legendX","legendY","legendSize","legendSwatchWidth","legendGap","exportWidth","exportDpi"];
  const text = ["mode","stackOrder","classOrder","colorD1","colorD2","colorAmbiguous","borderColor","canvasColor","axisColor","valueColor","title","yLabel","fontFamily","valueFormat"];
  const toggles = ["showAmbiguous","showAxes","showTicks","showTickLabels","showXLabels","showYLabel","showSegmentLabels","showTotals","showLegend"];
  const outputs = ["canvasWidth","canvasHeight","marginLeft","marginRight","marginTop","marginBottom","barWidth","barOpacity","borderWidth","titleSize","titleX","titleY","axisFontSize","yLabelSize","tickLength","axisWidth","valueSize","legendX","legendY","legendSize","legendSwatchWidth","legendGap"];
  function sync() {
    numeric.forEach(id => { if ($(id)) $(id).value = state[id]; });
    text.forEach(id => { if ($(id)) $(id).value = state[id]; });
    toggles.forEach(id => { if ($(id)) $(id).checked = state[id]; });
    $("yMin").value = state.yMin ?? ""; $("yMax").value = state.yMax ?? "";
    outputs.forEach(id => { if ($(`${id}Value`)) $(`${id}Value`).value = state[id]; });
  }
  function bind() {
    numeric.forEach(id => { const node = $(id); node.addEventListener(node.type === "range" ? "input" : "change", () => { state[id] = +node.value; if ($(`${id}Value`)) $(`${id}Value`).value = node.value; render(); }); });
    text.forEach(id => { const node = $(id); node.addEventListener(node.tagName === "SELECT" || node.type === "color" ? "change" : "input", () => { state[id] = node.value; if (id === "mode" && state.yLabel === (node.value === "percent" ? "Cells (n)" : "Cells (%)")) state.yLabel = node.value === "percent" ? "Cells (%)" : "Cells (n)"; sync(); render(); }); });
    toggles.forEach(id => $(id).onchange = () => { state[id] = $(id).checked; render(); });
    ["yMin","yMax"].forEach(id => $(id).oninput = () => { state[id] = $(id).value === "" ? null : +$(id).value; render(); });
    $("saveDefault").onclick = () => {
      localStorage.setItem(DEFAULT_STORAGE_KEY, JSON.stringify(state));
      render();
      $("compositionStatus").textContent = "Current parameters saved as the browser default. They will be restored after reopening or refreshing this page.";
    };
    $("resetAll").onclick = () => { state = savedDefault(); sync(); render(); };
    $("exportSvg").onclick = exportSVG; $("exportPng").onclick = exportPNG; $("exportCsv").onclick = exportCSV;
    $("exportLayout").onclick = () => download("Mouse_E_strict_T_identity_composition_layout.json", new Blob([JSON.stringify(state, null, 2)], {type: "application/json"}));
    $("importLayout").onchange = async event => { const file = event.target.files[0]; if (!file) return; try { state = {...structuredClone(DEFAULTS), ...JSON.parse(await file.text())}; sync(); render(); } catch { $("compositionStatus").textContent = "The selected layout JSON could not be read."; } event.target.value = ""; };
  }
  function serializedSVG() {
    const clone = $("compositionSvg").cloneNode(true), serializer = new XMLSerializer();
    clone.setAttribute("xmlns", NS); clone.setAttribute("width", `${state.exportWidth}in`); clone.setAttribute("height", `${(state.exportWidth * state.canvasHeight / state.canvasWidth).toFixed(4)}in`); clone.style.background = state.canvasColor;
    return `<?xml version="1.0" encoding="UTF-8"?>\n${serializer.serializeToString(clone)}`;
  }
  function download(name, blob) { const a = document.createElement("a"), url = URL.createObjectURL(blob); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 500); }
  function exportSVG() { download("Mouse_E_strict_T_identity_composition.svg", new Blob([serializedSVG()], {type: "image/svg+xml"})); }
  function exportPNG() {
    const width = Math.round(state.exportWidth * state.exportDpi), height = Math.round(width * state.canvasHeight / state.canvasWidth), img = new Image(), url = URL.createObjectURL(new Blob([serializedSVG()], {type: "image/svg+xml"}));
    img.onload = () => { const canvas = document.createElement("canvas"); canvas.width = width; canvas.height = height; const ctx = canvas.getContext("2d"); ctx.fillStyle = state.canvasColor; ctx.fillRect(0, 0, width, height); ctx.drawImage(img, 0, 0, width, height); canvas.toBlob(blob => download(`Mouse_E_strict_T_identity_composition_${state.exportDpi}dpi.png`, blob), "image/png"); URL.revokeObjectURL(url); }; img.src = url;
  }
  function exportCSV() {
    const header = ["E_class","stable_D1_n","stable_D1_percent_current_denominator","stable_D2_n","stable_D2_percent_current_denominator","T_identity_ambiguous_n","T_identity_ambiguous_percent_current_denominator","Displayed_total_n","Ambiguous_included"], body = rows.map(r => [r.group,r.D1,displayPct(r,"D1").toFixed(6),r.D2,displayPct(r,"D2").toFixed(6),r.Ambiguous,state.showAmbiguous?displayPct(r,"Ambiguous").toFixed(6):"NA",displayTotal(r),state.showAmbiguous]);
    download("Mouse_E_strict_T_identity_composition_data.csv", new Blob([[header,...body].map(row => row.join(",")).join("\n")], {type: "text/csv"}));
  }
  fetch("data/e_strict_t_identity.csv").then(response => { if (!response.ok) throw new Error("strict T-identity table unavailable"); return response.text(); }).then(text => {
    rows = summarize(parseCSV(text)); bind(); sync(); render();
  }).catch(error => { $("compositionStatus").textContent = `Data could not be loaded: ${error.message}`; });
})();
