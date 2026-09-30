(() => {
  "use strict";

  const NS = "http://www.w3.org/2000/svg";
  const STORAGE_KEY = "mouse-e-stability-statistics-default-v2";
  const CLASSES = ["E1", "E2", "E3", "E4", "E5"];
  const FACTORY = {
    mode: "percent",
    canvasWidth: 1160,
    canvasHeight: 650,
    panelGap: 46,
    outerMargin: 54,
    barWidth: 0.58,
    d1Color: "#ef3b9a",
    d2Color: "#00d7a7",
    ambiguousColor: "#9aa0a6",
    excludedColor: "#c9cdd2",
    backgroundColor: "#ffffff",
    textColor: "#171a1f",
    fontFamily: "Arial",
    fontSize: 11,
    titleSize: 14,
    showValues: true,
    showLegend: true,
    includeAmbiguous: true,
    exportWidth: 7.2,
    exportDpi: 900
  };

  const $ = id => document.getElementById(id);
  const el = (name, attrs = {}, textValue = "") => {
    const node = document.createElementNS(NS, name);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    if (textValue !== "") node.textContent = textValue;
    return node;
  };
  const clean = value => value.replace(/^\uFEFF/, "").replace(/^"|"$/g, "").replace(/""/g, '"');
  const norm = value => String(value || "").trim().toLowerCase();
  const isTrue = value => ["true", "1", "yes"].includes(norm(value));
  const isStableT = value => value === "D1" || value === "D2";

  function parseCSV(text) {
    const rows = text.trim().split(/\r?\n/);
    const header = rows.shift().split(",").map(clean);
    return rows.map(line => {
      const fields = line.split(",").map(clean);
      return Object.fromEntries(fields.map((field, index) => [header[index], field]));
    });
  }

  function defaults() {
    try {
      const saved = JSON.parse(localStorage.getItem(STORAGE_KEY) || "null");
      const legacy = JSON.parse(localStorage.getItem("mouse-identity-stability-statistics-default-v1") || "null");
      const candidate = saved || legacy;
      if (!candidate || typeof candidate !== "object") return structuredClone(FACTORY);
      const mapped = {...candidate};
      if ("includeTUnstable" in mapped && !("includeAmbiguous" in mapped)) mapped.includeAmbiguous = mapped.includeTUnstable;
      if ("stableColor" in mapped && !("d1Color" in mapped)) mapped.d1Color = mapped.stableColor;
      if ("oodColor" in mapped && !("d2Color" in mapped)) mapped.d2Color = mapped.oodColor;
      return {...structuredClone(FACTORY), ...mapped};
    } catch {
      return structuredClone(FACTORY);
    }
  }

  let state = defaults();
  let cells = [];

  function getStats() {
    const total = cells.length;
    const consensus = cells.filter(d => isTrue(d.HC_GC_consensus));
    const nonConsensus = total - consensus.length;
    const strict = consensus.filter(d => isStableT(d.T_identity));
    const ambiguous = consensus.filter(d => !isStableT(d.T_identity));
    const byClass = CLASSES.map(cls => {
      const rows = consensus.filter(d => d.E_class === cls);
      const d1 = rows.filter(d => d.T_identity === "D1").length;
      const d2 = rows.filter(d => d.T_identity === "D2").length;
      const amb = rows.length - d1 - d2;
      return {cls, rows, d1, d2, amb, total: rows.length, stable: d1 + d2};
    });
    return {
      total,
      consensus,
      nonConsensus,
      strict,
      ambiguous,
      d1: strict.filter(d => d.T_identity === "D1").length,
      d2: strict.filter(d => d.T_identity === "D2").length,
      byClass
    };
  }

  function pct(value, total) {
    return total ? 100 * value / total : 0;
  }

  function addText(svg, x, y, value, attrs = {}) {
    svg.append(el("text", {x, y, "font-family": state.fontFamily, "font-size": state.fontSize, fill: state.textColor, ...attrs}, value));
  }

  function panelFrame(svg, x, y, w, h, title) {
    svg.append(el("rect", {x, y, width: w, height: h, fill: "none", stroke: "#c8cdd3", "stroke-width": 0.7}));
    addText(svg, x + 4, y - 10, title, {"font-size": state.titleSize, "font-weight": 600});
  }

  function formatValue(value, total) {
    return state.mode === "percent" ? `${pct(value, total).toFixed(1)}%` : String(value);
  }

  function stackTotal(segments, explicitTotal) {
    return state.mode === "percent" ? 100 : (explicitTotal || segments.reduce((sum, d) => sum + d.value, 0));
  }

  function stackedBar(svg, x, y, w, h, segments, explicitTotal) {
    const totalRaw = explicitTotal || segments.reduce((sum, d) => sum + d.value, 0);
    const totalDisplay = stackTotal(segments, explicitTotal);
    let acc = 0;
    segments.forEach(segment => {
      const displayValue = state.mode === "percent" ? pct(segment.value, totalRaw) : segment.value;
      const segmentHeight = h * displayValue / totalDisplay;
      const segmentY = y + h - acc - segmentHeight;
      svg.append(el("rect", {x, y: segmentY, width: w, height: segmentHeight, fill: segment.color}));
      if (state.showValues && segmentHeight > state.fontSize + 6) {
        const textColor = segment.label === "Ambiguous" || segment.label === "Not consensus" ? state.textColor : "#ffffff";
        addText(svg, x + w / 2, segmentY + segmentHeight / 2 + state.fontSize * 0.33, formatValue(segment.value, totalRaw), {
          "text-anchor": "middle",
          "font-size": Math.max(8, state.fontSize - 1),
          fill: textColor
        });
      }
      acc += segmentHeight;
    });
  }

  function legend(svg, x, y, items) {
    if (!state.showLegend) return;
    items.forEach((item, index) => {
      const lx = x + index * 122;
      svg.append(el("rect", {x: lx, y: y - 9, width: 13, height: 9, fill: item.color}));
      addText(svg, lx + 18, y, item.label, {"font-size": Math.max(8, state.fontSize - 1)});
    });
  }

  function drawAxis(svg, x, y, w, h, label) {
    svg.append(el("line", {x1: x, y1: y + h, x2: x + w, y2: y + h, stroke: state.textColor, "stroke-width": 0.8}));
    svg.append(el("line", {x1: x, y1: y, x2: x, y2: y + h, stroke: state.textColor, "stroke-width": 0.8}));
    addText(svg, x - 28, y + h / 2, label, {"text-anchor": "middle", transform: `rotate(-90 ${x - 28} ${y + h / 2})`});
  }

  function render() {
    const s = getStats();
    const svg = $("stabilitySvg");
    const W = state.canvasWidth;
    const H = state.canvasHeight;
    const margin = state.outerMargin;
    const gap = state.panelGap;
    svg.innerHTML = "";
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.setAttribute("width", W);
    svg.setAttribute("height", H);
    svg.style.background = state.backgroundColor;
    svg.append(el("rect", {x: 0, y: 0, width: W, height: H, fill: state.backgroundColor}));

    const panelW = (W - 2 * margin - gap) / 2;
    const panelH = (H - 2 * margin - gap) / 2;
    const panels = [
      {x: margin, y: margin, w: panelW, h: panelH, title: "a  E-only cohort flow"},
      {x: margin + panelW + gap, y: margin, w: panelW, h: panelH, title: "b  HC-GC E consensus"},
      {x: margin, y: margin + panelH + gap, w: panelW, h: panelH, title: "c  Strict T identity by E class"},
      {x: margin + panelW + gap, y: margin + panelH + gap, w: panelW, h: panelH, title: "d  Ambiguous T identity by E class"}
    ];
    panels.forEach(panel => panelFrame(svg, panel.x, panel.y, panel.w, panel.h, panel.title));

    drawCohortFlow(svg, panels[0], s);
    drawConsensus(svg, panels[1], s);
    drawClassComposition(svg, panels[2], s);
    drawAmbiguous(svg, panels[3], s);

    $("summaryLine").textContent = `E cells ${s.total} · HC-GC E consensus ${s.consensus.length} · strict T stable ${s.strict.length} (D1 ${s.d1}, D2 ${s.d2}) · ambiguous ${s.ambiguous.length}`;
    $("statusLine").textContent = `Current E-only flow: ${s.total} active E cells -> ${s.consensus.length} HC-GC E consensus cells (${pct(s.consensus.length, s.total).toFixed(2)}%) -> ${s.strict.length} strict T-stable cells (${pct(s.strict.length, s.consensus.length).toFixed(2)}%). Ambiguous T identity: ${s.ambiguous.length}.`;
    updateTables(s);
  }

  function drawCohortFlow(svg, panel, s) {
    const baseY = panel.y + panel.h - 42;
    const barH = panel.h - 94;
    const barW = Math.min(66, panel.w * 0.14);
    const xs = [panel.x + panel.w * 0.20, panel.x + panel.w * 0.50, panel.x + panel.w * 0.80];
    const steps = [
      {label: "Active E cells", n: s.total, color: state.d1Color},
      {label: "HC-GC consensus", n: s.consensus.length, color: state.d2Color},
      {label: "Strict T stable", n: s.strict.length, color: "#4c78a8"}
    ];
    steps.forEach((step, index) => {
      const h = barH * step.n / s.total;
      svg.append(el("rect", {x: xs[index] - barW / 2, y: baseY - h, width: barW, height: h, fill: step.color}));
      if (state.showValues) addText(svg, xs[index], baseY - h - 7, String(step.n), {"text-anchor": "middle", "font-weight": 600});
      addText(svg, xs[index], baseY + 17, step.label, {"text-anchor": "middle", "font-size": Math.max(8, state.fontSize - 1)});
      if (index < 2) addText(svg, (xs[index] + xs[index + 1]) / 2, panel.y + panel.h * 0.52, "->", {"text-anchor": "middle", "font-size": state.titleSize});
    });
  }

  function drawConsensus(svg, panel, s) {
    const barW = Math.min(92, panel.w * 0.22);
    const x = panel.x + panel.w * 0.5 - barW / 2;
    const baseY = panel.y + panel.h - 42;
    const barH = panel.h - 94;
    stackedBar(svg, x, baseY - barH, barW, barH, [
      {label: "HC-GC consensus", value: s.consensus.length, color: state.d2Color},
      {label: "Not consensus", value: s.nonConsensus, color: state.excludedColor}
    ], s.total);
    addText(svg, panel.x + panel.w / 2, baseY + 18, `${s.total} active E cells`, {"text-anchor": "middle"});
    legend(svg, panel.x + 18, panel.y + 27, [
      {label: "HC-GC E consensus", color: state.d2Color},
      {label: "Not retained", color: state.excludedColor}
    ]);
  }

  function drawClassComposition(svg, panel, s) {
    const chartX = panel.x + 54;
    const chartY = panel.y + 46;
    const chartW = panel.w - 84;
    const chartH = panel.h - 94;
    const slot = chartW / CLASSES.length;
    const barW = slot * state.barWidth;
    drawAxis(svg, chartX, chartY, chartW, chartH, state.mode === "percent" ? "Cells %" : "Cells");
    s.byClass.forEach((row, index) => {
      const x = chartX + index * slot + (slot - barW) / 2;
      const segments = [
        {label: "D1", value: row.d1, color: state.d1Color},
        {label: "D2", value: row.d2, color: state.d2Color}
      ];
      if (state.includeAmbiguous) segments.push({label: "Ambiguous", value: row.amb, color: state.ambiguousColor});
      stackedBar(svg, x, chartY, barW, chartH, segments, state.includeAmbiguous ? row.total : row.stable);
      addText(svg, x + barW / 2, chartY + chartH + 18, row.cls, {"text-anchor": "middle"});
    });
    legend(svg, panel.x + 18, panel.y + 27, [
      {label: "D1", color: state.d1Color},
      {label: "D2", color: state.d2Color},
      {label: "Ambiguous", color: state.ambiguousColor}
    ]);
  }

  function drawAmbiguous(svg, panel, s) {
    const chartX = panel.x + 54;
    const chartY = panel.y + 46;
    const chartW = panel.w - 84;
    const chartH = panel.h - 94;
    const slot = chartW / CLASSES.length;
    const barW = slot * state.barWidth;
    const maxAmb = Math.max(1, ...s.byClass.map(row => state.mode === "percent" ? pct(row.amb, row.total) : row.amb));
    drawAxis(svg, chartX, chartY, chartW, chartH, state.mode === "percent" ? "Ambiguous %" : "Ambiguous n");
    s.byClass.forEach((row, index) => {
      const displayValue = state.mode === "percent" ? pct(row.amb, row.total) : row.amb;
      const h = chartH * displayValue / maxAmb;
      const x = chartX + index * slot + (slot - barW) / 2;
      svg.append(el("rect", {x, y: chartY + chartH - h, width: barW, height: h, fill: state.ambiguousColor}));
      if (state.showValues) addText(svg, x + barW / 2, chartY + chartH - h - 6, state.mode === "percent" ? `${displayValue.toFixed(1)}%` : String(row.amb), {"text-anchor": "middle", "font-size": Math.max(8, state.fontSize - 1)});
      addText(svg, x + barW / 2, chartY + chartH + 18, row.cls, {"text-anchor": "middle"});
    });
  }

  function updateTables(s) {
    const summaryRows = [
      ["Active E cohort", "All cells with E classification inputs", s.total, 100, "Starting population"],
      ["Active E cohort", "Not HC-GC consensus", s.nonConsensus, pct(s.nonConsensus, s.total), "Excluded from current E-consensus analyses"],
      ["HC-GC E consensus", "Retained", s.consensus.length, pct(s.consensus.length, s.total), "HC and GC E-class assignments agree"],
      ["HC-GC E consensus", "Strict T stable", s.strict.length, pct(s.strict.length, s.consensus.length), "Stable D1 or D2 transcriptomic identity"],
      ["HC-GC E consensus", "T ambiguous", s.ambiguous.length, pct(s.ambiguous.length, s.consensus.length), "Not assigned stable D1 or D2 by the strict rule"],
      ["Strict T-stable subset", "D1", s.d1, pct(s.d1, s.strict.length), "Strict D1 identity"],
      ["Strict T-stable subset", "D2", s.d2, pct(s.d2, s.strict.length), "Strict D2 identity"]
    ];
    $("summaryTable").querySelector("tbody").innerHTML = summaryRows.map(row => `<tr><td>${row[0]}</td><td>${row[1]}</td><td>${row[2]}</td><td>${row[3].toFixed(2)}</td><td>${row[4]}</td></tr>`).join("");
    $("jointTable").querySelector("tbody").innerHTML = s.byClass.map(row => {
      return `<tr><td>${row.cls}</td><td>${row.d1}</td><td>${row.d2}</td><td>${row.amb}</td><td>${row.total}</td><td>${pct(row.d1, row.total).toFixed(2)}</td><td>${pct(row.d2, row.total).toFixed(2)}</td><td>${pct(row.amb, row.total).toFixed(2)}</td></tr>`;
    }).join("");
    const sourceRows = [
      ["HC_E", "Hierarchical-clustering E class", "HC-GC E consensus gate"],
      ["GC_E", "Graph-clustering E class", "HC-GC E consensus gate"],
      ["HC_GC_consensus", "Whether HC_E and GC_E are concordant", "Defines the 450-cell E consensus cohort"],
      ["E_class", "Final E1-E5 consensus class", "E-class grouping"],
      ["RPCA_major_class", "Transcriptomic D1/D2 major class", "Strict T identity rule"],
      ["T_identity", "Final strict T identity", "D1, D2, or Ambiguous"],
      ["major_class_agreement", "RPCA/centroid D1-D2 agreement flag", "Strict T identity rule"]
    ];
    $("kTable").querySelector("tbody").innerHTML = sourceRows.map(row => `<tr><td>${row[0]}</td><td>${row[1]}</td><td>${row[2]}</td></tr>`).join("");
  }

  const numeric = ["canvasWidth", "canvasHeight", "panelGap", "outerMargin", "barWidth", "fontSize", "titleSize", "exportWidth", "exportDpi"];
  const textFields = ["mode", "d1Color", "d2Color", "ambiguousColor", "excludedColor", "backgroundColor", "textColor", "fontFamily"];
  const toggles = ["showValues", "showLegend", "includeAmbiguous"];

  function sync() {
    numeric.forEach(id => {
      $(id).value = state[id];
      const output = $(`${id}Value`);
      if (output) output.value = state[id];
    });
    textFields.forEach(id => { $(id).value = state[id]; });
    toggles.forEach(id => { $(id).checked = state[id]; });
  }

  function bind() {
    numeric.forEach(id => {
      $(id).addEventListener($(id).type === "range" ? "input" : "change", () => {
        state[id] = +$(id).value;
        const output = $(`${id}Value`);
        if (output) output.value = $(id).value;
        render();
      });
    });
    textFields.forEach(id => {
      $(id).addEventListener($(id).tagName === "SELECT" || $(id).type === "color" ? "change" : "input", () => {
        state[id] = $(id).value;
        render();
      });
    });
    toggles.forEach(id => {
      $(id).onchange = () => {
        state[id] = $(id).checked;
        render();
      };
    });
    $("saveDefault").onclick = () => {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
      render();
      $("statusLine").textContent = "Current E-only figure parameters saved as the browser default.";
    };
    $("resetAll").onclick = () => {
      state = defaults();
      sync();
      render();
    };
    $("exportSvg").onclick = exportSVG;
    $("exportPng").onclick = exportPNG;
    $("exportCsv").onclick = () => download("Mouse_E_stability_cell_data.csv", new Blob([toCSV(cells)], {type: "text/csv"}));
    $("exportSummary").onclick = exportSummary;
    $("exportLayout").onclick = () => download("Mouse_E_stability_layout.json", new Blob([JSON.stringify(state, null, 2)], {type: "application/json"}));
    $("importLayout").onchange = async event => {
      const file = event.target.files[0];
      if (!file) return;
      try {
        const imported = JSON.parse(await file.text());
        if ("includeTUnstable" in imported && !("includeAmbiguous" in imported)) imported.includeAmbiguous = imported.includeTUnstable;
        state = {...structuredClone(FACTORY), ...imported};
        sync();
        render();
      } catch {
        $("statusLine").textContent = "The selected layout JSON could not be read.";
      }
      event.target.value = "";
    };
  }

  function serializedSVG() {
    const clone = $("stabilitySvg").cloneNode(true);
    const serializer = new XMLSerializer();
    clone.setAttribute("xmlns", NS);
    clone.setAttribute("width", `${state.exportWidth}in`);
    clone.setAttribute("height", `${(state.exportWidth * state.canvasHeight / state.canvasWidth).toFixed(4)}in`);
    return `<?xml version="1.0" encoding="UTF-8"?>\n${serializer.serializeToString(clone)}`;
  }

  function download(name, blob) {
    const link = document.createElement("a");
    const url = URL.createObjectURL(blob);
    link.href = url;
    link.download = name;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 500);
  }

  function exportSVG() {
    download("Mouse_E_stability_statistics.svg", new Blob([serializedSVG()], {type: "image/svg+xml"}));
  }

  function exportPNG() {
    const width = Math.round(state.exportWidth * state.exportDpi);
    const height = Math.round(width * state.canvasHeight / state.canvasWidth);
    const image = new Image();
    const url = URL.createObjectURL(new Blob([serializedSVG()], {type: "image/svg+xml"}));
    image.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = width;
      canvas.height = height;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = state.backgroundColor;
      ctx.fillRect(0, 0, width, height);
      ctx.drawImage(image, 0, 0, width, height);
      canvas.toBlob(blob => download(`Mouse_E_stability_statistics_${state.exportDpi}dpi.png`, blob), "image/png");
      URL.revokeObjectURL(url);
    };
    image.src = url;
  }

  function toCSV(data) {
    if (!data.length) return "";
    const header = Object.keys(data[0]);
    return [header, ...data.map(row => header.map(key => row[key]))]
      .map(row => row.map(value => `"${String(value ?? "").replace(/"/g, '""')}"`).join(","))
      .join("\n");
  }

  function exportSummary() {
    const s = getStats();
    const lines = [
      "population,status,n,total,percent",
      `Active_E,all,${s.total},${s.total},100.000000`,
      `Active_E,HC_GC_E_consensus,${s.consensus.length},${s.total},${pct(s.consensus.length, s.total).toFixed(6)}`,
      `Active_E,not_HC_GC_consensus,${s.nonConsensus},${s.total},${pct(s.nonConsensus, s.total).toFixed(6)}`,
      `HC_GC_E_consensus,strict_T_stable,${s.strict.length},${s.consensus.length},${pct(s.strict.length, s.consensus.length).toFixed(6)}`,
      `HC_GC_E_consensus,T_ambiguous,${s.ambiguous.length},${s.consensus.length},${pct(s.ambiguous.length, s.consensus.length).toFixed(6)}`,
      `Strict_T_stable,D1,${s.d1},${s.strict.length},${pct(s.d1, s.strict.length).toFixed(6)}`,
      `Strict_T_stable,D2,${s.d2},${s.strict.length},${pct(s.d2, s.strict.length).toFixed(6)}`
    ];
    s.byClass.forEach(row => {
      lines.push(`${row.cls},D1,${row.d1},${row.total},${pct(row.d1, row.total).toFixed(6)}`);
      lines.push(`${row.cls},D2,${row.d2},${row.total},${pct(row.d2, row.total).toFixed(6)}`);
      lines.push(`${row.cls},Ambiguous,${row.amb},${row.total},${pct(row.amb, row.total).toFixed(6)}`);
    });
    download("Mouse_E_stability_summary.csv", new Blob([lines.join("\n")], {type: "text/csv"}));
  }

  fetch("data/e_stability_all493.csv")
    .then(response => {
      if (!response.ok) throw new Error("E stability table unavailable");
      return response.text();
    })
    .then(text => {
      cells = parseCSV(text);
      bind();
      sync();
      render();
    })
    .catch(error => {
      $("statusLine").textContent = `Data could not be loaded: ${error.message}`;
    });
})();
