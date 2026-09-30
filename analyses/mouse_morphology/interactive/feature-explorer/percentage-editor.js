(() => {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const DEFAULT_STORAGE_KEY = "mouse-morph-percentage-editor-default-v1";
  const CLASSES = ["M1", "M2", "M3", "M4"];
  const GROUPS = ["D1", "D2", "Ambiguous"];
  const DEFAULT_STATE = {
    barOrder: "M1,M2,M3,M4", stackOrder: "D1,D2,Ambiguous", includeAmbiguous: true,
    chartTitle: "D1/D2 and ambiguous composition across M classes", xTitle: "Morphological class", yTitle: "Cells (%)",
    labelM1: "M1", labelM2: "M2", labelM3: "M3", labelM4: "M4",
    showTitle: false, showXAxisTitle: true, showYAxisTitle: true,
    barWidth: 48, barGap: 18, segmentOpacity: 1, segmentStroke: 1, segmentStrokeColor: "#ffffff",
    showPercent: true, showCounts: false, showSampleN: false, percentDecimals: 1, valueSize: 10, valueColor: "#111827",
    showXAxis: true, showYAxis: true, showXLabels: true, showYLabels: true, showGrid: false,
    yMin: 0, yMax: 100, yTickIntervals: 5, axisWidth: 1.25, tickLength: 4, gridWidth: 0.5,
    axisColor: "#111827", gridColor: "#D0D5DD",
    canvasWidth: 420, canvasHeight: 420, marginLeft: 62, marginRight: 24, marginTop: 48, marginBottom: 64,
    fontFamily: "Arial", axisFontSize: 11, titleFontSize: 14, canvasColor: "#ffffff", plotColor: "#ffffff",
    colorD1: "#FF379B", colorD2: "#00EEB3", colorAmbiguous: "#BDBDBD", showLegend: true, legendPosition: "top",
    legendFontSize: 11, legendSwatch: 12, legendGap: 18,
    exportWidth: 2.1, exportDpi: 900
  };
  let state = structuredClone(DEFAULT_STATE);
  let rows = [];
  let composition = [];
  const $ = id => document.getElementById(id);

  function svgNode(name, attrs = {}, text = "") {
    const node = document.createElementNS(NS, name);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    if (text !== "") node.textContent = text;
    return node;
  }

  function parseCSV(text) {
    const lines = text.trim().split(/\r?\n/), header = lines.shift().split(",");
    return lines.map(line => Object.fromEntries(header.map((key, index) => [key, (line.split(",")[index] ?? "").trim()])));
  }

  function calculateComposition() {
    composition = CLASSES.map(mClass => {
      const cells = rows.filter(row => row.M_class === mClass && GROUPS.includes(row.T_identity));
      const counts = Object.fromEntries(GROUPS.map(group => [group, cells.filter(row => row.T_identity === group).length]));
      const activeGroups = state.includeAmbiguous ? GROUPS : ["D1", "D2"];
      const total = activeGroups.reduce((sum, group) => sum + counts[group], 0);
      const percentages = Object.fromEntries(GROUPS.map(group => [group, activeGroups.includes(group) && total ? counts[group] / total * 100 : 0]));
      return {mClass, counts, total, fullTotal: cells.length, percentages};
    });
  }

  function orderedComposition() {
    const requested = state.barOrder.split(",").map(value => value.trim().toUpperCase()).filter(value => CLASSES.includes(value));
    const order = [...new Set(requested)];
    return order.length ? order.map(id => composition.find(row => row.mClass === id)).filter(Boolean) : composition;
  }

  function labelFor(mClass) { return state[`label${mClass}`] || mClass; }
  function colorFor(group) { return state[`color${group}`]; }
  function pct(value) { return `${value.toFixed(state.percentDecimals)}%`; }

  function renderLegend(svg, width, height) {
    if (!state.showLegend) return;
    const order = state.stackOrder.split(",").filter(group => GROUPS.includes(group) && (state.includeAmbiguous || group !== "Ambiguous"));
    const itemWidths = order.map(group => state.legendSwatch + 6 + group.length * state.legendFontSize * 0.65);
    const totalWidth = itemWidths.reduce((sum, width) => sum + width, 0) + state.legendGap * (order.length - 1);
    const horizontal = ["top", "bottom"].includes(state.legendPosition);
    let x = state.legendPosition === "left" ? 8 : state.legendPosition === "right" ? width - state.marginRight + 7 : (width - totalWidth) / 2;
    let y = state.legendPosition === "bottom" ? height - 18 : state.legendPosition === "top" ? 20 : state.marginTop + 10;
    order.forEach((group, index) => {
      const gx = horizontal ? x : x;
      const gy = horizontal ? y : y + index * (state.legendSwatch + state.legendGap);
      svg.append(svgNode("rect", {x: gx, y: gy - state.legendSwatch + 2, width: state.legendSwatch, height: state.legendSwatch, fill: colorFor(group)}));
      svg.append(svgNode("text", {x: gx + state.legendSwatch + 6, y: gy + 2, "font-size": state.legendFontSize, fill: state.axisColor}, group));
      if (horizontal) x += itemWidths[index] + state.legendGap;
    });
  }

  function render() {
    const svg = $("pctSvg"), data = orderedComposition();
    const width = Math.max(120, state.canvasWidth), height = Math.max(120, state.canvasHeight);
    const left = Math.min(width - 20, state.marginLeft), right = Math.max(left + 10, width - state.marginRight);
    const top = Math.min(height - 20, state.marginTop), bottom = Math.max(top + 10, height - state.marginBottom);
    const plotWidth = right - left, plotHeight = bottom - top;
    const yMin = Math.min(+state.yMin, +state.yMax - 1), yMax = Math.max(+state.yMax, yMin + 1);
    const y = value => bottom - (value - yMin) / (yMax - yMin) * plotHeight;
    const totalBarsWidth = data.length * state.barWidth + Math.max(0, data.length - 1) * state.barGap;
    const startX = left + (plotWidth - totalBarsWidth) / 2;
    const stackOrder = state.stackOrder.split(",").filter(group => GROUPS.includes(group) && (state.includeAmbiguous || group !== "Ambiguous"));
    const expectedGroups = state.includeAmbiguous ? 3 : 2;
    if (stackOrder.length !== expectedGroups) stackOrder.splice(0, stackOrder.length, "D1", "D2", ...(state.includeAmbiguous ? ["Ambiguous"] : []));

    svg.innerHTML = "";
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    svg.setAttribute("width", `${state.exportWidth}in`);
    svg.setAttribute("height", `${(state.exportWidth * height / width).toFixed(4)}in`);
    svg.style.fontFamily = state.fontFamily;
    svg.style.background = state.canvasColor;
    svg.append(svgNode("title", {}, "D1, D2 and ambiguous identity composition across Mouse M1 to M4 morphology classes"));
    svg.append(svgNode("desc", {}, data.map(row => `${row.mClass}: D1 ${row.counts.D1} (${pct(row.percentages.D1)}), D2 ${row.counts.D2} (${pct(row.percentages.D2)}), ambiguous ${row.counts.Ambiguous} (${pct(row.percentages.Ambiguous)})`).join("; ")));
    svg.append(svgNode("rect", {x: 0, y: 0, width, height, fill: state.canvasColor}));
    svg.append(svgNode("rect", {x: left, y: top, width: plotWidth, height: plotHeight, fill: state.plotColor}));

    if (state.showTitle) svg.append(svgNode("text", {x: width / 2, y: Math.max(state.titleFontSize, top - 15), "text-anchor": "middle", "font-size": state.titleFontSize, fill: state.axisColor}, state.chartTitle));

    const intervals = Math.max(1, +state.yTickIntervals);
    for (let index = 0; index <= intervals; index++) {
      const value = yMin + (yMax - yMin) * index / intervals, yy = y(value);
      if (state.showGrid && index > 0) svg.append(svgNode("line", {x1: left, y1: yy, x2: right, y2: yy, stroke: state.gridColor, "stroke-width": state.gridWidth}));
      if (state.showYAxis) svg.append(svgNode("line", {x1: left - state.tickLength, y1: yy, x2: left, y2: yy, stroke: state.axisColor, "stroke-width": state.axisWidth}));
      if (state.showYLabels) svg.append(svgNode("text", {x: left - state.tickLength - 4, y: yy + state.axisFontSize * 0.34, "text-anchor": "end", "font-size": state.axisFontSize, fill: state.axisColor}, `${Math.round(value)}%`));
    }

    if (state.showXAxis) svg.append(svgNode("line", {x1: left, y1: bottom, x2: right, y2: bottom, stroke: state.axisColor, "stroke-width": state.axisWidth}));
    if (state.showYAxis) svg.append(svgNode("line", {x1: left, y1: top, x2: left, y2: bottom, stroke: state.axisColor, "stroke-width": state.axisWidth}));

    data.forEach((row, index) => {
      const x = startX + index * (state.barWidth + state.barGap);
      let cumulative = 0;
      stackOrder.forEach(group => {
        const value = row.percentages[group], low = cumulative, high = cumulative + value;
        const yHigh = y(high), yLow = y(low), rectY = Math.min(yHigh, yLow), rectHeight = Math.abs(yLow - yHigh);
        svg.append(svgNode("rect", {x, y: rectY, width: state.barWidth, height: rectHeight, fill: colorFor(group), "fill-opacity": state.segmentOpacity, stroke: state.segmentStrokeColor, "stroke-width": state.segmentStroke, "data-group": group, "data-class": row.mClass}));
        const labels = [];
        if (state.showPercent) labels.push(pct(value));
        if (state.showCounts) labels.push(`n=${row.counts[group]}`);
        if (labels.length && rectHeight >= state.valueSize + 3) svg.append(svgNode("text", {x: x + state.barWidth / 2, y: (yHigh + yLow) / 2 + state.valueSize * 0.34, "text-anchor": "middle", "font-size": state.valueSize, fill: state.valueColor}, labels.join(" · ")));
        cumulative = high;
      });
      if (state.showXAxis) svg.append(svgNode("line", {x1: x + state.barWidth / 2, y1: bottom, x2: x + state.barWidth / 2, y2: bottom + state.tickLength, stroke: state.axisColor, "stroke-width": state.axisWidth}));
      if (state.showXLabels) svg.append(svgNode("text", {x: x + state.barWidth / 2, y: bottom + state.axisFontSize + state.tickLength + 4, "text-anchor": "middle", "font-size": state.axisFontSize, fill: state.axisColor}, labelFor(row.mClass)));
      if (state.showSampleN) svg.append(svgNode("text", {x: x + state.barWidth / 2, y: y(100) - 7, "text-anchor": "middle", "font-size": state.axisFontSize, fill: state.axisColor}, `n=${row.total}`));
    });

    if (state.showXAxisTitle) svg.append(svgNode("text", {x: (left + right) / 2, y: height - 12, "text-anchor": "middle", "font-size": state.axisFontSize, fill: state.axisColor}, state.xTitle));
    if (state.showYAxisTitle) svg.append(svgNode("text", {x: 16, y: (top + bottom) / 2, "text-anchor": "middle", "font-size": state.axisFontSize, fill: state.axisColor, transform: `rotate(-90 16 ${(top + bottom) / 2})`}, state.yTitle));
    renderLegend(svg, width, height);

    const strictN=rows.filter(row=>row.T_identity==="D1"||row.T_identity==="D2").length, ambiguousN=rows.length-strictN;
    $("pctSummary").textContent = `${rows.length} HC–GC consensus · ${strictN} strict D1/D2 · ${ambiguousN} ambiguous/non-strict · ${data.length} bars`;
    $("pctStatus").textContent = state.includeAmbiguous ? "All 181 HC–GC morphology-consensus cells are included; percentages are normalized within M class." : "Ambiguous/non-strict cells are hidden; D1 and D2 are renormalized to 100% within each M class.";
    $("dataBody").innerHTML = data.map(row => `<tr><td>${labelFor(row.mClass)}</td><td>${row.counts.D1}</td><td>${pct(row.percentages.D1)}</td><td>${row.counts.D2}</td><td>${pct(row.percentages.D2)}</td><td>${row.counts.Ambiguous}</td><td>${pct(row.percentages.Ambiguous)}</td><td>${row.total}</td></tr>`).join("");
  }

  const numeric = ["barWidth","barGap","segmentOpacity","segmentStroke","percentDecimals","valueSize","yMin","yMax","yTickIntervals","axisWidth","tickLength","gridWidth","canvasWidth","canvasHeight","marginLeft","marginRight","marginTop","marginBottom","axisFontSize","titleFontSize","legendFontSize","legendSwatch","legendGap","exportWidth","exportDpi"];
  const toggles = ["includeAmbiguous","showTitle","showXAxisTitle","showYAxisTitle","showPercent","showCounts","showSampleN","showXAxis","showYAxis","showXLabels","showYLabels","showGrid","showLegend"];
  const texts = ["barOrder","stackOrder","chartTitle","xTitle","yTitle","labelM1","labelM2","labelM3","labelM4","segmentStrokeColor","valueColor","axisColor","gridColor","fontFamily","canvasColor","plotColor","colorD1","colorD2","colorAmbiguous","legendPosition"];
  const outputDecimals = new Set(["segmentOpacity","segmentStroke","axisWidth","gridWidth","exportWidth"]);

  function outputValue(id) {
    const output = $(`${id}Value`); if (!output) return;
    const value = +state[id]; output.value = outputDecimals.has(id) ? value.toFixed(2).replace(/0+$/, "").replace(/\.$/, "") : String(value);
  }

  function syncControls() {
    numeric.forEach(id => { if ($(id)) { $(id).value = state[id]; outputValue(id); } });
    toggles.forEach(id => { if ($(id)) $(id).checked = state[id]; });
    texts.forEach(id => { if ($(id)) $(id).value = state[id]; });
  }

  function serializedSVG() {
    const copy = $("pctSvg").cloneNode(true); copy.setAttribute("xmlns", NS);
    return new XMLSerializer().serializeToString(copy);
  }

  function download(name, blob) {
    const link = document.createElement("a"); link.href = URL.createObjectURL(blob); link.download = name; link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
  }

  function savedDefault(){
    try { const value=localStorage.getItem(DEFAULT_STORAGE_KEY); return value ? {...structuredClone(DEFAULT_STATE),...JSON.parse(value)} : null; }
    catch { return null; }
  }

  function exportPNG() {
    const svg = $("pctSvg"), box = svg.viewBox.baseVal;
    const width = Math.round(state.exportWidth * state.exportDpi), height = Math.round(width * box.height / box.width);
    if (width * height > 8e7) { $("pctStatus").textContent = "Requested PNG exceeds the safe 80-megapixel browser limit."; return; }
    const image = new Image(), url = URL.createObjectURL(new Blob([serializedSVG()], {type: "image/svg+xml"}));
    image.onload = () => {
      const canvas = document.createElement("canvas"); canvas.width = width; canvas.height = height;
      const context = canvas.getContext("2d"); context.fillStyle = state.canvasColor; context.fillRect(0, 0, width, height); context.drawImage(image, 0, 0, width, height);
      URL.revokeObjectURL(url); canvas.toBlob(blob => download(`Mouse_M_D1-D2_percentage_${state.exportDpi}dpi.png`, blob), "image/png");
    };
    image.src = url;
  }

  function bindControls() {
    numeric.forEach(id => { if (!$(id)) return; $(id).oninput = event => { state[id] = +event.target.value; outputValue(id); render(); }; });
    toggles.forEach(id => { $(id).onchange = event => { state[id] = event.target.checked; if(id==="includeAmbiguous")calculateComposition(); render(); }; });
    texts.forEach(id => { $(id).oninput = event => { state[id] = event.target.value; render(); }; });
    $("saveAsDefault").onclick = () => { try { localStorage.setItem(DEFAULT_STORAGE_KEY,JSON.stringify(state)); $("pctStatus").textContent="Current settings were saved as this browser's default."; } catch { $("pctStatus").textContent="The browser did not allow the default settings to be saved."; } };
    $("resetAll").onclick = () => { state = savedDefault() || structuredClone(DEFAULT_STATE); calculateComposition(); syncControls(); render(); };
    $("restoreFactory").onclick = () => { localStorage.removeItem(DEFAULT_STORAGE_KEY); state=structuredClone(DEFAULT_STATE); calculateComposition(); syncControls(); render(); $("pctStatus").textContent="Factory defaults restored; the saved browser default was removed."; };
    $("exportSvg").onclick = () => download("Mouse_M_D1-D2_percentage.svg", new Blob([serializedSVG()], {type: "image/svg+xml"}));
    $("exportPng").onclick = exportPNG;
    $("exportCsv").onclick = () => {
      const data = orderedComposition(), lines = [["M_class","D1_n","D1_percent","D2_n","D2_percent","Ambiguous_n","Ambiguous_percent","displayed_total_n","full_consensus_total_n"], ...data.map(row => [row.mClass,row.counts.D1,row.percentages.D1,row.counts.D2,row.percentages.D2,row.counts.Ambiguous,row.percentages.Ambiguous,row.total,row.fullTotal])];
      download("Mouse_M_D1-D2_percentage_data.csv", new Blob([lines.map(line => line.join(",")).join("\n")], {type: "text/csv"}));
    };
    $("exportLayout").onclick = () => download("Mouse_M_D1-D2_percentage_layout.json", new Blob([JSON.stringify(state, null, 2)], {type: "application/json"}));
    $("importLayout").onchange = async event => {
      const file = event.target.files[0]; if (!file) return;
      try { state = {...structuredClone(DEFAULT_STATE), ...JSON.parse(await file.text())}; calculateComposition(); syncControls(); render(); }
      catch { $("pctStatus").textContent = "The selected layout JSON could not be read."; }
      event.target.value = "";
    };
    syncControls();
  }

  try {
    if (!window.MOUSE_M_RRR_DATA?.assignments || !window.MOUSE_M_STRICT_IDENTITY_DATA) throw new Error("embedded consensus assignments or strict identity audit are unavailable");
    const auditById=new Map(parseCSV(window.MOUSE_M_STRICT_IDENTITY_DATA).map(row=>[row.MSN_unique_ID,row]));
    rows = parseCSV(window.MOUSE_M_RRR_DATA.assignments)
      .filter(row => row.HC_GC_consensus === "True" && CLASSES.includes(row.M_class) && row.M_class === row.GC_as_M_class)
      .map(row => { const audit=auditById.get(row.MSN_unique_ID)||{}; const strict=audit.strict_D1D2_pass==="True"&&["D1","D2"].includes(audit.strict_D1D2_identity); return {...row,T_identity:strict?audit.strict_D1D2_identity:"Ambiguous",strictPass:strict}; });
    state=savedDefault()||structuredClone(DEFAULT_STATE); calculateComposition(); bindControls(); render();
  } catch (error) {
    $("pctStatus").textContent = `Data could not be loaded: ${error.message}`;
  }
})();
