(() => {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const DEFAULT_STORAGE_KEY = "mouse-morph-d1d2-editor-default-v1";
  const M_CLASSES = ["M1", "M2", "M3", "M4"];
  const G_CLASSES = ["G1", "G2", "G3", "G4"];
  const DEFAULT_STATE = {
    identity: "Total", chartTitle: "Total morphology cells", xTitle: "t-SNE 1", yTitle: "t-SNE 2", autoTitle: true,
    showTitle: true, showBackground: true, showLegend: true,
    pointSize: 3, pointOpacity: 1, backgroundSize: 2, backgroundOpacity: 0.45,
    ellipseLevel: 0.85, ellipseWidth: 1, ellipseFillOpacity: 0, showEllipses: true,
    showPointStroke: false, pointStrokeColor: "#ffffff", backgroundColor: "#BDBDBD",
    xMin: -1, xMax: 1, yMin: -1, yMax: 1, rangePadding: 0.04, lockAutoRange: true, equalAspect: true,
    showXAxis: true, showYAxis: true, showAxisTitles: true, showTicks: false, tickCount: 5,
    axisWidth: 1, axisColor: "#111827",
    canvasWidth: 520, canvasHeight: 520, marginLeft: 58, marginRight: 24, marginTop: 42, marginBottom: 58,
    fontFamily: "Arial", fontSize: 11, titleSize: 14, canvasColor: "#ffffff",
    colorM1: "#00468B", colorM2: "#42B540", colorM3: "#ED0000", colorM4: "#0099B4",
    exportWidth: 2.2, exportDpi: 900
  };
  let state = structuredClone(DEFAULT_STATE), cells = [], consensusCells = [], strictCells = [], strictConsensusCells = [], lastPlotted = [], mCountRows = [];
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

  function extent(values) { return [Math.min(...values), Math.max(...values)]; }
  function color(className) { return state[`colorM${className.slice(1)}`]; }
  function isHighlighted(cell) {
    if (state.identity === "GCMerged") return true;
    if (state.identity === "Total") return true;
    if (state.identity === "All") return cell.consensus;
    if (state.identity === "StrictBoth") return cell.strictPass;
    if (state.identity === "StrictMorphConsensus") return cell.strictPass && cell.consensus;
    if (state.identity === "StrictD1") return cell.strictPass && cell.consensus && cell.strictIdentity === "D1";
    if (state.identity === "StrictD2") return cell.strictPass && cell.consensus && cell.strictIdentity === "D2";
    return false;
  }

  function cohortFor(mode) {
    if (mode === "GCMerged") return cells;
    if (mode === "Total") return cells;
    if (mode === "All") return cells.filter(cell => cell.consensus);
    if (mode === "StrictBoth") return cells.filter(cell => cell.strictPass);
    if (mode === "StrictMorphConsensus") return cells.filter(cell => cell.strictPass && cell.consensus);
    if (mode === "StrictD1") return cells.filter(cell => cell.strictPass && cell.consensus && cell.strictIdentity === "D1");
    if (mode === "StrictD2") return cells.filter(cell => cell.strictPass && cell.consensus && cell.strictIdentity === "D2");
    return [];
  }

  function updateMCounts() {
    const conditions = [["Total", "1. Total morphology cells"], ["All", "2. HC–GC morphology consensus"], ["StrictMorphConsensus", "3. Strict D1 + D2 within morphology consensus"], ["StrictD1", "3a. Strict D1"], ["StrictD2", "3b. Strict D2"], ["StrictBoth", "Strict identity incl. morphology-discordant"], ["GCMerged", "GC merged"]];
    mCountRows = conditions.map(([mode, label]) => {
      const cohort = cohortFor(mode), counts = Object.fromEntries(M_CLASSES.map(name => [name, 0]));
      cohort.forEach(cell => { counts[cell.mClass] += 1; });
      return {mode, label, counts, total:cohort.length};
    });
    $("mCountsBody").innerHTML = "";
    mCountRows.forEach(row => {
      const tr = document.createElement("tr");
      if (row.mode === state.identity) tr.className = "current";
      [row.label, ...M_CLASSES.map(name => row.counts[name]), row.total].forEach(value => { const td=document.createElement("td"); td.textContent=value; tr.append(td); });
      $("mCountsBody").append(tr);
    });
    $("mCountsAudit").textContent="Strict identity = Stable_D1/Stable_D2 with agree_RPCA = TRUE and agree_Pearson = TRUE. M1–M4 are the frozen NPC3 / HC K=4 classes.";
    const strictD1 = strictCells.filter(cell => cell.strictIdentity === "D1").length;
    const strictD2 = strictCells.filter(cell => cell.strictIdentity === "D2").length;
    $("identityAudit").textContent=`Audit: 187 morphology cells matched to the identity table; ${strictCells.length} passed strict D1/D2 (${strictD1} D1, ${strictD2} D2), and ${strictConsensusCells.length} additionally passed HC–GC morphology consensus.`;
  }

  function covarianceEllipse(group) {
    if (group.length < 3) return [];
    const mx = group.reduce((sum, cell) => sum + cell.x, 0) / group.length;
    const my = group.reduce((sum, cell) => sum + cell.y, 0) / group.length;
    let a = 0, b = 0, c = 0;
    group.forEach(cell => { const x = cell.x - mx, y = cell.y - my; a += x * x; b += x * y; c += y * y; });
    a /= group.length - 1; b /= group.length - 1; c /= group.length - 1;
    const trace = a + c, disc = Math.sqrt((a - c) ** 2 + 4 * b ** 2);
    const l1 = Math.max(0, (trace + disc) / 2), l2 = Math.max(0, (trace - disc) / 2);
    const theta = 0.5 * Math.atan2(2 * b, a - c), radius = Math.sqrt(-2 * Math.log(1 - state.ellipseLevel));
    return Array.from({length: 161}, (_, index) => {
      const angle = index / 160 * Math.PI * 2, u = radius * Math.sqrt(l1) * Math.cos(angle), v = radius * Math.sqrt(l2) * Math.sin(angle);
      return [mx + u * Math.cos(theta) - v * Math.sin(theta), my + u * Math.sin(theta) + v * Math.cos(theta)];
    });
  }

  function autoDomains() {
    const [xmin, xmax] = extent(consensusCells.map(cell => cell.x)), [ymin, ymax] = extent(consensusCells.map(cell => cell.y));
    const xp = (xmax - xmin || 1) * state.rangePadding, yp = (ymax - ymin || 1) * state.rangePadding;
    return {xMin: xmin - xp, xMax: xmax + xp, yMin: ymin - yp, yMax: ymax + yp};
  }

  function domains(plotWidth, plotHeight) {
    let d = state.lockAutoRange ? autoDomains() : {xMin:+state.xMin,xMax:+state.xMax,yMin:+state.yMin,yMax:+state.yMax};
    if (!(d.xMax > d.xMin)) d.xMax = d.xMin + 1;
    if (!(d.yMax > d.yMin)) d.yMax = d.yMin + 1;
    if (state.equalAspect) {
      const cx = (d.xMin + d.xMax) / 2, cy = (d.yMin + d.yMax) / 2;
      let xSpan = d.xMax - d.xMin, ySpan = d.yMax - d.yMin;
      const target = plotWidth / Math.max(1, plotHeight);
      if (xSpan / ySpan < target) xSpan = ySpan * target; else ySpan = xSpan / target;
      d = {xMin:cx-xSpan/2,xMax:cx+xSpan/2,yMin:cy-ySpan/2,yMax:cy+ySpan/2};
    }
    return d;
  }

  function render() {
    if (!consensusCells.length) return;
    const gcMode = state.identity === "GCMerged";
    const plotCells = cells;
    const classes = gcMode ? G_CLASSES : M_CLASSES;
    const classOf = cell => gcMode ? cell.gcMergedClass : cell.mClass;
    if (state.autoTitle) {
      const titles = {Total:"Total morphology cells", StrictD1:"Strict D1", StrictD2:"Strict D2", StrictBoth:"Strict D1 + D2", StrictMorphConsensus:"Strict D1/D2 + HC–GC consensus", All:"HC–GC consensus", GCMerged:"GC merged"};
      state.chartTitle = titles[state.identity] || state.identity;
    }
    const svg = $("distSvg"), width = state.canvasWidth, height = state.canvasHeight;
    const left = +state.marginLeft, right = width - +state.marginRight, top = +state.marginTop, bottom = height - +state.marginBottom;
    const plotWidth = Math.max(10, right-left), plotHeight = Math.max(10, bottom-top), d = domains(plotWidth, plotHeight);
    const X = value => left + (value-d.xMin)/(d.xMax-d.xMin)*plotWidth;
    const Y = value => bottom - (value-d.yMin)/(d.yMax-d.yMin)*plotHeight;
    const selected = cell => gcMode || isHighlighted(cell);
    const highlighted = plotCells.filter(selected), background = plotCells.filter(cell => !selected(cell));
    lastPlotted = plotCells.map(cell => ({...cell, highlighted:selected(cell)}));

    svg.innerHTML = ""; svg.setAttribute("viewBox", `0 0 ${width} ${height}`); svg.setAttribute("width", `${state.exportWidth}in`); svg.setAttribute("height", `${(state.exportWidth*height/width).toFixed(4)}in`);
    svg.style.fontFamily = state.fontFamily; svg.style.background = state.canvasColor;
    svg.append(svgNode("title", {}, `${state.chartTitle} distribution in morphology t-SNE`));
    const strictMode = state.identity.startsWith("Strict");
    const ellipseReference = strictMode ? strictConsensusCells : consensusCells;
    svg.append(svgNode("desc", {}, gcMode ? `${plotCells.length} final morphology cells coloured by merged GC class. G1-G4 ellipses are fitted from the GC assignments.` : `${highlighted.length} cells highlighted among ${plotCells.length} displayed final morphology cells. M-class ellipses are fitted from ${ellipseReference.length} reference cells.`));
    svg.append(svgNode("rect", {x:0,y:0,width,height,fill:state.canvasColor}));

    if (state.showBackground) background.forEach(cell => svg.append(svgNode("circle", {cx:X(cell.x),cy:Y(cell.y),r:state.backgroundSize,fill:state.backgroundColor,"fill-opacity":state.backgroundOpacity})));
    highlighted.forEach(cell => svg.append(svgNode("circle", {cx:X(cell.x),cy:Y(cell.y),r:state.pointSize,fill:color(classOf(cell)),"fill-opacity":state.pointOpacity,stroke:state.showPointStroke?state.pointStrokeColor:"none","stroke-width":state.showPointStroke?0.5:0})));

    if (state.showEllipses) classes.forEach(className => {
      const reference = gcMode ? cells.filter(cell => cell.gcMergedClass === className) : ellipseReference.filter(cell => cell.mClass === className);
      const points = covarianceEllipse(reference);
      svg.append(svgNode("polygon", {points:points.map(([x,y])=>`${X(x)},${Y(y)}`).join(" "),fill:color(className),"fill-opacity":state.ellipseFillOpacity,stroke:color(className),"stroke-width":state.ellipseWidth}));
    });

    const ticks = Array.from({length:Math.max(2,+state.tickCount)},(_,i)=>i/(Math.max(2,+state.tickCount)-1));
    if (state.showXAxis) svg.append(svgNode("line", {x1:left,y1:bottom,x2:right,y2:bottom,stroke:state.axisColor,"stroke-width":state.axisWidth}));
    if (state.showYAxis) svg.append(svgNode("line", {x1:left,y1:top,x2:left,y2:bottom,stroke:state.axisColor,"stroke-width":state.axisWidth}));
    if (state.showTicks) ticks.forEach(t => {
      const xv=d.xMin+t*(d.xMax-d.xMin), yv=d.yMin+t*(d.yMax-d.yMin), xx=X(xv), yy=Y(yv);
      svg.append(svgNode("line",{x1:xx,y1:bottom,x2:xx,y2:bottom+4,stroke:state.axisColor,"stroke-width":state.axisWidth}));
      svg.append(svgNode("text",{x:xx,y:bottom+state.fontSize+7,"text-anchor":"middle","font-size":state.fontSize,fill:state.axisColor},xv.toFixed(1)));
      svg.append(svgNode("line",{x1:left-4,y1:yy,x2:left,y2:yy,stroke:state.axisColor,"stroke-width":state.axisWidth}));
      svg.append(svgNode("text",{x:left-7,y:yy+state.fontSize*.34,"text-anchor":"end","font-size":state.fontSize,fill:state.axisColor},yv.toFixed(1)));
    });
    if (state.showAxisTitles) {
      svg.append(svgNode("text",{x:(left+right)/2,y:height-10,"text-anchor":"middle","font-size":state.fontSize,fill:state.axisColor},state.xTitle));
      svg.append(svgNode("text",{x:14,y:(top+bottom)/2,"text-anchor":"middle","font-size":state.fontSize,fill:state.axisColor,transform:`rotate(-90 14 ${(top+bottom)/2})`},state.yTitle));
    }
    if (state.showTitle) svg.append(svgNode("text",{x:(left+right)/2,y:Math.max(state.titleSize,top-13),"text-anchor":"middle","font-size":state.titleSize,fill:state.axisColor},state.chartTitle));
    if (state.showLegend) {
      const start=right-150,y=top-17;
      classes.forEach((className,index)=>{const x=start+index*38;svg.append(svgNode("circle",{cx:x,cy:y,r:3,fill:color(className)}));svg.append(svgNode("text",{x:x+6,y:y+3,"font-size":state.fontSize,fill:state.axisColor},className));});
    }

    if (state.lockAutoRange) { Object.assign(state,d); ["xMin","xMax","yMin","yMax"].forEach(id=>{$(id).value=state[id].toFixed(2);}); }
    $("chartTitle").value=state.chartTitle;
    $("distSummary").textContent=gcMode ? `${plotCells.length} final cells · GC-fitted ${Math.round(state.ellipseLevel*100)}% ellipses` : `${highlighted.length} highlighted · ${plotCells.length} displayed · ellipse reference n=${ellipseReference.length}`;
    $("distStatus").textContent=gcMode ? "Points and ellipse geometry use the merged GC assignments (G1–G4) for all 187 final morphology cells." : strictMode ? "Strict D1/D2 comes from the RPCA–Pearson audit; all 187 morphology cells remain visible and the M ellipses use the 168 strict HC–GC-consensus cells." : "All 187 cells are displayed; M-class ellipse geometry uses the 181 HC–GC-consensus cells.";
    const counts=classes.map(className=>`${className}: ${highlighted.filter(cell=>classOf(cell)===className).length}/${(gcMode?plotCells:ellipseReference).filter(cell=>classOf(cell)===className).length}`).join(" · ");
    $("distAudit").textContent=`Displayed / ellipse-reference cells — ${counts}`;
    updateMCounts();
  }

  const numeric=["pointSize","pointOpacity","backgroundSize","backgroundOpacity","ellipseLevel","ellipseWidth","ellipseFillOpacity","xMin","xMax","yMin","yMax","rangePadding","tickCount","axisWidth","canvasWidth","canvasHeight","marginLeft","marginRight","marginTop","marginBottom","fontSize","titleSize","exportWidth","exportDpi"];
  const toggles=["autoTitle","showTitle","showBackground","showLegend","showEllipses","showPointStroke","lockAutoRange","equalAspect","showXAxis","showYAxis","showAxisTitles","showTicks"];
  const texts=["identity","chartTitle","xTitle","yTitle","pointStrokeColor","backgroundColor","axisColor","fontFamily","canvasColor","colorM1","colorM2","colorM3","colorM4"];
  const decimals=new Set(["pointSize","pointOpacity","backgroundSize","backgroundOpacity","ellipseLevel","ellipseWidth","ellipseFillOpacity","rangePadding","axisWidth","exportWidth"]);
  function outputValue(id){const out=$(`${id}Value`);if(!out)return;const v=+state[id];out.value=id==="ellipseLevel"?`${Math.round(v*100)}%`:decimals.has(id)?v.toFixed(2).replace(/0+$/,"").replace(/\.$/,""):String(v);}
  function sync(){numeric.forEach(id=>{if($(id)){ $(id).value=state[id];outputValue(id);}});toggles.forEach(id=>{$(id).checked=state[id];});texts.forEach(id=>{$(id).value=state[id];});}
  function serializedSVG(){const copy=$("distSvg").cloneNode(true);copy.setAttribute("xmlns",NS);return new XMLSerializer().serializeToString(copy);}
  function download(name,blob){const a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
  function savedDefault(){try{const value=localStorage.getItem(DEFAULT_STORAGE_KEY);return value?{...structuredClone(DEFAULT_STATE),...JSON.parse(value)}:null;}catch{return null;}}
  function announce(message){$("distStatus").textContent=message;}
  function exportPNG(){const svg=$("distSvg"),box=svg.viewBox.baseVal,w=Math.round(state.exportWidth*state.exportDpi),h=Math.round(w*box.height/box.width);if(w*h>8e7){$("distStatus").textContent="Requested PNG exceeds the safe 80-megapixel browser limit.";return;}const img=new Image(),url=URL.createObjectURL(new Blob([serializedSVG()],{type:"image/svg+xml"}));img.onload=()=>{const c=document.createElement("canvas");c.width=w;c.height=h;const ctx=c.getContext("2d");ctx.fillStyle=state.canvasColor;ctx.fillRect(0,0,w,h);ctx.drawImage(img,0,0,w,h);URL.revokeObjectURL(url);c.toBlob(blob=>download(`Mouse_${state.identity}_morphology_distribution_${state.exportDpi}dpi.png`,blob),"image/png")};img.src=url;}
  function bind(){numeric.forEach(id=>{if(!$(id))return;$(id).oninput=e=>{state[id]=+e.target.value;outputValue(id);render();};});toggles.forEach(id=>{$(id).onchange=e=>{state[id]=e.target.checked;render();};});texts.forEach(id=>{$(id).oninput=e=>{state[id]=e.target.value;render();};});$("saveAsDefault").onclick=()=>{try{localStorage.setItem(DEFAULT_STORAGE_KEY,JSON.stringify(state));announce("Current settings were saved as this browser's default.");}catch{announce("The browser did not allow the default settings to be saved.");}};$("resetAll").onclick=()=>{state=savedDefault()||structuredClone(DEFAULT_STATE);sync();render();};$("restoreFactory").onclick=()=>{localStorage.removeItem(DEFAULT_STORAGE_KEY);state=structuredClone(DEFAULT_STATE);sync();render();announce("Factory defaults restored; the saved browser default was removed.");};$("exportSvg").onclick=()=>download(`Mouse_${state.identity}_morphology_distribution.svg`,new Blob([serializedSVG()],{type:"image/svg+xml"}));$("exportPng").onclick=exportPNG;$("exportCsv").onclick=()=>{const header="MSN_unique_ID,M_consensus_class,GC_merged_class,strict_D1D2_identity,strict_D1D2_pass,HC_GC_morph_consensus,agree_RPCA,agree_Pearson,tSNE1,tSNE2,highlighted\n",body=lastPlotted.map(c=>[c.id,c.mClass,c.gcMergedClass,c.strictIdentity,c.strictPass,c.consensus,c.agreeRPCA,c.agreePearson,c.x,c.y,c.highlighted].join(",")).join("\n");download(`Mouse_${state.identity}_morphology_distribution_cells.csv`,new Blob([header+body],{type:"text/csv"}));};$("exportMCounts").onclick=()=>{const header="Condition,M1,M2,M3,M4,Total\n",body=mCountRows.map(row=>[row.label,...M_CLASSES.map(name=>row.counts[name]),row.total].join(",")).join("\n");download("Mouse_morphology_conditions_Mclass_counts.csv",new Blob([header+body],{type:"text/csv"}));};$("exportLayout").onclick=()=>download("Mouse_morphology_distribution_layout.json",new Blob([JSON.stringify(state,null,2)],{type:"application/json"}));$("importLayout").onchange=async e=>{const f=e.target.files[0];if(!f)return;try{state={...structuredClone(DEFAULT_STATE),...JSON.parse(await f.text())};sync();render();}catch{$("distStatus").textContent="The selected layout JSON could not be read.";}e.target.value="";};sync();}

  try {
    if(!window.MOUSE_M_RRR_DATA?.assignments)throw new Error("embedded morphology assignments are unavailable");
    if(!window.MOUSE_M_STRICT_IDENTITY_DATA)throw new Error("strict RPCA–Pearson identity audit is unavailable");
    const auditById=new Map(parseCSV(window.MOUSE_M_STRICT_IDENTITY_DATA).map(row=>[row.MSN_unique_ID,row]));
    cells=parseCSV(window.MOUSE_M_RRR_DATA.assignments).map(row=>{const audit=auditById.get(row.MSN_unique_ID)||{};return {id:row.MSN_unique_ID,mClass:row.M_class,gcClass:row.GC_as_M_class,gcMergedClass:row.GC_matched_class,consensus:row.HC_GC_consensus==="True",identity:row.Final_consensus_CellType,strictIdentity:audit.strict_D1D2_identity||"",strictPass:audit.strict_D1D2_pass==="True",auditPresent:audit.audit_present==="True",agreeRPCA:audit.agree_RPCA==="TRUE",agreePearson:audit.agree_Pearson==="TRUE",markerMajor:audit.marker_major||"",pearsonMajor:audit.pearson_major||"",bootstrapTop1:+audit.bootstrap_top1_fraction,bootstrapMargin:+audit.bootstrap_margin,x:+row.tSNE1,y:+row.tSNE2};}).filter(cell=>Number.isFinite(cell.x)&&Number.isFinite(cell.y)&&M_CLASSES.includes(cell.mClass)&&G_CLASSES.includes(cell.gcMergedClass));
    consensusCells=cells.filter(cell=>cell.consensus&&cell.mClass===cell.gcClass);
    strictCells=cells.filter(cell=>cell.strictPass&&["D1","D2"].includes(cell.strictIdentity));
    strictConsensusCells=strictCells.filter(cell=>cell.consensus&&cell.mClass===cell.gcClass);
    state=savedDefault()||structuredClone(DEFAULT_STATE);
    bind();render();
  } catch(error){$("distStatus").textContent=`Data could not be loaded: ${error.message}`;}
})();
