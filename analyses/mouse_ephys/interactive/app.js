const levels = ["E1", "E2", "E3", "E4", "E5"];
const defaults = {
  columns: 3, pointSize: 4, opacity: .5, pCutoff: .05,
  showPoints: true, showBoxes: true, showP: true,
  colors: {E1:"#4E79A7",E2:"#F28E2B",E3:"#59A14F",E4:"#2AA6B8",E5:"#B07AA1"}
};
const featureInfo = [
  ["E_Holding.MP..mV.","Holding MP","mV"],["E_Input.resistance..MOhm.","Input resistance","MΩ"],
  ["E_Membrane.time.constant..ms.","Membrane tau","ms"],["E_Rheobase..pA.","Rheobase","pA"],
  ["E_Sag.ratio","Sag ratio","ratio"],["E_Sag.time..s.","Sag time","s"],
  ["E_AP.threshold..mV.","AP threshold","mV"],["E_AP.amplitude..mV.","AP amplitude","mV"],
  ["E_AP.width..ms.","AP width","ms"],["E_Upstroke.to.downstroke.ratio","Up/down ratio","ratio"],
  ["E_Afterhyperpolarization..mV.","AHP","mV"],["E_Max.number.of.APs","Max APs","count"],
  ["E_Latency..ms.","Latency","ms"],["E_Latency....20pA.current..ms.","Latency @ +20 pA","ms"],
  ["E_ISI.adaptation.index","ISI adaptation","index"],["E_ISI.coefficient.of.variation","ISI CV","CV"],
  ["E_AP.amplitude.adaptation.index","AP amplitude adaptation","index"],["E_AP.coefficient.of.variation","AP CV","CV"]
].map(([key,name,unit])=>({key,name,unit}));
let state = structuredClone(defaults);
state.visible = Object.fromEntries(featureInfo.map(f=>[f.key,true]));
state.order = featureInfo.map(f=>f.key);
let raw = [], pairs = [], omnibus = [];

const $ = id => document.getElementById(id);
const svgNS = "http://www.w3.org/2000/svg";
function el(name, attrs={}, text="") { const n=document.createElementNS(svgNS,name); for(const [k,v] of Object.entries(attrs)) n.setAttribute(k,v); if(text) n.textContent=text; return n; }
function quantile(a,q){const s=[...a].sort((x,y)=>x-y),p=(s.length-1)*q,b=Math.floor(p),r=p-b;return s[b+1]===undefined?s[b]:s[b]+r*(s[b+1]-s[b]);}
function hash(s){let h=2166136261;for(let i=0;i<s.length;i++){h^=s.charCodeAt(i);h=Math.imul(h,16777619);}return (h>>>0)/4294967295;}
function fmt(v){if(Math.abs(v)<.001||Math.abs(v)>=10000)return v.toExponential(1).replace("e-","e−");return Number(v.toPrecision(3)).toString();}
function ticks(min,max,n=3){const span=max-min||1, rawStep=span/n, mag=10**Math.floor(Math.log10(rawStep)), norm=rawStep/mag;const step=(norm<1.5?1:norm<3?2:norm<7?5:10)*mag;const out=[];for(let v=Math.ceil(min/step)*step;v<=max+step*1e-6;v+=step)out.push(v);return out;}
function assignLevels(rows){const bins=[];return [...rows].sort((a,b)=>(levels.indexOf(a.Group_2)-levels.indexOf(a.Group_1))-(levels.indexOf(b.Group_2)-levels.indexOf(b.Group_1))).map(r=>{const a=levels.indexOf(r.Group_1),b=levels.indexOf(r.Group_2);let l=0;while((bins[l]||[]).some(([x,y])=>a<=y&&b>=x))l++;(bins[l]??=[]).push([a,b]);return {...r,level:l};});}

function renderPlot(feature){
  const W=360,H=310,m={l:58,r:12,b:38,t:12};
  const sig=state.showP?assignLevels(pairs.filter(p=>p.Feature===feature.key&&p.P_adj_Holm_within_feature<state.pCutoff)):[];
  const levelsN=sig.length?Math.max(...sig.map(p=>p.level))+1:0;
  const annotationH=levelsN?Math.min(118,levelsN*15+8):4;
  const plotTop=m.t+annotationH, plotBottom=H-m.b;
  const values=raw.map(r=>+r[feature.key]), vmin=Math.min(...values), vmax=Math.max(...values), pad=(vmax-vmin||1)*.06;
  const ymin=vmin-pad,ymax=vmax+pad, y=v=>plotBottom-(v-ymin)/(ymax-ymin)*(plotBottom-plotTop);
  const x=i=>m.l+(i+.5)*(W-m.l-m.r)/5;
  const svg=el("svg",{viewBox:`0 0 ${W} ${H}`,role:"img","aria-label":`${feature.name} by E-type`});
  svg.append(el("line",{x1:m.l,y1:plotTop,x2:m.l,y2:plotBottom,stroke:"#111827","stroke-width":1}));
  svg.append(el("line",{x1:m.l,y1:plotBottom,x2:W-m.r,y2:plotBottom,stroke:"#111827","stroke-width":1}));
  for(const t of ticks(vmin,vmax)){svg.append(el("line",{x1:m.l-4,y1:y(t),x2:m.l,y2:y(t),stroke:"#111827"}));const tx=el("text",{x:m.l-7,y:y(t)+3,"text-anchor":"end","font-size":9},fmt(t));svg.append(tx);}
  levels.forEach((g,i)=>svg.append(el("text",{x:x(i),y:H-17,"text-anchor":"middle","font-size":10},g)));
  const ylabel=el("text",{x:14,y:(plotTop+plotBottom)/2,"text-anchor":"middle","font-size":10,transform:`rotate(-90 14 ${(plotTop+plotBottom)/2})`},feature.unit==="ratio"||feature.unit==="index"||feature.unit==="CV"?feature.name:`${feature.name} (${feature.unit})`);svg.append(ylabel);
  levels.forEach((g,i)=>{
    const rows=raw.filter(r=>r.E_type===g), vals=rows.map(r=>+r[feature.key]);
    if(state.showPoints) rows.forEach(r=>svg.append(el("circle",{cx:x(i)+(hash(r.MSN_unique_ID+feature.key)-.5)*27,cy:y(+r[feature.key]),r:state.pointSize/2,fill:state.colors[g],"fill-opacity":state.opacity})));
    if(state.showBoxes){
      const q1=quantile(vals,.25),med=quantile(vals,.5),q3=quantile(vals,.75),iqr=q3-q1,lo=Math.max(Math.min(...vals),q1-1.5*iqr),hi=Math.min(Math.max(...vals),q3+1.5*iqr),notch=Math.min(iqr/2,1.57*iqr/Math.sqrt(vals.length)),w=25,nw=10,cx=x(i);
      svg.append(el("line",{x1:cx,y1:y(lo),x2:cx,y2:y(hi),stroke:"#111827","stroke-width":.8}));
      svg.append(el("line",{x1:cx-8,y1:y(lo),x2:cx+8,y2:y(lo),stroke:"#111827","stroke-width":.8}));
      svg.append(el("line",{x1:cx-8,y1:y(hi),x2:cx+8,y2:y(hi),stroke:"#111827","stroke-width":.8}));
      const d=`M${cx-w},${y(q3)} L${cx-w},${y(med+notch)} L${cx-nw},${y(med)} L${cx-w},${y(med-notch)} L${cx-w},${y(q1)} L${cx+w},${y(q1)} L${cx+w},${y(med-notch)} L${cx+nw},${y(med)} L${cx+w},${y(med+notch)} L${cx+w},${y(q3)} Z`;
      svg.append(el("path",{d,fill:"white","fill-opacity":.86,stroke:"#111827","stroke-width":.8}));
      svg.append(el("line",{x1:cx-nw,y1:y(med),x2:cx+nw,y2:y(med),stroke:"#111827","stroke-width":1}));
    }
  });
  sig.forEach(p=>{const yv=m.t+6+p.level*15,a=levels.indexOf(p.Group_1),b=levels.indexOf(p.Group_2);svg.append(el("line",{x1:x(a),y1:yv+5,x2:x(b),y2:yv+5,stroke:"#111827","stroke-width":.75}));svg.append(el("text",{x:(x(a)+x(b))/2,y:yv+2,"text-anchor":"middle","font-size":8},fmt(+p.P_adj_Holm_within_feature)));});
  return svg;
}

function render(){
  $("plots").style.setProperty("--columns",state.columns);
  const plotRoot=$("plots");plotRoot.innerHTML="";
  const map=Object.fromEntries(featureInfo.map(f=>[f.key,f]));
  state.order.filter(k=>state.visible[k]).forEach(key=>{
    const f=map[key], card=document.createElement("article");card.className="plot-card";card.draggable=true;card.dataset.key=key;
    const head=document.createElement("div");head.className="card-head";const h=document.createElement("h2");h.textContent=f.name;
    const btn=document.createElement("button");btn.textContent="SVG";btn.type="button";btn.onclick=()=>downloadSvg(card.querySelector("svg"),f.name);
    head.append(h,btn);card.append(head,renderPlot(f));plotRoot.append(card);
  });
  $("status").textContent=`Showing ${state.order.filter(k=>state.visible[k]).length} features · drag cards to reorder`;
  bindDrag();
}
function downloadSvg(svg,name){const copy=svg.cloneNode(true);copy.setAttribute("xmlns",svgNS);const blob=new Blob([new XMLSerializer().serializeToString(copy)],{type:"image/svg+xml"});const a=document.createElement("a");a.href=URL.createObjectURL(blob);a.download=`${name.replaceAll(" ","_")}_E1-E5.svg`;a.click();URL.revokeObjectURL(a.href);}
function bindDrag(){let active;document.querySelectorAll(".plot-card").forEach(card=>{card.ondragstart=()=>{active=card.dataset.key;card.classList.add("dragging")};card.ondragend=()=>card.classList.remove("dragging");card.ondragover=e=>e.preventDefault();card.ondrop=e=>{e.preventDefault();const target=card.dataset.key,a=state.order.indexOf(active),b=state.order.indexOf(target);state.order.splice(a,1);state.order.splice(b,0,active);render();};});}
function buildControls(){
  $("colorControls").innerHTML="";levels.forEach(g=>{const label=document.createElement("label"),input=document.createElement("input");input.type="color";input.value=state.colors[g];input.oninput=()=>{state.colors[g]=input.value;render()};label.append(input,g);$("colorControls").append(label);});
  $("featureControls").innerHTML="";featureInfo.forEach(f=>{const label=document.createElement("label"),input=document.createElement("input");input.type="checkbox";input.checked=state.visible[f.key];input.onchange=()=>{state.visible[f.key]=input.checked;render()};label.append(input,f.name);$("featureControls").append(label);});
}
function connect(){
  [["columns","change",e=>state.columns=+e.target.value],["pointSize","input",e=>{state.pointSize=+e.target.value;$("pointSizeValue").value=state.pointSize.toFixed(1)}],["opacity","input",e=>{state.opacity=+e.target.value;$("opacityValue").value=state.opacity.toFixed(2)}],["pCutoff","change",e=>state.pCutoff=+e.target.value],["showPoints","change",e=>state.showPoints=e.target.checked],["showBoxes","change",e=>state.showBoxes=e.target.checked],["showP","change",e=>state.showP=e.target.checked]].forEach(([id,ev,fn])=>$(id).addEventListener(ev,e=>{fn(e);render()}));
  $("allFeatures").onclick=()=>{Object.keys(state.visible).forEach(k=>state.visible[k]=true);buildControls();render()};
  $("noFeatures").onclick=()=>{Object.keys(state.visible).forEach(k=>state.visible[k]=false);buildControls();render()};
  $("reset").onclick=()=>{state={...structuredClone(defaults),visible:Object.fromEntries(featureInfo.map(f=>[f.key,true])),order:featureInfo.map(f=>f.key)};sync();buildControls();render()};
  $("saveConfig").onclick=()=>{const a=document.createElement("a");a.href=URL.createObjectURL(new Blob([JSON.stringify(state,null,2)],{type:"application/json"}));a.download="Etype_plot_settings.json";a.click();URL.revokeObjectURL(a.href)};
  $("loadConfig").onchange=async e=>{try{state={...state,...JSON.parse(await e.target.files[0].text())};sync();buildControls();render()}catch{$("status").textContent="Could not read settings file."}};
}
function sync(){for(const id of ["columns","pointSize","opacity","pCutoff"])$(id).value=state[id];for(const id of ["showPoints","showBoxes","showP"])$(id).checked=state[id];$("pointSizeValue").value=(+state.pointSize).toFixed(1);$("opacityValue").value=(+state.opacity).toFixed(2);}
function registerWebMcp(){
  const context=document.modelContext;if(!context?.registerTool)return;
  const lifecycle=new AbortController();
  const preferences=()=>({columns:state.columns,pointSize:state.pointSize,opacity:state.opacity,pCutoff:state.pCutoff,showPoints:state.showPoints,showBoxes:state.showBoxes,showP:state.showP});
  context.registerTool({name:"read_plot_preferences",title:"Read plot preferences",description:"Read the current visible plot settings.",inputSchema:{type:"object",properties:{},additionalProperties:false},annotations:{readOnlyHint:true,untrustedContentHint:false},execute:()=>preferences()},{signal:lifecycle.signal});
  context.registerTool({name:"configure_plot_view",title:"Configure plot view",description:"Update the visible E-type plot layout and rendering preferences.",inputSchema:{type:"object",properties:{columns:{type:"integer",minimum:1,maximum:3},pointSize:{type:"number",minimum:1,maximum:8},opacity:{type:"number",minimum:.1,maximum:1},pCutoff:{type:"number",enum:[.05,.01,.001]},showPoints:{type:"boolean"},showBoxes:{type:"boolean"},showP:{type:"boolean"}},additionalProperties:false},annotations:{readOnlyHint:false,untrustedContentHint:false},execute:input=>{if(!input||typeof input!=="object")throw new TypeError("Settings must be an object");for(const [key,value] of Object.entries(input)){if(!(key in preferences()))throw new TypeError(`Unknown setting: ${key}`);state[key]=value;}sync();render();return preferences();}},{signal:lifecycle.signal});
}

try {
  [raw,pairs,omnibus]=await Promise.all(["data/raw.json","data/pairwise.json","data/omnibus.json"].map(u=>fetch(u).then(r=>{if(!r.ok)throw new Error(u);return r.json()})));
  buildControls();connect();sync();render();registerWebMcp();
} catch (error) { $("status").textContent=`Data could not be loaded: ${error.message}`; }
