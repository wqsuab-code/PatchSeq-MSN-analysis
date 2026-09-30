const NS = 'http://www.w3.org/2000/svg';
const STORAGE_KEY = 'macaque-m4-radar-frozen-zrange-default-v1';
const classes = ['M1', 'M2', 'M3', 'M4'];
const features = [
  {id:'basal_dendrite_max_path_distance', label:'BD-MaxPD'},
  {id:'basal_dendrite_bias_medial', label:'BD-Bias-Med'},
  {id:'basal_dendrite_max_branch_order', label:'BD-MaxBO'},
  {id:'basal_dendrite_num_branches', label:'BD-NBranch'},
  {id:'basal_dendrite_total_length', label:'BD-TotLen'},
  {id:'basal_dendrite_stem_exit_MedialLateral', label:'BD-StemExit-ML'},
  {id:'basal_dendrite_stem_exit_ventral', label:'BD-StemExit-Ven'},
  {id:'basal_dendrite_stem_exit_dorsal', label:'BD-StemExit-Dor'},
  {id:'basal_dendrite_max_euclidean_distance', label:'BD-MaxED'},
  {id:'basal_dendrite_extent_dorsal', label:'BDE-Dor'}
];
const expected = {M1:43,M2:42,M3:20,M4:12};
const defaults = {
  canvasWidth:1160, canvasHeight:300, sideMargin:34, panelGap:14, radarRadius:92, centerY:156,
  showTitles:true, titleSize:15,
  classes:Object.fromEntries(classes.map(c=>[c,{enabled:true,title:`${c} (n=${expected[c]})`}])) ,
  colors:{M1:'#1F77B4',M2:'#D9A400',M3:'#8C564B',M4:'#E377C2'},
  scaleMode:'zscore', zMin:-3, zMax:3, centerStatistic:'median', showCells:true, cellWidth:0.35, cellOpacity:0.04,
  showIqr:true, iqrOpacity:0.12, centerWidth:2.0, centerFillOpacity:0.0, curveTension:0.82,
  ringCount:6, gridWidth:0.7, gridOpacity:0.42, labelSize:11, labelDistance:18,
  labelOrientation:'tangent', showRadialValues:false, radialValueSize:9,
  featureOrder:features.map(f=>f.id), featureLabels:Object.fromEntries(features.map(f=>[f.id,f.label])),
  centerColor:'#202020', gridColor:'#A8A8A8', textColor:'#303030', canvasColor:'#FFFFFF',
  exportWidthIn:5.2, exportDpi:600
};

const $ = id => document.getElementById(id);
const svg = $('radarSvg');
const el = (tag, attrs={}, text='') => {
  const node = document.createElementNS(NS, tag);
  for (const [key,value] of Object.entries(attrs)) node.setAttribute(key,value);
  if (text) node.textContent = text;
  return node;
};
function mergeState(value={}) {
  return {
    ...structuredClone(defaults), ...value,
    colors:{...defaults.colors,...(value.colors||{})},
    classes:{...structuredClone(defaults.classes),...(value.classes||{})},
    featureLabels:{...defaults.featureLabels,...(value.featureLabels||{})},
    featureOrder:Array.isArray(value.featureOrder) && value.featureOrder.length===features.length ? value.featureOrder : [...defaults.featureOrder]
  };
}
function savedDefault(){try{return mergeState(JSON.parse(localStorage.getItem(STORAGE_KEY)||'{}'))}catch{return structuredClone(defaults)}}
let state=savedDefault(), rows=[], sourceCsv='';

function parseCsv(text){
  const [header,...lines]=text.trim().split(/\r?\n/), keys=header.split(',');
  return lines.map(line=>Object.fromEntries(line.split(',').map((v,i)=>[keys[i],v])));
}
function quantile(values,p){const a=[...values].sort((x,y)=>x-y),i=(a.length-1)*p,l=Math.floor(i),h=Math.ceil(i);return a[l]+(a[h]-a[l])*(i-l)}
function statistic(values,kind){return kind==='mean'?values.reduce((a,b)=>a+b,0)/values.length:quantile(values,.5)}
function polar(cx,cy,r,index,n){const angle=-Math.PI/2+index*2*Math.PI/n;return [cx+Math.cos(angle)*r,cy+Math.sin(angle)*r,angle]}
function closedPath(points,tension=0){
  if (!points.length) return '';
  if (tension<=0) return points.map((p,i)=>(i?'L':'M')+p[0].toFixed(2)+' '+p[1].toFixed(2)).join(' ')+' Z';
  const n=points.length, k=tension/6; let d=`M${points[0][0].toFixed(2)} ${points[0][1].toFixed(2)}`;
  for(let i=0;i<n;i++){
    const p0=points[(i-1+n)%n],p1=points[i],p2=points[(i+1)%n],p3=points[(i+2)%n];
    const c1=[p1[0]+(p2[0]-p0[0])*k,p1[1]+(p2[1]-p0[1])*k];
    const c2=[p2[0]-(p3[0]-p1[0])*k,p2[1]-(p3[1]-p1[1])*k];
    d+=` C${c1[0].toFixed(2)} ${c1[1].toFixed(2)} ${c2[0].toFixed(2)} ${c2[1].toFixed(2)} ${p2[0].toFixed(2)} ${p2[1].toFixed(2)}`;
  }
  return d+' Z';
}
function featurePoints(values,cx,cy,radius){return values.map((v,i)=>polar(cx,cy,Math.max(0,Math.min(1,v))*radius,i,values.length))}
function featureValue(row,id){if(state.scaleMode==='percentile')return +row[id+'__scaled'];const span=state.zMax-state.zMin;return Math.max(0,Math.min(1,(+row[id+'__z']-state.zMin)/span))}
function featureValues(row){return state.featureOrder.map(id=>featureValue(row,id))}

function classControls(){
  $('classControls').innerHTML=classes.map(c=>`<div class="panel-row"><span><input data-class-enabled="${c}" type="checkbox"><strong>${c}</strong></span><input data-class-title="${c}" type="text"><button data-class-only="${c}" type="button">Only</button></div>`).join('');
  document.querySelectorAll('[data-class-enabled]').forEach(n=>n.onchange=()=>{state.classes[n.dataset.classEnabled].enabled=n.checked;render()});
  document.querySelectorAll('[data-class-title]').forEach(n=>n.oninput=()=>{state.classes[n.dataset.classTitle].title=n.value;render()});
  document.querySelectorAll('[data-class-only]').forEach(n=>n.onclick=()=>{classes.forEach(c=>state.classes[c].enabled=c===n.dataset.classOnly);sync();render()});
}
function featureControls(){
  $('featureControls').innerHTML=state.featureOrder.map((id,i)=>`<div class="feature-row"><span>${i+1}</span><input data-feature-label="${id}" value="${state.featureLabels[id]}"><button data-feature-up="${id}" type="button" aria-label="Move ${id} up">↑</button><button data-feature-down="${id}" type="button" aria-label="Move ${id} down">↓</button></div>`).join('');
  document.querySelectorAll('[data-feature-label]').forEach(n=>n.oninput=()=>{state.featureLabels[n.dataset.featureLabel]=n.value;render()});
  const move=(id,delta)=>{const i=state.featureOrder.indexOf(id),j=Math.max(0,Math.min(state.featureOrder.length-1,i+delta));if(i!==j){[state.featureOrder[i],state.featureOrder[j]]=[state.featureOrder[j],state.featureOrder[i]];featureControls();render()}};
  document.querySelectorAll('[data-feature-up]').forEach(n=>n.onclick=()=>move(n.dataset.featureUp,-1));
  document.querySelectorAll('[data-feature-down]').forEach(n=>n.onclick=()=>move(n.dataset.featureDown,1));
}
function colorControls(){
  $('colorControls').innerHTML=classes.map(c=>`<label>${c}<input data-class-color="${c}" type="color"></label>`).join('');
  document.querySelectorAll('[data-class-color]').forEach(n=>n.oninput=()=>{state.colors[n.dataset.classColor]=n.value;render()});
}
function bindControls(){
  classControls(); featureControls(); colorControls();
  const booleans=new Set(['showTitles','showCells','showIqr','showRadialValues']);
  const strings=new Set(['scaleMode','centerStatistic','labelOrientation','centerColor','gridColor','textColor','canvasColor']);
  const ids=['canvasWidth','canvasHeight','sideMargin','panelGap','radarRadius','centerY','showTitles','titleSize','scaleMode','zMin','zMax','centerStatistic','showCells','cellWidth','cellOpacity','showIqr','iqrOpacity','centerWidth','centerFillOpacity','curveTension','ringCount','gridWidth','gridOpacity','labelSize','labelDistance','labelOrientation','showRadialValues','radialValueSize','centerColor','gridColor','textColor','canvasColor','exportWidthIn','exportDpi'];
  ids.forEach(id=>{const n=$(id);if(!n)return;n.addEventListener(n.type==='range'?'input':'change',()=>{state[id]=booleans.has(id)?n.checked:strings.has(id)?n.value:+n.value;render()})});
  $('saveDefault').onclick=()=>{localStorage.setItem(STORAGE_KEY,JSON.stringify(state));$('radarStatus').textContent='Current radar settings saved as the browser default.'};
  $('restoreDefault').onclick=()=>{state=savedDefault();classControls();featureControls();colorControls();sync();render()};
  $('resetAll').onclick=()=>{state=structuredClone(defaults);classControls();featureControls();colorControls();sync();render()};
  $('exportSvg').onclick=exportSvg; $('exportPng').onclick=exportPng; $('exportPdf').onclick=exportPdf; $('exportCsv').onclick=exportCsv; $('exportLayout').onclick=exportLayout; $('importLayout').onchange=importLayout;
}
function sync(){
  for(const [key,value] of Object.entries(state)){
    const n=$(key); if(!n||typeof value==='object')continue;
    if(n.type==='checkbox')n.checked=value; else n.value=value;
    const out=$(key+'Value'); if(out)out.value=Number.isInteger(value)?String(value):(+value).toFixed(2).replace(/0+$/,'').replace(/\.$/,'');
  }
  document.querySelectorAll('[data-class-enabled]').forEach(n=>n.checked=state.classes[n.dataset.classEnabled].enabled);
  document.querySelectorAll('[data-class-title]').forEach(n=>n.value=state.classes[n.dataset.classTitle].title);
  document.querySelectorAll('[data-class-color]').forEach(n=>n.value=state.colors[n.dataset.classColor]);
}
function render(){
  if(!rows.length)return;
  sync();
  if(state.scaleMode==='zscore' && !(state.zMin < state.zMax)){$('radarStatus').textContent='Z-score lower bound must be smaller than the upper bound.';return;}
  const visible=classes.filter(c=>state.classes[c].enabled); if(!visible.length){$('radarStatus').textContent='Select at least one M class.';return}
  const width=state.canvasWidth,height=state.canvasHeight,n=visible.length;
  svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.replaceChildren();svg.style.background=state.canvasColor;
  svg.append(el('rect',{x:0,y:0,width,height,fill:state.canvasColor}));
  const usable=width-2*state.sideMargin-state.panelGap*(n-1),panelW=usable/n;
  visible.forEach((group,panelIndex)=>{
    const groupRows=rows.filter(r=>r.M_class===group),cx=state.sideMargin+panelIndex*(panelW+state.panelGap)+panelW/2,cy=state.centerY,r=Math.min(state.radarRadius,panelW/2-state.labelDistance-8,height-state.centerY-state.labelDistance-8);
    const g=el('g',{'data-class':group});svg.append(g);
    for(let ring=1;ring<=state.ringCount;ring++)g.append(el('circle',{cx,cy,r:r*ring/state.ringCount,fill:'none',stroke:state.gridColor,'stroke-width':state.gridWidth,'stroke-opacity':state.gridOpacity}));
    state.featureOrder.forEach((id,i)=>{const [x,y,angle]=polar(cx,cy,r,i,state.featureOrder.length);g.append(el('line',{x1:cx,y1:cy,x2:x,y2:y,stroke:state.gridColor,'stroke-width':state.gridWidth,'stroke-opacity':state.gridOpacity}));const [lx,ly]=polar(cx,cy,r+state.labelDistance,i,state.featureOrder.length);let rotation=0;if(state.labelOrientation==='tangent'){rotation=angle*180/Math.PI+90;if(rotation>90&&rotation<270)rotation-=180;if(rotation<-90)rotation+=180}g.append(el('text',{x:lx,y:ly,'text-anchor':'middle','dominant-baseline':'middle','font-size':state.labelSize,fill:state.textColor,transform:`rotate(${rotation} ${lx} ${ly})`},state.featureLabels[id]));});
    if(state.showRadialValues)for(let ring=1;ring<=state.ringCount;ring++){const fraction=ring/state.ringCount,value=state.scaleMode==='zscore'?state.zMin+(state.zMax-state.zMin)*fraction:fraction;g.append(el('text',{x:cx+4,y:cy-r*fraction+state.radialValueSize/3,'font-size':state.radialValueSize,fill:state.textColor,'fill-opacity':.7},value.toFixed(1)))}
    const color=state.colors[group];
    if(state.showCells)for(const row of groupRows){const path=el('path',{d:closedPath(featurePoints(featureValues(row),cx,cy,r),state.curveTension),fill:'none',stroke:color,'stroke-width':state.cellWidth,'stroke-opacity':state.cellOpacity,'vector-effect':'non-scaling-stroke'});path.append(el('title',{},row.cell_label));g.append(path)}
    const columns=state.featureOrder.map(id=>groupRows.map(row=>featureValue(row,id)));
    const q25=columns.map(v=>quantile(v,.25)),q75=columns.map(v=>quantile(v,.75)),center=columns.map(v=>statistic(v,state.centerStatistic));
    if(state.showIqr){const outer=featurePoints(q75,cx,cy,r),inner=featurePoints(q25,cx,cy,r).reverse();g.append(el('path',{d:closedPath(outer,state.curveTension)+' '+closedPath(inner,state.curveTension),fill:color,'fill-opacity':state.iqrOpacity,'fill-rule':'evenodd',stroke:'none'}))}
    const centerPath=closedPath(featurePoints(center,cx,cy,r),state.curveTension);g.append(el('path',{d:centerPath,fill:color,'fill-opacity':state.centerFillOpacity,stroke:'none'}));g.append(el('path',{d:centerPath,fill:'none',stroke:state.centerColor,'stroke-width':state.centerWidth,'vector-effect':'non-scaling-stroke'}));
    if(state.showTitles)g.append(el('text',{x:cx,y:Math.max(state.titleSize+4,cy-r-state.labelDistance-state.labelSize-16),'text-anchor':'middle','font-size':state.titleSize,'font-weight':600,style:`fill:${color}`},state.classes[group].title));
  });
  const pxW=Math.round(state.exportWidthIn*state.exportDpi),pxH=Math.round(pxW*height/width);$('exportPixels').value=`${pxW.toLocaleString()} × ${pxH.toLocaleString()} px`;
  $('radarStatus').textContent=`Mouse-E-style radar · frozen ten morphology indicators · ${state.scaleMode==='zscore'?`Z-score ${state.zMin} to ${state.zMax}`:'percentile 0 to 1'} · ${visible.join(', ')} · ${rows.length} frozen consensus cells · ${state.centerStatistic} with ${state.showIqr?'IQR':'no interval'}`;
}

function serialized(){const copy=svg.cloneNode(true);copy.setAttribute('xmlns',NS);copy.setAttribute('width',state.canvasWidth);copy.setAttribute('height',state.canvasHeight);return new XMLSerializer().serializeToString(copy)}
function download(blob,name){const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}
function exportSvg(){download(new Blob([serialized()],{type:'image/svg+xml'}),'Macaque_M1-M4_top10_radar.svg')}
async function rasterCanvas(){const w=Math.round(state.exportWidthIn*state.exportDpi),h=Math.round(w*state.canvasHeight/state.canvasWidth);if(w*h>8e7)throw Error('Requested image exceeds 80 megapixels');const image=new Image(),url=URL.createObjectURL(new Blob([serialized()],{type:'image/svg+xml'}));await new Promise((ok,no)=>{image.onload=ok;image.onerror=no;image.src=url});const c=document.createElement('canvas');c.width=w;c.height=h;const x=c.getContext('2d',{alpha:false});x.fillStyle=state.canvasColor;x.fillRect(0,0,w,h);x.drawImage(image,0,0,w,h);URL.revokeObjectURL(url);return c}
async function exportPng(){try{const c=await rasterCanvas(),blob=await new Promise(ok=>c.toBlob(ok,'image/png'));download(blob,`Macaque_M1-M4_top10_radar_${state.exportDpi}dpi.png`)}catch(error){$('radarStatus').textContent=error.message}}
function pdfBlob(jpeg,w,h,widthIn){const enc=new TextEncoder(),parts=[],offsets=[0],push=x=>parts.push(typeof x==='string'?enc.encode(x):x),pageW=widthIn*72,pageH=pageW*h/w;push('%PDF-1.4\n%âãÏÓ\n');let length=parts.reduce((a,b)=>a+b.length,0);const obj=(id,body)=>{offsets[id]=length;const a=enc.encode(`${id} 0 obj\n${body}\nendobj\n`);parts.push(a);length+=a.length};obj(1,'<< /Type /Catalog /Pages 2 0 R >>');obj(2,'<< /Type /Pages /Kids [3 0 R] /Count 1 >>');obj(3,`<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${pageW.toFixed(3)} ${pageH.toFixed(3)}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>`);offsets[4]=length;push(`4 0 obj\n<< /Type /XObject /Subtype /Image /Width ${w} /Height ${h} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpeg.length} >>\nstream\n`);push(jpeg);push('\nendstream\nendobj\n');length=parts.reduce((a,b)=>a+b.length,0);const content=`q\n${pageW.toFixed(3)} 0 0 ${pageH.toFixed(3)} 0 0 cm\n/Im0 Do\nQ`;obj(5,`<< /Length ${content.length} >>\nstream\n${content}\nendstream`);const xref=length;push(`xref\n0 6\n0000000000 65535 f \n${offsets.slice(1).map(o=>String(o).padStart(10,'0')+' 00000 n ').join('\n')}\ntrailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF`);return new Blob(parts,{type:'application/pdf'})}
async function exportPdf(){try{const c=await rasterCanvas(),blob=await new Promise(ok=>c.toBlob(ok,'image/jpeg',.98)),jpeg=new Uint8Array(await blob.arrayBuffer());download(pdfBlob(jpeg,c.width,c.height,state.exportWidthIn),'Macaque_M1-M4_top10_radar.pdf')}catch(error){$('radarStatus').textContent=error.message}}
function exportCsv(){download(new Blob([sourceCsv],{type:'text/csv'}),'Macaque_M1-M4_radar_frozen117_cells.csv')}
function exportLayout(){download(new Blob([JSON.stringify(state,null,2)],{type:'application/json'}),'Macaque_M1-M4_radar_layout.json')}
async function importLayout(event){const file=event.target.files[0];if(!file)return;try{state=mergeState(JSON.parse(await file.text()));classControls();featureControls();colorControls();sync();render()}catch{$('radarStatus').textContent='The selected radar layout could not be read.'}event.target.value=''}

async function init(){
  bindControls();sync();
  sourceCsv=await fetch('data/macaque_m_radar_cells.csv?v=zrange-20260923').then(r=>{if(!r.ok)throw Error(r.status);return r.text()});
  rows=parseCsv(sourceCsv);
  if(rows.length!==117)throw Error(`Expected 117 cells, found ${rows.length}`);
  for(const c of classes)if(rows.filter(r=>r.M_class===c).length!==expected[c])throw Error(`Unexpected ${c} count`);
  render();
}
init().catch(error=>$('radarStatus').textContent='Unable to load Macaque M radar data: '+error.message);
