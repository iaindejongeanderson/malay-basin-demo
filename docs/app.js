/* ------------------------------------------------------------------ *
 * Malay Basin screening cut-off explorer
 * After de Jonge-Anderson et al. (2025), IJGGC 143, 104347.
 *
 * Zones are reclassified in the browser from packed property grids, so
 * any combination of cut-offs can be explored - not just the published one.
 * ------------------------------------------------------------------ */
'use strict';

/* ============================ 1. DATA ============================= */

/* Per-aquifer constants. thickness/NTG are the published petrophysical values
   (Table 2); pubArea/pubCap are what the paper reports at its own cut-offs;
   repro is how closely this tool's reclassification matches the archived
   optimal-zone raster at the published cut-offs. */
const AQUIFERS = [
  { id:'B', age:'Pliocene',      h:162, hSd:12,  ntg:0.17, ntgSd:0.09, pubArea:0,     pubCap:null,             repro:99.82 },
  { id:'D', age:'Late Miocene',  h:287, hSd:262, ntg:0.16, ntgSd:0.09, pubArea:3348,  pubCap:[0.52,0.14,0.02], repro:94.26 },
  { id:'E', age:'Late Miocene',  h:354, hSd:280, ntg:0.27, ntgSd:0.14, pubArea:13894, pubCap:[3.99,1.14,0.18], repro:94.15 },
  { id:'F', age:'Mid Miocene',   h:449, hSd:415, ntg:0.12, ntgSd:0.06, pubArea:18108, pubCap:[2.90,0.84,0.13], repro:95.91 },
  { id:'H', age:'Mid Miocene',   h:393, hSd:294, ntg:0.12, ntgSd:0.10, pubArea:22290, pubCap:[3.93,1.04,0.14], repro:95.46 },
  { id:'I', age:'Early Miocene', h:610, hSd:264, ntg:0.13, ntgSd:0.08, pubArea:24924, pubCap:[5.22,1.67,0.33], repro:98.05 },
  { id:'J', age:'Early Miocene', h:272, hSd:118, ntg:0.42, ntgSd:0.17, pubArea:12898, pubCap:[3.67,1.22,0.27], repro:97.52 },
  { id:'K', age:'Early Miocene', h:383, hSd:176, ntg:0.44, ntgSd:0.13, pubArea:10643, pubCap:[4.28,1.52,0.37], repro:97.87 }
];

/* The published screening cut-offs. */
const PUBLISHED = {
  porOpt:10, porSub:6,          /* % */
  denOpt:300, denSub:100,       /* kg/m3 */
  fauOpt:2, fauSub:1,           /* fault setback class, see FAULT_LABELS */
  needSC:true,                  /* require supercritical CO2 */
  overOpt:true, overSub:false   /* exclude overpressured ground */
};

/* Fault setback classes as carried by the archived fault rasters. */
const FAULT_LABELS = ['any', '>0.5 km', '>2 km', '>10 km'];
const FAULT_LONG = [
  'no setback applied',
  'outside the 0.5 km fault buffer',
  'outside the 2 km fault buffer',
  'outside the 10 km fault aggregation distance'
];

const SWIRR = 0.27, SWIRR_SD = 0.05, EFF = 0.02, EFF_SD = 0.01;
const TRIALS = 10000, SEED = 42;

/* Packed-byte quantisation, must match the exporter. */
const POR_STEP = 0.2, DEN_STEP = 4;

/* ========================= 2. SMALL MATHS ========================= */

function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;
  let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;
  return((t^t>>>14)>>>0)/4294967296;};}

function truncNormal(rng,mu,sd,lo,hi){
  if(!(sd>0))return Math.min(hi,Math.max(lo,mu));
  for(let i=0;i<500;i++){
    const u1=1-rng(),u2=rng();
    const z=Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2);
    const v=mu+sd*z;
    if(v>=lo&&v<=hi)return v;
  }
  return Math.min(hi,Math.max(lo,mu));
}

function pctOf(s,p){if(!s.length)return 0;
  const i=p*(s.length-1),lo=Math.floor(i),hi=Math.min(lo+1,s.length-1);
  return s[lo]+(i-lo)*(s[hi]-s[lo]);}

const fmtGt = v => v===0?'0':v>=100?v.toFixed(0):v>=10?v.toFixed(1):v>=1?v.toFixed(2)
  :v>=0.01?v.toFixed(3):v>=0.0001?v.toFixed(4):v.toExponential(1);
const fmtPct = v => v>0 && v<10 ? v.toFixed(1) : v.toFixed(0);
const fmtInt = v => Math.round(v).toLocaleString('en-GB');

/* ============ 3. PROJECTION: Kertau 1968 / UTM zone 48N ============ */
const PROJ=(function(){
  const a=6377304.063,f=1/300.8017,e2=2*f-f*f,k0=0.9996;
  const lon0=105*Math.PI/180,FE=500000,ep2=e2/(1-e2);
  const e1=(1-Math.sqrt(1-e2))/(1+Math.sqrt(1-e2));
  function M(phi){return a*((1-e2/4-3*e2*e2/64-5*e2*e2*e2/256)*phi
    -(3*e2/8+3*e2*e2/32+45*e2*e2*e2/1024)*Math.sin(2*phi)
    +(15*e2*e2/256+45*e2*e2*e2/1024)*Math.sin(4*phi)
    -(35*e2*e2*e2/3072)*Math.sin(6*phi));}
  function forward(lonDeg,latDeg){
    const phi=latDeg*Math.PI/180,lam=lonDeg*Math.PI/180;
    const N=a/Math.sqrt(1-e2*Math.sin(phi)**2);
    const T=Math.tan(phi)**2,C=ep2*Math.cos(phi)**2,A1=(lam-lon0)*Math.cos(phi);
    return [FE+k0*N*(A1+(1-T+C)*A1**3/6+(5-18*T+T*T+72*C-58*ep2)*A1**5/120),
      k0*(M(phi)+N*Math.tan(phi)*(A1*A1/2+(5-T+9*C+4*C*C)*A1**4/24
        +(61-58*T+T*T+600*C-330*ep2)*A1**6/720))];
  }
  function inverse(E,N){
    const x=E-FE,y=N;
    const mu=(y/k0)/(a*(1-e2/4-3*e2*e2/64-5*e2*e2*e2/256));
    const phi1=mu+(3*e1/2-27*e1**3/32)*Math.sin(2*mu)+(21*e1*e1/16-55*e1**4/32)*Math.sin(4*mu)
      +(151*e1**3/96)*Math.sin(6*mu)+(1097*e1**4/512)*Math.sin(8*mu);
    const C1=ep2*Math.cos(phi1)**2,T1=Math.tan(phi1)**2;
    const N1=a/Math.sqrt(1-e2*Math.sin(phi1)**2);
    const R1=a*(1-e2)/Math.pow(1-e2*Math.sin(phi1)**2,1.5);
    const D=x/(N1*k0);
    const phi=phi1-(N1*Math.tan(phi1)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*ep2)*D**4/24
      +(61+90*T1+298*C1+45*T1*T1-252*ep2-3*C1*C1)*D**6/720);
    const lam=lon0+(D-(1+2*T1+C1)*D**3/6
      +(5-2*C1+28*T1-3*C1*C1+8*ep2+24*T1*T1)*D**5/120)/Math.cos(phi1);
    return [lam*180/Math.PI,phi*180/Math.PI];
  }
  return {forward,inverse};
})();

/* ===================== 4. GRID LOAD / DECODE ====================== */
/* Packed RGB: R = porosity/0.2, G = density/4,
   B bits 0-1 fault class, bit2 gas phase, bit3 overpressured, bit4 valid. */

const grids = {};
const pending = {};

function loadGrid(id){
  if(grids[id]) return Promise.resolve(grids[id]);
  if(pending[id]) return pending[id];
  pending[id] = new Promise(resolve=>{
    const done = () => {
      const m = window.MALAY_CRIT && window.MALAY_CRIT[id];
      if(!m){ resolve(null); return; }
      const img = new Image();
      img.onload = () => {
        const cv = document.createElement('canvas');
        cv.width = m.w; cv.height = m.h;
        const ctx = cv.getContext('2d',{willReadFrequently:true});
        ctx.drawImage(img,0,0);
        const px = ctx.getImageData(0,0,m.w,m.h).data;
        const n = m.w*m.h;
        const por=new Uint8Array(n), den=new Uint8Array(n), flg=new Uint8Array(n);
        let footprint=0;
        for(let i=0;i<n;i++){
          por[i]=px[i*4]; den[i]=px[i*4+1]; flg[i]=px[i*4+2];
          if(flg[i]&16) footprint++;
        }
        const g={id,w:m.w,h:m.h,originX:m.originX,originY:m.originY,cell:m.cell,
                 por,den,flg,footprint,cellKm2:m.cell*m.cell/1e6,
                 cls:new Uint8Array(n),base:null,baseKey:null};
        grids[id]=g; resolve(g);
      };
      img.onerror = () => resolve(null);
      img.src = m.png;
    };
    if(window.MALAY_CRIT && window.MALAY_CRIT[id]) { done(); return; }
    const s=document.createElement('script');
    s.src='data/crit_'+id+'.js';
    s.onload=done;
    s.onerror=()=>resolve(null);
    document.head.appendChild(s);
  });
  return pending[id];
}

/* ==================== 5. THE CLASSIFICATION ======================= */
/* This is the whole point of the app: one pass, run on every control change. */

function classify(g,c){
  const pOpt=Math.round(c.porOpt/POR_STEP), pSub=Math.round(c.porSub/POR_STEP);
  const dOpt=Math.round(c.denOpt/DEN_STEP), dSub=Math.round(c.denSub/DEN_STEP);
  const por=g.por, den=g.den, flg=g.flg, cls=g.cls, n=g.w*g.h;

  let nOpt=0,nSub=0,nNon=0;
  let sPor=0,sPor2=0,sDen=0,sDen2=0;

  for(let i=0;i<n;i++){
    const b=flg[i];
    if(!(b&16)){ cls[i]=0; continue; }
    const p=por[i], d=den[i], fc=b&3;
    const sc = !c.needSC || !(b&4);
    const over = (b&8)!==0;

    if(p>=pOpt && d>=dOpt && fc>=c.fauOpt && sc && !(c.overOpt&&over)){
      cls[i]=3; nOpt++;
      const pv=p*POR_STEP, dv=d*DEN_STEP;
      sPor+=pv; sPor2+=pv*pv; sDen+=dv; sDen2+=dv*dv;
    } else if(p>=pSub && d>=dSub && fc>=c.fauSub && sc && !(c.overSub&&over)){
      cls[i]=2; nSub++;
    } else { cls[i]=1; nNon++; }
  }

  const k=g.cellKm2;
  const mPor = nOpt? sPor/nOpt : 0, mDen = nOpt? sDen/nOpt : 0;
  return {
    nOpt,nSub,nNon,
    optArea:nOpt*k, subArea:nSub*k, nonArea:nNon*k, footArea:g.footprint*k,
    por: mPor, porSd: nOpt? Math.sqrt(Math.max(0,sPor2/nOpt-mPor*mPor)) : 0,
    den: mDen, denSd: nOpt? Math.sqrt(Math.max(0,sDen2/nOpt-mDen*mDen)) : 0
  };
}

/* Area that would be optimal if only one criterion were dropped / applied alone. */
function breakdown(g,c){
  const pOpt=Math.round(c.porOpt/POR_STEP), dOpt=Math.round(c.denOpt/DEN_STEP);
  const por=g.por,den=g.den,flg=g.flg,n=g.w*g.h;
  let fPor=0,fDen=0,fFau=0,fPha=0,fOver=0;
  for(let i=0;i<n;i++){
    const b=flg[i];
    if(!(b&16)) continue;
    if(por[i]<pOpt) fPor++;
    if(den[i]<dOpt) fDen++;
    if((b&3)<c.fauOpt) fFau++;
    if(c.needSC && (b&4)) fPha++;
    if(c.overOpt && (b&8)) fOver++;
  }
  const k=g.cellKm2;
  return [
    {key:'por',  label:'Porosity',      area:fPor*k},
    {key:'den',  label:'CO₂ density', area:fDen*k},
    {key:'fau',  label:'Fault setback', area:fFau*k},
    {key:'pha',  label:'CO₂ phase',   area:fPha*k},
    {key:'over', label:'Overpressure',  area:fOver*k}
  ];
}

/* ======================= 6. CAPACITY ============================== */

function capacity(stats,a,trials,seed){
  if(!stats.optArea) return {p10:0,p50:0,p90:0,mean:0,samples:null,sorted:[]};
  const rng=mulberry32(seed>>>0);
  const out=new Float64Array(trials);
  const A=stats.optArea*1e6;
  for(let i=0;i<trials;i++){
    const h=truncNormal(rng,a.h,a.hSd,0,Infinity);
    const ntg=truncNormal(rng,a.ntg,a.ntgSd,0,1);
    const por=truncNormal(rng,stats.por,stats.porSd,0,100)/100;
    const sw=truncNormal(rng,SWIRR,SWIRR_SD,0,1);
    const e=truncNormal(rng,EFF,EFF_SD,0,1);
    const rho=truncNormal(rng,stats.den,stats.denSd,0,Infinity);
    out[i]=A*h*ntg*por*(1-sw)*e*rho/1e12;
  }
  const sorted=Float64Array.from(out).sort();
  let s=0; for(let i=0;i<trials;i++) s+=out[i];
  return {p10:pctOf(sorted,0.90),p50:pctOf(sorted,0.50),p90:pctOf(sorted,0.10),
          mean:s/trials,samples:out,sorted};
}

/* ========================== 7. STATE ============================== */

const state={
  aquifer:'I',
  cuts:Object.assign({},PUBLISHED),
  sweep:'por',
  grid:null, stats:null, cap:null, bd:null,
  cam:null
};

const isPublished = () => Object.keys(PUBLISHED).every(k=>state.cuts[k]===PUBLISHED[k]);

/* ===================== 8. CANVAS PLUMBING ======================== */

function cssVar(n){return getComputedStyle(document.documentElement).getPropertyValue(n).trim();}
function hexToRgb(h){h=h.replace('#','');if(h.length===3)h=h.split('').map(c=>c+c).join('');
  return [parseInt(h.slice(0,2),16),parseInt(h.slice(2,4),16),parseInt(h.slice(4,6),16)];}

function fitCanvas(cv,cssH){
  const dpr=window.devicePixelRatio||1;
  const w=cv.clientWidth||cv.parentElement.clientWidth;
  const h=cssH||cv.clientHeight;
  cv.width=Math.max(1,Math.round(w*dpr));
  cv.height=Math.max(1,Math.round(h*dpr));
  cv.style.height=h+'px';
  const ctx=cv.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.clearRect(0,0,w,h);
  return {ctx,w,h};
}

function roundTopRect(ctx,x,y,w,h,r){
  r=Math.max(0,Math.min(r,w/2,h));
  ctx.beginPath();
  ctx.moveTo(x,y+h);ctx.lineTo(x,y+r);ctx.quadraticCurveTo(x,y,x+r,y);
  ctx.lineTo(x+w-r,y);ctx.quadraticCurveTo(x+w,y,x+w,y+r);ctx.lineTo(x+w,y+h);
  ctx.closePath();ctx.fill();
}
function roundRightRect(ctx,x,y,w,h,r){
  r=Math.max(0,Math.min(r,w,h/2));
  ctx.beginPath();
  ctx.moveTo(x,y);ctx.lineTo(x+w-r,y);ctx.quadraticCurveTo(x+w,y,x+w,y+r);
  ctx.lineTo(x+w,y+h-r);ctx.quadraticCurveTo(x+w,y+h,x+w-r,y+h);
  ctx.lineTo(x,y+h);ctx.closePath();ctx.fill();
}

function bindTip(cv,tipEl,hitFn){
  cv.addEventListener('mousemove',e=>{
    const r=cv.getBoundingClientRect();
    const hit=hitFn(e.clientX-r.left,e.clientY-r.top,r.width,r.height);
    if(!hit){tipEl.classList.remove('on');return;}
    tipEl.innerHTML=hit.html;tipEl.classList.add('on');
    const tw=tipEl.offsetWidth,th=tipEl.offsetHeight;
    let x=hit.x-tw/2,y=hit.y-th-10;
    x=Math.max(2,Math.min(r.width-tw-2,x));
    if(y<2)y=hit.y+14;
    tipEl.style.left=x+'px';tipEl.style.top=y+'px';
  });
  cv.addEventListener('mouseleave',()=>tipEl.classList.remove('on'));
}

/* ========================== 9. MAP =============================== */

const mapCv=document.getElementById('map');
const readoutEl=document.getElementById('readout');

function resetCam(){
  const g=state.grid; if(!g)return;
  const holder=mapCv.parentElement;
  const W=g.w*g.cell, H=g.h*g.cell;
  const k=Math.min(holder.clientWidth/W,holder.clientHeight/H)*0.94;
  state.cam={cx:g.originX+W/2, cy:g.originY-H/2, k};
}
function w2s(x,y,w,h){const c=state.cam;return [w/2+(x-c.cx)*c.k, h/2-(y-c.cy)*c.k];}
function s2w(sx,sy,w,h){const c=state.cam;return [c.cx+(sx-w/2)/c.k, c.cy-(sy-h/2)/c.k];}

function buildBase(g){
  const key=(document.documentElement.getAttribute('data-theme')||'auto')+
    (window.matchMedia('(prefers-color-scheme: dark)').matches?'d':'l')+'|'+JSON.stringify(state.cuts);
  if(g.base&&g.baseKey===key) return g.base;
  const col={1:hexToRgb(cssVar('--zone-nonviable')),
             2:hexToRgb(cssVar('--zone-suboptimal')),
             3:hexToRgb(cssVar('--zone-optimal'))};
  const cv=g.base||document.createElement('canvas');
  cv.width=g.w;cv.height=g.h;
  const ctx=cv.getContext('2d');
  const img=ctx.createImageData(g.w,g.h);
  const d=img.data;
  for(let i=0,n=g.w*g.h;i<n;i++){
    const c=g.cls[i];
    if(!c){d[i*4+3]=0;continue;}
    const rgb=col[c];
    d[i*4]=rgb[0];d[i*4+1]=rgb[1];d[i*4+2]=rgb[2];d[i*4+3]=255;
  }
  ctx.putImageData(img,0,0);
  g.base=cv;g.baseKey=key;
  return cv;
}

function drawMap(){
  const holder=mapCv.parentElement;
  const cssW=holder.clientWidth, cssH=holder.clientHeight;
  const dpr=window.devicePixelRatio||1;
  mapCv.width=Math.round(cssW*dpr);mapCv.height=Math.round(cssH*dpr);
  const ctx=mapCv.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  ctx.clearRect(0,0,cssW,cssH);
  const g=state.grid; if(!g||!state.cam) return;

  const [sx0,sy0]=w2s(g.originX,g.originY,cssW,cssH);
  const dw=g.w*g.cell*state.cam.k, dh=g.h*g.cell*state.cam.k;
  ctx.imageSmoothingEnabled = state.cam.k*g.cell < 1;
  ctx.drawImage(buildBase(g),sx0,sy0,dw,dh);

  drawGraticule(ctx,cssW,cssH);
  drawScaleBar(ctx,cssW,cssH);
}

function drawGraticule(ctx,w,h){
  const [lonA,latA]=PROJ.inverse(...s2w(0,h,w,h));
  const [lonB,latB]=PROJ.inverse(...s2w(w,0,w,h));
  const span=Math.max(lonB-lonA,latB-latA);
  const step=span>4?1:span>2?0.5:span>0.8?0.25:0.1;
  ctx.save();
  ctx.strokeStyle=cssVar('--grid');ctx.fillStyle=cssVar('--muted');
  ctx.lineWidth=1;ctx.font='10px system-ui, sans-serif';
  for(let lon=Math.ceil(lonA/step)*step;lon<=lonB;lon+=step){
    ctx.beginPath();
    for(let i=0;i<=20;i++){
      const lat=latA+(latB-latA)*i/20;
      const [sx,sy]=w2s(...PROJ.forward(lon,lat),w,h);
      i?ctx.lineTo(sx,sy):ctx.moveTo(sx,sy);
    }
    ctx.stroke();
    const [lx,ly]=w2s(...PROJ.forward(lon,latB),w,h);
    ctx.fillText(lon.toFixed(step<1?2:0)+'°E',lx+3,Math.max(ly+11,11));
  }
  for(let lat=Math.ceil(latA/step)*step;lat<=latB;lat+=step){
    ctx.beginPath();
    for(let i=0;i<=20;i++){
      const lon=lonA+(lonB-lonA)*i/20;
      const [sx,sy]=w2s(...PROJ.forward(lon,lat),w,h);
      i?ctx.lineTo(sx,sy):ctx.moveTo(sx,sy);
    }
    ctx.stroke();
    const [lx,ly]=w2s(...PROJ.forward(lonA,lat),w,h);
    ctx.fillText(lat.toFixed(step<1?2:0)+'°N',Math.max(lx+3,4),ly-3);
  }
  ctx.restore();
}

function drawScaleBar(ctx,w,h){
  const metres=Math.min(140,w*0.3)/state.cam.k;
  const nice=[1e3,2e3,5e3,1e4,2e4,5e4,1e5,2e5];
  let pick=nice[0]; for(const n of nice) if(n<=metres) pick=n;
  const px=pick*state.cam.k, x=w-px-14, y=h-16;
  ctx.save();
  ctx.strokeStyle=cssVar('--text-primary');ctx.fillStyle=cssVar('--text-secondary');
  ctx.lineWidth=2;
  ctx.beginPath();ctx.moveTo(x,y-4);ctx.lineTo(x,y);ctx.lineTo(x+px,y);ctx.lineTo(x+px,y-4);ctx.stroke();
  ctx.font='11px system-ui, sans-serif';ctx.textAlign='center';
  ctx.fillText((pick/1000)+' km',x+px/2,y-7);
  ctx.restore();
}

/* map interaction */
let drag=null;
mapCv.addEventListener('pointerdown',e=>{
  mapCv.setPointerCapture(e.pointerId);
  const r=mapCv.getBoundingClientRect();
  drag={sx:e.clientX-r.left,sy:e.clientY-r.top,cx:state.cam.cx,cy:state.cam.cy};
});
mapCv.addEventListener('pointermove',e=>{
  const r=mapCv.getBoundingClientRect();
  const sx=e.clientX-r.left, sy=e.clientY-r.top;
  if(drag){
    state.cam.cx=drag.cx-(sx-drag.sx)/state.cam.k;
    state.cam.cy=drag.cy+(sy-drag.sy)/state.cam.k;
    drawMap(); return;
  }
  updateReadout(sx,sy,r.width,r.height);
});
mapCv.addEventListener('pointerup',()=>{drag=null;});
mapCv.addEventListener('pointerleave',()=>{drag=null;renderReadout();});
mapCv.addEventListener('wheel',e=>{
  e.preventDefault();
  const r=mapCv.getBoundingClientRect();
  const sx=e.clientX-r.left, sy=e.clientY-r.top;
  const [wx,wy]=s2w(sx,sy,r.width,r.height);
  zoomAbout(wx,wy,Math.exp(-e.deltaY*0.0015),sx,sy,r.width,r.height);
},{passive:false});

function zoomAbout(wx,wy,f,sx,sy,w,h){
  const g=state.grid,c=state.cam;
  const base=Math.min(w/(g.w*g.cell),h/(g.h*g.cell))*0.94;
  c.k=Math.max(base*0.8,Math.min(base*40,c.k*f));
  c.cx=wx-(sx-w/2)/c.k;
  c.cy=wy+(sy-h/2)/c.k;
  drawMap();
}

const ZONE_NAME={0:'outside aquifer',1:'non-viable',2:'sub-optimal',3:'optimal'};

function renderReadout(){
  const s=state.stats; if(!s){readoutEl.innerHTML='';return;}
  readoutEl.innerHTML =
    `<span><b>${fmtInt(s.optArea)} km&#178;</b> optimal</span>` +
    `<span>${fmtInt(s.subArea)} km&#178; sub-optimal</span>` +
    `<span>${fmtInt(s.footArea)} km&#178; aquifer footprint</span>`;
}

function updateReadout(sx,sy,w,h){
  const g=state.grid; if(!g) return;
  const [wx,wy]=s2w(sx,sy,w,h);
  const px=Math.floor((wx-g.originX)/g.cell), py=Math.floor((g.originY-wy)/g.cell);
  renderReadout();
  if(px<0||py<0||px>=g.w||py>=g.h) return;
  const i=py*g.w+px, b=g.flg[i];
  const [lon,lat]=PROJ.inverse(wx,wy);
  let s=`<span>${lat.toFixed(2)}°N ${lon.toFixed(2)}°E</span>`;
  if(b&16){
    s+=`<span><b>${ZONE_NAME[g.cls[i]]}</b></span>`+
       `<span>&phi; ${(g.por[i]*POR_STEP).toFixed(1)} % &middot; `+
       `&rho; ${fmtInt(g.den[i]*DEN_STEP)} kg/m&#179; &middot; `+
       `${FAULT_LABELS[b&3]} &middot; ${(b&4)?'gas':'supercritical'}`+
       `${(b&8)?' &middot; overpressured':''}</span>`;
  } else {
    s+=`<span>outside aquifer</span>`;
  }
  readoutEl.innerHTML=readoutEl.innerHTML+s;
}

/* ====================== 10. CUT-OFF CONTROLS ====================== */

const CRITERIA=[
  {key:'por', name:'Porosity', unit:'%', min:0, max:40, step:0.5,
   optKey:'porOpt', subKey:'porSub', pub:'≥ 10 % / ≥ 6 %'},
  {key:'den', name:'CO₂ density', unit:'kg/m³', min:0, max:900, step:10,
   optKey:'denOpt', subKey:'denSub', pub:'≥ 300 / ≥ 100 kg/m³'}
];

function renderControls(){
  const host=document.getElementById('critControls');
  host.innerHTML='';

  for(const c of CRITERIA){
    const d=document.createElement('div');
    d.className='crit';
    d.innerHTML=
      `<div class="crit-hd"><span class="nm">${c.name}</span><span class="pub">published ${c.pub}</span></div>`+
      row('Optimal',   c, c.optKey, '--zone-optimal')+
      row('Sub-opt.',  c, c.subKey, '--zone-suboptimal');
    host.appendChild(d);
  }

  /* fault setback */
  const fd=document.createElement('div');
  fd.className='crit';
  fd.innerHTML=
    `<div class="crit-hd"><span class="nm">Fault setback</span><span class="pub">published &gt;2 km / &gt;0.5 km</span></div>`+
    segRow('Optimal','fauOpt','--zone-optimal')+
    segRow('Sub-opt.','fauSub','--zone-suboptimal');
  host.appendChild(fd);

  /* phase + overpressure */
  const od=document.createElement('div');
  od.className='crit';
  od.innerHTML=
    `<div class="crit-hd"><span class="nm">Phase &amp; pressure</span><span class="pub">published as ticked</span></div>`+
    `<label class="chk"><input type="checkbox" data-flag="needSC" ${state.cuts.needSC?'checked':''}>
       Require supercritical CO&#8322; (both classes)</label>`+
    `<label class="chk"><input type="checkbox" data-flag="overOpt" ${state.cuts.overOpt?'checked':''}>
       Exclude overpressured ground from <b>optimal</b></label>`+
    `<label class="chk"><input type="checkbox" data-flag="overSub" ${state.cuts.overSub?'checked':''}>
       Exclude overpressured ground from <b>sub-optimal</b></label>`;
  host.appendChild(od);

  host.querySelectorAll('input[type="range"],input[type="number"]').forEach(el=>{
    el.addEventListener('input',()=>{
      const k=el.dataset.k;
      let v=parseFloat(el.value);
      if(!isFinite(v))return;
      state.cuts[k]=v;
      /* sub-optimal must stay no stricter than optimal */
      if(k==='porOpt'&&state.cuts.porSub>v) state.cuts.porSub=v;
      if(k==='porSub'&&v>state.cuts.porOpt) state.cuts.porOpt=v;
      if(k==='denOpt'&&state.cuts.denSub>v) state.cuts.denSub=v;
      if(k==='denSub'&&v>state.cuts.denOpt) state.cuts.denOpt=v;
      syncControlValues();
      recompute(true);
    });
  });
  host.querySelectorAll('[data-fault]').forEach(b=>{
    b.addEventListener('click',()=>{
      const k=b.dataset.faultKey, v=parseInt(b.dataset.fault,10);
      state.cuts[k]=v;
      if(k==='fauOpt'&&state.cuts.fauSub>v) state.cuts.fauSub=v;
      if(k==='fauSub'&&v>state.cuts.fauOpt) state.cuts.fauOpt=v;
      renderControls(); recompute();
    });
  });
  host.querySelectorAll('[data-flag]').forEach(cb=>{
    cb.addEventListener('change',()=>{ state.cuts[cb.dataset.flag]=cb.checked; recompute(); });
  });

  function row(label,c,key,swatch){
    return `<div class="crow">
      <span class="lbl"><span class="sw" style="background:var(${swatch})"></span>${label}</span>
      <input type="range" data-k="${key}" min="${c.min}" max="${c.max}" step="${c.step}" value="${state.cuts[key]}">
      <input type="number" data-k="${key}" min="${c.min}" max="${c.max}" step="${c.step}" value="${state.cuts[key]}">
    </div>`;
  }
  function segRow(label,key,swatch){
    let b='';
    for(let i=0;i<4;i++){
      b+=`<button type="button" data-fault="${i}" data-fault-key="${key}"
           aria-pressed="${state.cuts[key]===i}" title="${FAULT_LONG[i]}">${FAULT_LABELS[i]}</button>`;
    }
    return `<div class="crow" style="grid-template-columns:76px 1fr">
      <span class="lbl"><span class="sw" style="background:var(${swatch})"></span>${label}</span>
      <span class="seg">${b}</span></div>`;
  }
}

/* keep paired range/number inputs in step without a full re-render */
function syncControlValues(){
  document.querySelectorAll('#critControls input[data-k]').forEach(el=>{
    const v=state.cuts[el.dataset.k];
    if(parseFloat(el.value)!==v) el.value=v;
  });
}

/* ========================== 11. CHARTS =========================== */

let bdHit=()=>null, histHit=()=>null, swHit=()=>null;

function drawBreakdown(){
  const cv=document.getElementById('breakdown');
  const {ctx,w,h}=fitCanvas(cv,190);
  const bd=state.bd, st=state.stats;
  if(!bd||!st){return;}
  const rows=bd.slice().sort((a,b)=>b.area-a.area);
  const foot=st.footArea||1;
  const padL=104,padR=52,padT=6,padB=24;
  const plotW=w-padL-padR;
  const rowH=(h-padT-padB)/rows.length;
  const max=Math.max(foot*0.02,...rows.map(r=>r.area));

  ctx.font='11px system-ui, sans-serif';
  const bars=[];
  for(let i=0;i<rows.length;i++){
    const r=rows[i];
    const y=padT+i*rowH+rowH*0.2, bh=rowH*0.6;
    const bw=(r.area/max)*plotW;
    ctx.fillStyle=cssVar('--series-1');
    roundRightRect(ctx,padL,y,Math.max(1,bw),bh,3);
    ctx.fillStyle=cssVar('--text-secondary');
    ctx.textAlign='right';ctx.textBaseline='middle';
    ctx.fillText(r.label,padL-9,y+bh/2);
    ctx.textAlign='left';
    ctx.fillStyle=cssVar('--muted');
    ctx.fillText(fmtPct(r.area/foot*100)+' %',padL+bw+7,y+bh/2);
    bars.push({r,y,bh,bw});
  }
  ctx.strokeStyle=cssVar('--axis');ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(padL+0.5,padT);ctx.lineTo(padL+0.5,h-padB+2);ctx.stroke();
  ctx.fillStyle=cssVar('--text-secondary');
  ctx.textAlign='left';ctx.textBaseline='top';
  ctx.fillText('Rejected area, one criterion at a time — footprint '+fmtInt(foot)+' km²',padL,h-padB+7);

  bdHit=(mx,my)=>{
    for(const b of bars){
      if(my>=b.y&&my<=b.y+b.bh&&mx>=padL&&mx<=padL+Math.max(6,b.bw)){
        return {x:padL+b.bw/2,y:b.y,
          html:`<b>${b.r.label}</b><br>rejects ${fmtInt(b.r.area)} km² `+
               `(${(b.r.area/foot*100).toFixed(1)} % of footprint)`};
      }
    }
    return null;
  };
}

function drawHistogram(){
  const cv=document.getElementById('hist');
  const {ctx,w,h}=fitCanvas(cv,180);
  const res=state.cap;
  if(!res||!res.samples){
    ctx.fillStyle=cssVar('--muted');ctx.font='13px system-ui, sans-serif';
    ctx.textAlign='center';ctx.textBaseline='middle';
    ctx.fillText('No optimal zone at these cut-offs',w/2,h/2);
    histHit=()=>null;return;
  }
  const padL=40,padR=10,padT=14,padB=28;
  const plotW=w-padL-padR,plotH=h-padT-padB;
  const hi=Math.max(res.p10*1.35,pctOf(res.sorted,0.985))||1;
  const NB=40,bw=hi/NB;
  const bins=new Array(NB).fill(0);
  for(let i=0;i<res.samples.length;i++){
    const b=Math.floor(res.samples[i]/bw);
    if(b>=0&&b<NB)bins[b]++;
  }
  const maxC=Math.max(...bins)||1;
  const X=v=>padL+(v/hi)*plotW, Y=c=>padT+plotH-(c/maxC)*plotH;

  ctx.strokeStyle=cssVar('--grid');ctx.lineWidth=1;
  ctx.fillStyle=cssVar('--muted');ctx.font='10px system-ui, sans-serif';
  ctx.textAlign='right';ctx.textBaseline='middle';
  for(let i=0;i<=3;i++){
    const c=maxC*i/3,y=Math.round(Y(c))+0.5;
    ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(w-padR,y);ctx.stroke();
    ctx.fillText((c/TRIALS*100).toFixed(0)+'%',padL-6,y);
  }
  ctx.fillStyle=cssVar('--series-1');
  const pxw=plotW/NB;
  for(let i=0;i<NB;i++){
    if(!bins[i])continue;
    const y=Y(bins[i]);
    roundTopRect(ctx,padL+i*pxw+1,y,Math.max(1,pxw-2),padT+plotH-y,3);
  }
  ctx.strokeStyle=cssVar('--axis');
  ctx.beginPath();ctx.moveTo(padL,padT+plotH+0.5);ctx.lineTo(w-padR,padT+plotH+0.5);ctx.stroke();
  ctx.fillStyle=cssVar('--muted');ctx.textAlign='center';ctx.textBaseline='top';
  for(let i=0;i<=4;i++) ctx.fillText(fmtGt(hi*i/4),X(hi*i/4),padT+plotH+5);
  ctx.fillStyle=cssVar('--text-secondary');
  ctx.fillText('Storage capacity (Gt)',padL+plotW/2,padT+plotH+16);

  ctx.save();
  ctx.strokeStyle=cssVar('--text-primary');ctx.fillStyle=cssVar('--text-primary');
  ctx.lineWidth=1.5;ctx.font='600 10px system-ui, sans-serif';ctx.textBaseline='alphabetic';
  for(const m of [{v:res.p90,t:'P90'},{v:res.p50,t:'P50'},{v:res.p10,t:'P10'}]){
    if(m.v>hi)continue;
    const x=Math.round(X(m.v))+0.5;
    ctx.beginPath();ctx.moveTo(x,padT-3);ctx.lineTo(x,padT+plotH);ctx.stroke();
    const right=x>w-padR-34;
    ctx.textAlign=right?'right':'left';
    ctx.fillText(m.t,x+(right?-3:3),padT+4);
  }
  ctx.restore();

  histHit=(mx,my)=>{
    if(mx<padL||mx>w-padR||my<padT||my>padT+plotH)return null;
    const i=Math.floor((mx-padL)/pxw);
    if(i<0||i>=NB||!bins[i])return null;
    return {x:padL+(i+0.5)*pxw,y:Y(bins[i]),
      html:`<b>${fmtGt(i*bw)} – ${fmtGt((i+1)*bw)} Gt</b><br>${(bins[i]/TRIALS*100).toFixed(1)} % of trials`};
  };
}

/* sweep: recompute area + capacity across one criterion's range */
function computeSweep(){
  const g=state.grid,a=AQUIFERS.find(x=>x.id===state.aquifer);
  if(!g)return null;
  const key=state.sweep;
  let vals;
  if(key==='por') vals=Array.from({length:21},(_,i)=>i*2);            /* 0-40 % */
  else if(key==='den') vals=Array.from({length:19},(_,i)=>i*50);      /* 0-900 */
  else vals=[0,1,2,3];

  const saveCls=g.cls;
  g.cls=new Uint8Array(g.w*g.h);                 /* scratch, keep the live one intact */
  const pts=[];
  for(const v of vals){
    const c=Object.assign({},state.cuts);
    if(key==='por'){ c.porOpt=v; if(c.porSub>v)c.porSub=v; }
    else if(key==='den'){ c.denOpt=v; if(c.denSub>v)c.denSub=v; }
    else { c.fauOpt=v; if(c.fauSub>v)c.fauSub=v; }
    const st=classify(g,c);
    const cap=capacity(st,a,2000,SEED);
    pts.push({v,area:st.optArea,cap:cap.p50});
  }
  g.cls=saveCls;
  return {key,pts};
}

function drawSweep(){
  const sw=state.sweepData;
  const areaCv=document.getElementById('sweepArea');
  const capCv=document.getElementById('sweepCap');
  if(!sw){fitCanvas(areaCv,150);fitCanvas(capCv,150);swHit=()=>null;return;}

  const cur = state.sweep==='por'?state.cuts.porOpt:state.sweep==='den'?state.cuts.denOpt:state.cuts.fauOpt;
  const xlabel = state.sweep==='por'?'Optimal porosity cut-off (%)'
    : state.sweep==='den'?'Optimal CO₂ density cut-off (kg/m³)'
    : 'Optimal fault setback class';
  const hits=[];

  panel(areaCv,'Optimal area (km²)', p=>p.area, v=>fmtInt(v), false);
  panel(capCv, 'Capacity P50 (Gt)',      p=>p.cap,  v=>fmtGt(v),  true);

  function panel(cv,title,get,fmt,showX){
    const {ctx,w,h}=fitCanvas(cv,150);
    const padL=52,padR=12,padT=20,padB=showX?30:12;
    const plotW=w-padL-padR,plotH=h-padT-padB;
    const max=Math.max(...sw.pts.map(get))||1;
    const vmin=sw.pts[0].v, vmax=sw.pts[sw.pts.length-1].v;
    const X=v=>padL+((v-vmin)/((vmax-vmin)||1))*plotW;
    const Y=y=>padT+plotH-(y/max)*plotH;

    ctx.strokeStyle=cssVar('--grid');ctx.lineWidth=1;
    ctx.fillStyle=cssVar('--muted');ctx.font='10px system-ui, sans-serif';
    ctx.textAlign='right';ctx.textBaseline='middle';
    for(let i=0;i<=3;i++){
      const v=max*i/3,y=Math.round(Y(v))+0.5;
      ctx.beginPath();ctx.moveTo(padL,y);ctx.lineTo(w-padR,y);ctx.stroke();
      ctx.fillText(fmt(v),padL-6,y);
    }
    ctx.fillStyle=cssVar('--text-secondary');
    ctx.textAlign='left';ctx.textBaseline='top';
    ctx.font='11px system-ui, sans-serif';
    ctx.fillText(title,padL,2);

    /* marker for the value currently set */
    const xc=Math.round(X(cur))+0.5;
    if(cur>=vmin&&cur<=vmax){
      ctx.strokeStyle=cssVar('--series-2');ctx.lineWidth=1.5;
      ctx.beginPath();ctx.moveTo(xc,padT);ctx.lineTo(xc,padT+plotH);ctx.stroke();
    }

    ctx.strokeStyle=cssVar('--series-1');ctx.lineWidth=2;
    ctx.lineJoin='round';ctx.beginPath();
    sw.pts.forEach((p,i)=>{const x=X(p.v),y=Y(get(p));i?ctx.lineTo(x,y):ctx.moveTo(x,y);});
    ctx.stroke();

    ctx.fillStyle=cssVar('--series-1');
    for(const p of sw.pts){
      const x=X(p.v),y=Y(get(p));
      ctx.beginPath();ctx.arc(x,y,2.5,0,Math.PI*2);ctx.fill();
      if(showX) hits.push({x,y,p});
    }

    ctx.strokeStyle=cssVar('--axis');ctx.lineWidth=1;
    ctx.beginPath();ctx.moveTo(padL,padT+plotH+0.5);ctx.lineTo(w-padR,padT+plotH+0.5);ctx.stroke();
    if(showX){
      ctx.fillStyle=cssVar('--muted');ctx.textAlign='center';ctx.textBaseline='top';
      const ticks = state.sweep==='fau'? sw.pts.map(p=>p.v) : sw.pts.filter((_,i)=>i%4===0).map(p=>p.v);
      for(const t of ticks)
        ctx.fillText(state.sweep==='fau'?FAULT_LABELS[t]:String(t),X(t),padT+plotH+5);
      ctx.fillStyle=cssVar('--text-secondary');
      ctx.fillText(xlabel,padL+plotW/2,padT+plotH+17);
    }
  }

  swHit=(mx,my)=>{
    let best=null,bd=1e9;
    for(const hh of hits){
      const d=Math.abs(hh.x-mx);
      if(d<bd&&d<14){bd=d;best=hh;}
    }
    if(!best)return null;
    const lab = state.sweep==='fau'?FAULT_LABELS[best.p.v]
      : best.p.v+(state.sweep==='por'?' %':' kg/m³');
    return {x:best.x,y:best.y,
      html:`<b>cut-off ${lab}</b><br>${fmtInt(best.p.area)} km² optimal<br>${fmtGt(best.p.cap)} Gt (P50)`};
  };
}

/* ======================== 12. RENDERING ========================== */

function renderStats(){
  const st=state.stats,cap=state.cap,a=AQUIFERS.find(x=>x.id===state.aquifer);
  const el=document.getElementById('stats');
  if(!st){el.innerHTML='';return;}
  const dArea = a.pubArea? (st.optArea-a.pubArea)/a.pubArea*100 : null;
  const dCap = (a.pubCap&&cap.p50)? (cap.p50-a.pubCap[1])/a.pubCap[1]*100 : null;
  const sign=v=>(v>=0?'+':'−')+Math.abs(v).toFixed(0)+' %';
  el.innerHTML=
    `<div class="stat hero"><div class="k">Optimal area</div>
      <div class="v">${fmtInt(st.optArea)}<span class="u"> km&#178;</span></div>
      <div class="d">${dArea===null?'no optimal zone published':sign(dArea)+' vs published '+fmtInt(a.pubArea)}</div></div>
     <div class="stat"><div class="k">Sub-optimal</div>
      <div class="v">${fmtInt(st.subArea)}<span class="u"> km&#178;</span></div>
      <div class="d">${fmtPct(st.subArea/st.footArea*100)} % of footprint</div></div>
     <div class="stat"><div class="k">Capacity P50</div>
      <div class="v">${fmtGt(cap.p50)}<span class="u"> Gt</span></div>
      <div class="d">${dCap===null?'&mdash;':sign(dCap)+' vs published '+a.pubCap[1]}</div></div>
     <div class="stat"><div class="k">P10 &ndash; P90</div>
      <div class="v" style="font-size:17px">${fmtGt(cap.p10)}&ndash;${fmtGt(cap.p90)}</div>
      <div class="d">high to low estimate (Gt)</div></div>`;

  const badge=document.getElementById('pubBadge');
  badge.className='badge '+(isPublished()?'pub':'mod');
  badge.textContent=isPublished()?'published cut-offs':'modified cut-offs';
}

function renderCapInputs(){
  const st=state.stats,a=AQUIFERS.find(x=>x.id===state.aquifer);
  const t=document.getElementById('capInputs');
  if(!st){t.innerHTML='';return;}
  t.innerHTML=
    `<thead><tr><th>Term</th><th>Value</th><th>Std. dev.</th><th>Source</th></tr></thead><tbody>
     <tr><td>Area, A</td><td>${fmtInt(st.optArea)} km&#178;</td><td>&mdash;</td><td>live from map</td></tr>
     <tr><td>Porosity, &phi;</td><td>${st.por.toFixed(1)} %</td><td>${st.porSd.toFixed(1)}</td><td>live, within optimal zone</td></tr>
     <tr><td>CO&#8322; density, &rho;</td><td>${fmtInt(st.den)} kg/m&#179;</td><td>${fmtInt(st.denSd)}</td><td>live, within optimal zone</td></tr>
     <tr><td>Thickness, h</td><td>${fmtInt(a.h)} m</td><td>${fmtInt(a.hSd)}</td><td>published, well petrophysics</td></tr>
     <tr><td>Net-to-gross</td><td>${a.ntg.toFixed(2)}</td><td>${a.ntgSd.toFixed(2)}</td><td>published, well petrophysics</td></tr>
     <tr><td>S<sub>wirr</sub></td><td>${SWIRR}</td><td>${SWIRR_SD}</td><td>literature; SD assumed</td></tr>
     <tr><td>Efficiency, E</td><td>${EFF}</td><td>${EFF_SD}</td><td>literature; SD fitted</td></tr>
     </tbody>`;
  document.getElementById('capInfo').textContent=fmtInt(TRIALS)+' trials';
}

function renderTabs(){
  const tabs=document.querySelector('.tabs');
  tabs.querySelectorAll('.tab').forEach(t=>t.remove());
  for(const a of AQUIFERS){
    const b=document.createElement('button');
    b.className='tab';b.type='button';b.setAttribute('role','tab');
    b.setAttribute('aria-selected',String(a.id===state.aquifer));
    b.innerHTML=`<span class="dot${a.pubArea?'':' none'}"></span>Group ${a.id}`;
    b.title=a.age+(a.pubArea?'':' — no optimal zones at the published cut-offs');
    b.addEventListener('click',()=>selectAquifer(a.id));
    tabs.appendChild(b);
  }
}

let sweepTimer=null, computeTimer=null;

function recompute(debounce){
  const go=()=>{
    const g=state.grid; if(!g)return;
    state.stats=classify(g,state.cuts);
    state.bd=breakdown(g,state.cuts);
    state.cap=capacity(state.stats,AQUIFERS.find(x=>x.id===state.aquifer),TRIALS,SEED);
    drawMap();
    renderStats();
    renderReadout();
    renderCapInputs();
    drawBreakdown();
    drawHistogram();
    clearTimeout(sweepTimer);
    sweepTimer=setTimeout(()=>{state.sweepData=computeSweep();drawSweep();},220);
  };
  clearTimeout(computeTimer);
  if(debounce) computeTimer=setTimeout(go,40); else go();
}

function selectAquifer(id){
  state.aquifer=id;
  renderTabs();
  const a=AQUIFERS.find(x=>x.id===id);
  document.getElementById('mapTitle').textContent=`Group ${id} — optimal injection zones`;
  document.getElementById('mapNote').innerHTML=
    `${a.age}. Zones are recomputed live from the property grids at 400 m. `+
    (a.pubArea? `The paper reports ${fmtInt(a.pubArea)} km² of optimal zone here.`
              : `The paper reports no optimal zones here — relax the CO₂ density cut-off to see why.`);
  const holder=mapCv.parentElement;
  holder.querySelectorAll('.loading').forEach(n=>n.remove());
  const load=document.createElement('div');
  load.className='loading';load.textContent='Loading grid…';
  load.style.cssText='position:absolute;inset:0;background:var(--map-bg)';
  holder.appendChild(load);

  loadGrid(id).then(g=>{
    load.remove();
    if(!g){ readoutEl.textContent='Could not load the grid for Group '+id+'.'; return; }
    state.grid=g; g.base=null; g.baseKey=null;
    resetCam();
    recompute();
  });
}

/* ========================= 13. WIRING ============================ */

document.getElementById('resetCuts').addEventListener('click',()=>{
  state.cuts=Object.assign({},PUBLISHED);
  renderControls(); recompute();
});
document.querySelectorAll('#sweepPick [data-sweep]').forEach(b=>{
  b.addEventListener('click',()=>{
    state.sweep=b.dataset.sweep;
    document.querySelectorAll('#sweepPick [data-sweep]').forEach(x=>
      x.setAttribute('aria-pressed',String(x===b)));
    document.getElementById('sweepNote').textContent=
      state.sweep==='por'?'Sweeping the optimal porosity cut-off with every other setting held where you left it.'
      :state.sweep==='den'?'Sweeping the optimal CO₂ density cut-off with every other setting held where you left it.'
      :'Sweeping the optimal fault setback with every other setting held where you left it.';
    state.sweepData=computeSweep(); drawSweep();
  });
});
document.getElementById('zIn').addEventListener('click',()=>{
  const r=mapCv.getBoundingClientRect();
  zoomAbout(state.cam.cx,state.cam.cy,1.5,r.width/2,r.height/2,r.width,r.height);
});
document.getElementById('zOut').addEventListener('click',()=>{
  const r=mapCv.getBoundingClientRect();
  zoomAbout(state.cam.cx,state.cam.cy,1/1.5,r.width/2,r.height/2,r.width,r.height);
});
document.getElementById('zRst').addEventListener('click',()=>{resetCam();drawMap();});

function applyTheme(t){
  if(t)document.documentElement.setAttribute('data-theme',t);
  else document.documentElement.removeAttribute('data-theme');
  const attr=document.documentElement.getAttribute('data-theme');
  const dark=attr==='dark'||(!attr&&window.matchMedia('(prefers-color-scheme: dark)').matches);
  document.getElementById('themeBtn').textContent=dark?'Light mode':'Dark mode';
  for(const k in grids){grids[k].base=null;grids[k].baseKey=null;}
  drawMap();drawBreakdown();drawHistogram();drawSweep();
}
document.getElementById('themeBtn').addEventListener('click',()=>{
  const attr=document.documentElement.getAttribute('data-theme');
  const dark=attr==='dark'||(!attr&&window.matchMedia('(prefers-color-scheme: dark)').matches);
  const next=dark?'light':'dark';
  try{localStorage.setItem('mb-theme',next);}catch(_){}
  applyTheme(next);
});

/* static prose */
document.getElementById('pubTable').innerHTML=
  `<thead><tr><th>Criterion</th><th>Optimal</th><th>Sub-optimal</th></tr></thead><tbody>
   <tr><td>Porosity</td><td>&ge; 10 % (&asymp; 100 mD)</td><td>&ge; 6 % (&asymp; 10 mD)</td></tr>
   <tr><td>CO&#8322; density</td><td>&ge; 300 kg/m&#179;</td><td>&ge; 100 kg/m&#179;</td></tr>
   <tr><td>CO&#8322; phase</td><td>supercritical</td><td>supercritical</td></tr>
   <tr><td>Distance to mapped fault</td><td>&ge; 2 km</td><td>&ge; 100 m</td></tr>
   <tr><td>Overpressure</td><td>excluded</td><td>excluded*</td></tr>
   </tbody>`;

document.getElementById('reproNote').innerHTML=
  `<p>Set to the published cut-offs, this tool's classification agrees with the archived optimal-zone
   rasters cell for cell as follows:</p>
   <p>` + AQUIFERS.map(a=>`Group ${a.id} <b>${a.repro.toFixed(1)} %</b>`).join(' &middot; ') + `</p>
   <p>The residual few per cent sits almost entirely on zone edges, where resampling the 200 m archive
   to 400 m and quantising porosity and density to 0.2 % and 4 kg/m&#179; move a boundary cell one way or
   the other. Areas land within a few per cent of the archived rasters.</p>`;

document.getElementById('quirkNote').innerHTML=
  `<p><b>1. Overpressure is not excluded from the sub-optimal class.</b> The summary table marks
   overpressure &ldquo;excluded&rdquo; for both classes, but
   <code>make_optimal_zone_map.py</code> tests it only on the optimal branch, and the archived rasters
   agree with the code: applying it to both drops the match against Group I from 97 % to 78 %. The
   checkbox above is therefore unticked for sub-optimal by default, reproducing what was actually run.</p>
   <p><b>2. The fault rasters are georeferenced to a different extent than they were used at.</b>
   Every fault grid is archived at the full 2000&times;2050 basin extent, but for Groups D, E, F and H
   the other property grids are clipped to a smaller window with a different origin. The published zone
   maps match a fault lookup done by <i>array index</i> rather than by map position &mdash; agreement is
   94&ndash;96 % that way against 86&ndash;92 % when the faults are placed by their own georeferencing.
   In other words the fault control effectively sits about 3.5 km east and 5.3 km north of where its
   header says, for those four aquifers only. This tool follows the published behaviour so that it
   reproduces the paper; Groups B, I, J and K are unaffected because all their grids share one extent.</p>`;

/* boot */
(function init(){
  let saved=null;
  try{saved=localStorage.getItem('mb-theme');}catch(_){}
  if(saved)document.documentElement.setAttribute('data-theme',saved);
  const attr=document.documentElement.getAttribute('data-theme');
  const dark=attr==='dark'||(!attr&&window.matchMedia('(prefers-color-scheme: dark)').matches);
  document.getElementById('themeBtn').textContent=dark?'Light mode':'Dark mode';

  bindTip(document.getElementById('breakdown'),document.getElementById('bdTip'),(x,y,w,h)=>bdHit(x,y,w,h));
  bindTip(document.getElementById('hist'),document.getElementById('histTip'),(x,y,w,h)=>histHit(x,y,w,h));
  bindTip(document.getElementById('sweepCap'),document.getElementById('swTip'),(x,y,w,h)=>swHit(x,y,w,h));

  document.getElementById('sweepNote').textContent=
    'Sweeping the optimal porosity cut-off with every other setting held where you left it.';

  renderControls();
  selectAquifer('I');

  let rt=null;
  window.addEventListener('resize',()=>{
    clearTimeout(rt);
    rt=setTimeout(()=>{drawMap();drawBreakdown();drawHistogram();drawSweep();},120);
  });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change',()=>{
    if(!document.documentElement.getAttribute('data-theme'))applyTheme(null);
  });
})();
