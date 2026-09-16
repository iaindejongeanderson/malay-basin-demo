/* ------------------------------------------------------------------ *
 * Malay Basin CO2 Storage Capacity Calculator
 * After de Jonge-Anderson et al. (2025), IJGGC 143, 104347.
 * ------------------------------------------------------------------ */
'use strict';

/* ============================ 1. DATA ============================= */

/* Table 2 of the paper. Area / porosity / CO2 density are measured within the
   optimal zones; thickness and NTG are fixed per group from well petrophysics.
   `pub` holds the published capacity in Gt as [P10, P50, P90] (exceedance
   convention: P10 is the HIGH estimate).
   `sub*` are statistics inside the sub-optimal zones, computed here from the
   archived grids; the paper maps those zones but does not include them in the
   volumetrics. */
const AQUIFERS = [
  { id:'B', age:'Pliocene',       area:0,     h:162, hSd:12,  ntg:0.17, ntgSd:0.09,
    por:29.8, porSd:1.0, rho:179, rhoSd:11, pub:null,
    subArea:2238,  subPor:29.8, subPorSd:1.0, subRho:179, subRhoSd:11, noOptimal:true },
  { id:'D', age:'Late Miocene',   area:3348,  h:287, hSd:262, ntg:0.16, ntgSd:0.09,
    por:22, porSd:1, rho:313, rhoSd:11, pub:[0.52,0.14,0.02],
    subArea:28030, subPor:26.0, subPorSd:2.2, subRho:246, subRhoSd:35 },
  { id:'E', age:'Late Miocene',   area:13894, h:354, hSd:280, ntg:0.27, ntgSd:0.14,
    por:21, porSd:2, rho:336, rhoSd:30, pub:[3.99,1.14,0.18],
    subArea:30610, subPor:24.4, subPorSd:3.9, subRho:347, subRhoSd:178 },
  { id:'F', age:'Mid Miocene',    area:18108, h:449, hSd:415, ntg:0.12, ntgSd:0.06,
    por:22, porSd:2, rho:340, rhoSd:41, pub:[2.90,0.84,0.13],
    subArea:41343, subPor:21.4, subPorSd:5.6, subRho:511, subRhoSd:210 },
  { id:'H', age:'Mid Miocene',    area:22290, h:393, hSd:294, ntg:0.12, ntgSd:0.10,
    por:23, porSd:2, rho:327, rhoSd:29, pub:[3.93,1.04,0.14],
    subArea:49319, subPor:21.1, subPorSd:6.2, subRho:504, subRhoSd:206 },
  { id:'I', age:'Early Miocene',  area:24924, h:610, hSd:264, ntg:0.13, ntgSd:0.08,
    por:20, porSd:3, rho:362, rhoSd:37, pub:[5.22,1.67,0.33],
    subArea:38271, subPor:19.4, subPorSd:7.9, subRho:451, subRhoSd:214 },
  { id:'J', age:'Early Miocene',  area:12898, h:272, hSd:118, ntg:0.42, ntgSd:0.17,
    por:18, porSd:3, rho:394, rhoSd:40, pub:[3.67,1.22,0.27],
    subArea:29108, subPor:13.2, subPorSd:5.8, subRho:612, subRhoSd:158 },
  { id:'K', age:'Early Miocene',  area:10643, h:383, hSd:176, ntg:0.44, ntgSd:0.13,
    por:17, porSd:4, rho:403, rhoSd:42, pub:[4.28,1.52,0.37],
    subArea:26549, subPor:12.1, subPorSd:5.0, subRho:645, subRhoSd:132 }
];

/* Global terms. The paper states the MEANS only (E = 2 %, Swirr = 27 %) and does
   not report their standard deviations. E_SD = 0.01 was chosen here because it
   reproduces the published P10-P50-P90 of Table 2 most closely across all seven
   aquifers; see the "Assumptions" note in the page. */
const GLOBAL_DEFAULT = { swirr:0.27, swirrSd:0.05, eff:0.02, effSd:0.01 };
const PUBLISHED_TOTAL = { p10:24.5, p50:7.6, p90:1.44 };

/* Whole-basin raster extent (union of all grids), EPSG:24548 metres. */
const WORLD = { x0:180000, x1:580000, y0:410000, y1:820000 };

const CLASS = { NODATA:0, NONVIABLE:1, SUBOPTIMAL:2, OPTIMAL:3 };

/* ========================= 2. SMALL MATHS ========================= */

function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

/* Standard normal CDF via Abramowitz & Stegun 7.1.26 erf. */
function normCdf(z) {
  const s = z < 0 ? -1 : 1, x = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * x);
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t
              - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return 0.5 * (1 + s * y);
}

/* Inverse standard normal CDF (Acklam). */
function normInv(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00];
  const pl = 0.02425;
  let _q, _r;
  if (p < pl) {
    _q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*_q+c[1])*_q+c[2])*_q+c[3])*_q+c[4])*_q+c[5]) /
           ((((d[0]*_q+d[1])*_q+d[2])*_q+d[3])*_q+1);
  }
  if (p > 1 - pl) {
    _q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*_q+c[1])*_q+c[2])*_q+c[3])*_q+c[4])*_q+c[5]) /
            ((((d[0]*_q+d[1])*_q+d[2])*_q+d[3])*_q+1);
  }
  _q = p - 0.5; _r = _q * _q;
  return (((((a[0]*_r+a[1])*_r+a[2])*_r+a[3])*_r+a[4])*_r+a[5])*_q /
         (((((b[0]*_r+b[1])*_r+b[2])*_r+b[3])*_r+b[4])*_r+1);
}

/* Quantile of N(mu,sd) truncated to [lo,hi]. */
function truncQuantile(p, mu, sd, lo, hi) {
  if (!(sd > 0)) return Math.min(hi, Math.max(lo, mu));
  const A = normCdf((lo - mu) / sd), B = normCdf((hi - mu) / sd);
  if (B - A < 1e-12) return Math.min(hi, Math.max(lo, mu));
  return Math.min(hi, Math.max(lo, mu + sd * normInv(A + p * (B - A))));
}

/* One draw from N(mu,sd) truncated to [lo,hi], by rejection. */
function truncNormal(rng, mu, sd, lo, hi) {
  if (!(sd > 0)) return Math.min(hi, Math.max(lo, mu));
  for (let i = 0; i < 500; i++) {
    const u1 = 1 - rng(), u2 = rng();
    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    const v = mu + sd * z;
    if (v >= lo && v <= hi) return v;
  }
  return Math.min(hi, Math.max(lo, mu));
}

/* Percentile of a pre-sorted array, linear interpolation. */
function pctOf(sorted, p) {
  if (!sorted.length) return 0;
  const i = p * (sorted.length - 1), lo = Math.floor(i), hi = Math.min(lo + 1, sorted.length - 1);
  return sorted[lo] + (i - lo) * (sorted[hi] - sorted[lo]);
}

const fmtGt = v => v === 0 ? '0'
  : v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2)
  : v >= 0.01 ? v.toFixed(3) : v.toExponential(1);
const fmtInt = v => Math.round(v).toLocaleString('en-GB');
const fmtKm2 = v => v >= 100 ? fmtInt(v) : v.toFixed(1);

/* ============ 3. PROJECTION: Kertau 1968 / UTM zone 48N ============ */
/* Everest 1948 ellipsoid, taken from the GeoTIFF GeoDoubleParams. */
const PROJ = (function () {
  const a = 6377304.063, invF = 300.8017, f = 1 / invF;
  const e2 = 2 * f - f * f, k0 = 0.9996, lon0 = 105 * Math.PI / 180, FE = 500000, FN = 0;
  const ep2 = e2 / (1 - e2);
  const e1 = (1 - Math.sqrt(1 - e2)) / (1 + Math.sqrt(1 - e2));

  function M(phi) {
    return a * ((1 - e2/4 - 3*e2*e2/64 - 5*e2*e2*e2/256) * phi
      - (3*e2/8 + 3*e2*e2/32 + 45*e2*e2*e2/1024) * Math.sin(2*phi)
      + (15*e2*e2/256 + 45*e2*e2*e2/1024) * Math.sin(4*phi)
      - (35*e2*e2*e2/3072) * Math.sin(6*phi));
  }

  function forward(lonDeg, latDeg) {
    const phi = latDeg * Math.PI / 180, lam = lonDeg * Math.PI / 180;
    const N = a / Math.sqrt(1 - e2 * Math.sin(phi) ** 2);
    const T = Math.tan(phi) ** 2, C = ep2 * Math.cos(phi) ** 2;
    const A1 = (lam - lon0) * Math.cos(phi);
    const E = FE + k0 * N * (A1 + (1 - T + C) * A1**3 / 6
            + (5 - 18*T + T*T + 72*C - 58*ep2) * A1**5 / 120);
    const Nn = FN + k0 * (M(phi) + N * Math.tan(phi) * (A1*A1/2 + (5 - T + 9*C + 4*C*C) * A1**4 / 24
            + (61 - 58*T + T*T + 600*C - 330*ep2) * A1**6 / 720));
    return [E, Nn];
  }

  function inverse(E, N) {
    const x = E - FE, y = N - FN;
    const mu = (y / k0) / (a * (1 - e2/4 - 3*e2*e2/64 - 5*e2*e2*e2/256));
    const phi1 = mu + (3*e1/2 - 27*e1**3/32) * Math.sin(2*mu)
               + (21*e1*e1/16 - 55*e1**4/32) * Math.sin(4*mu)
               + (151*e1**3/96) * Math.sin(6*mu) + (1097*e1**4/512) * Math.sin(8*mu);
    const C1 = ep2 * Math.cos(phi1) ** 2, T1 = Math.tan(phi1) ** 2;
    const N1 = a / Math.sqrt(1 - e2 * Math.sin(phi1) ** 2);
    const R1 = a * (1 - e2) / Math.pow(1 - e2 * Math.sin(phi1) ** 2, 1.5);
    const D = x / (N1 * k0);
    const phi = phi1 - (N1 * Math.tan(phi1) / R1) * (D*D/2
              - (5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*ep2) * D**4 / 24
              + (61 + 90*T1 + 298*C1 + 45*T1*T1 - 252*ep2 - 3*C1*C1) * D**6 / 720);
    const lam = lon0 + (D - (1 + 2*T1 + C1) * D**3 / 6
              + (5 - 2*C1 + 28*T1 - 3*C1*C1 + 8*ep2 + 24*T1*T1) * D**5 / 120) / Math.cos(phi1);
    return [lam * 180 / Math.PI, phi * 180 / Math.PI];
  }
  return { forward, inverse };
})();

/* ===================== 4. GRID DECODE / CACHE ===================== */

const gridCache = {};   /* id -> { w,h,codes:Uint8Array, base:Canvas, meta } */

function loadGrid(id) {
  if (gridCache[id]) return Promise.resolve(gridCache[id]);
  const meta = window.MALAY_GRIDS[id];
  if (!meta) return Promise.resolve(null);
  return new Promise(resolve => {
    const img = new Image();
    img.onload = () => {
      const w = meta.width, h = meta.height;
      const off = document.createElement('canvas');
      off.width = w; off.height = h;
      const octx = off.getContext('2d', { willReadFrequently: true });
      octx.drawImage(img, 0, 0);
      const px = octx.getImageData(0, 0, w, h).data;
      const codes = new Uint8Array(w * h);
      for (let i = 0, n = w * h; i < n; i++) codes[i] = px[i * 4];   /* R channel */
      const g = { id, w, h, codes, meta, base:null, baseTheme:null };
      gridCache[id] = g;
      resolve(g);
    };
    img.onerror = () => resolve(null);
    img.src = meta.png;
  });
}

function zoneOf(code) {
  if (code === 0) return CLASS.NODATA;
  if (code === 1) return CLASS.NONVIABLE;
  if (code === 2) return CLASS.SUBOPTIMAL;
  return CLASS.OPTIMAL;                    /* 3 = optimal, 4+ = optimal + cluster id */
}
const clusterOf = code => code >= 4 ? code - 3 : 0;

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
function hexToRgb(h) {
  h = h.replace('#', '');
  if (h.length === 3) h = h.split('').map(c => c + c).join('');
  return [parseInt(h.slice(0,2),16), parseInt(h.slice(2,4),16), parseInt(h.slice(4,6),16)];
}

/* Paint the class raster into an offscreen canvas, once per aquifer per theme. */
function buildBase(g) {
  const theme = document.documentElement.getAttribute('data-theme') || 'auto';
  const key = theme + '|' + (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'd' : 'l');
  if (g.base && g.baseTheme === key) return g.base;

  const col = {
    [CLASS.NONVIABLE]:  hexToRgb(cssVar('--zone-nonviable')),
    [CLASS.SUBOPTIMAL]: hexToRgb(cssVar('--zone-suboptimal')),
    [CLASS.OPTIMAL]:    hexToRgb(cssVar('--zone-optimal'))
  };
  const cv = document.createElement('canvas');
  cv.width = g.w; cv.height = g.h;
  const ctx = cv.getContext('2d');
  const img = ctx.createImageData(g.w, g.h);
  const d = img.data;
  for (let i = 0, n = g.w * g.h; i < n; i++) {
    const z = zoneOf(g.codes[i]);
    if (z === CLASS.NODATA) { d[i*4+3] = 0; continue; }
    const c = col[z];
    d[i*4] = c[0]; d[i*4+1] = c[1]; d[i*4+2] = c[2]; d[i*4+3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  g.base = cv; g.baseTheme = key;
  return cv;
}

/* Overlay canvas containing only the currently selected pixels. */
function buildSelLayer(g, test) {
  const cv = document.createElement('canvas');
  cv.width = g.w; cv.height = g.h;
  const ctx = cv.getContext('2d');
  const img = ctx.createImageData(g.w, g.h);
  const d = img.data;
  const c = hexToRgb(cssVar('--zone-optimal'));
  let count = 0;
  for (let i = 0, n = g.w * g.h; i < n; i++) {
    if (!test(g.codes[i], i)) { d[i*4+3] = 0; continue; }
    count++;
    d[i*4] = c[0]; d[i*4+1] = c[1]; d[i*4+2] = c[2]; d[i*4+3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  return { canvas: cv, count };
}

/* ========================== 5. APP STATE ========================== */

const state = {
  aquifer: 'I',
  mode: 'all',                 /* all | clus | rect */
  cluster: null,
  rect: null,                  /* {x0,y0,x1,y1} in world metres */
  selArea: null,               /* km2 from map, null = use params.area */
  params: null,
  global: Object.assign({}, GLOBAL_DEFAULT),
  trials: 10000,
  seed: 42,
  result: null,
  basin: null,
  cam: null,
  grid: null,
  selLayer: null
};

function defaultParams(id) {
  const a = AQUIFERS.find(x => x.id === id);
  return {
    area: a.area, h: a.h, hSd: a.hSd, ntg: a.ntg, ntgSd: a.ntgSd,
    por: a.por, porSd: a.porSd, rho: a.rho, rhoSd: a.rhoSd,
    swirr: state.global.swirr, swirrSd: state.global.swirrSd,
    eff: state.global.eff, effSd: state.global.effSd
  };
}

/* ======================= 6. MONTE CARLO =========================== */

/* Returns { samples(sorted), p10, p50, p90, mean, det } in Gt. */
function runMC(p, trials, seed) {
  const rng = mulberry32(seed >>> 0);
  const out = new Float64Array(trials);
  const A = p.area * 1e6;                      /* km2 -> m2 */
  for (let i = 0; i < trials; i++) {
    const h   = truncNormal(rng, p.h,   p.hSd,   0, Infinity);
    const ntg = truncNormal(rng, p.ntg, p.ntgSd, 0, 1);
    const por = truncNormal(rng, p.por, p.porSd, 0, 100) / 100;
    const sw  = truncNormal(rng, p.swirr, p.swirrSd, 0, 1);
    const e   = truncNormal(rng, p.eff, p.effSd, 0, 1);
    const rho = truncNormal(rng, p.rho, p.rhoSd, 0, Infinity);
    out[i] = A * h * ntg * por * (1 - sw) * e * rho / 1e12;   /* kg -> Gt */
  }
  const sorted = Float64Array.from(out).sort();
  let s = 0; for (let i = 0; i < trials; i++) s += out[i];
  return {
    samples: out, sorted,
    p10: pctOf(sorted, 0.90),      /* exceedance: P10 = high */
    p50: pctOf(sorted, 0.50),
    p90: pctOf(sorted, 0.10),      /* exceedance: P90 = low  */
    mean: s / trials,
    det: A * p.h * p.ntg * (p.por/100) * (1 - p.swirr) * p.eff * p.rho / 1e12
  };
}

function effectiveParams() {
  const p = Object.assign({}, state.params);
  if (state.selArea !== null) p.area = state.selArea;
  p.swirr = state.global.swirr; p.swirrSd = state.global.swirrSd;
  p.eff = state.global.eff;     p.effSd = state.global.effSd;
  return p;
}

/* ===================== 7. CANVAS PLUMBING ========================= */

function fitCanvas(cv, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const w = cv.clientWidth || cv.parentElement.clientWidth;
  const h = cssH || cv.clientHeight;
  cv.width = Math.max(1, Math.round(w * dpr));
  cv.height = Math.max(1, Math.round(h * dpr));
  cv.style.height = h + 'px';
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function roundTopRect(ctx, x, y, w, h, r) {
  r = Math.max(0, Math.min(r, w / 2, h));
  ctx.beginPath();
  ctx.moveTo(x, y + h);
  ctx.lineTo(x, y + r);
  ctx.quadraticCurveTo(x, y, x + r, y);
  ctx.lineTo(x + w - r, y);
  ctx.quadraticCurveTo(x + w, y, x + w, y + r);
  ctx.lineTo(x + w, y + h);
  ctx.closePath();
  ctx.fill();
}

/* ========================== 8. THE MAP ============================ */

const mapCv = document.getElementById('map');
const readoutEl = document.getElementById('readout');

function resetCam() {
  const holder = mapCv.parentElement;
  const w = holder.clientWidth, h = holder.clientHeight;
  const k = Math.min(w / (WORLD.x1 - WORLD.x0), h / (WORLD.y1 - WORLD.y0)) * 0.94;
  state.cam = { cx: (WORLD.x0 + WORLD.x1) / 2, cy: (WORLD.y0 + WORLD.y1) / 2, k };
}

function w2s(x, y, w, h) {
  const c = state.cam;
  return [w / 2 + (x - c.cx) * c.k, h / 2 - (y - c.cy) * c.k];
}
function s2w(sx, sy, w, h) {
  const c = state.cam;
  return [c.cx + (sx - w / 2) / c.k, c.cy - (sy - h / 2) / c.k];
}

function drawMap() {
  const holder = mapCv.parentElement;
  const cssW = holder.clientWidth, cssH = holder.clientHeight;
  const dpr = window.devicePixelRatio || 1;
  mapCv.width = Math.round(cssW * dpr);
  mapCv.height = Math.round(cssH * dpr);
  const ctx = mapCv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const g = state.grid;
  if (!g) return;
  const m = g.meta;

  const [sx0, sy0] = w2s(m.originX, m.originY, cssW, cssH);
  const dw = m.width * m.cell * state.cam.k;
  const dh = m.height * m.cell * state.cam.k;

  ctx.imageSmoothingEnabled = state.cam.k * m.cell < 1;
  const hasSel = !!state.selLayer;
  ctx.globalAlpha = hasSel ? 0.35 : 1;
  ctx.drawImage(buildBase(g), sx0, sy0, dw, dh);
  ctx.globalAlpha = 1;
  if (hasSel) ctx.drawImage(state.selLayer.canvas, sx0, sy0, dw, dh);

  drawGraticule(ctx, cssW, cssH);

  /* rectangle selection outline */
  if (state.rect) {
    const r = state.rect;
    const [ax, ay] = w2s(Math.min(r.x0, r.x1), Math.max(r.y0, r.y1), cssW, cssH);
    const [bx, by] = w2s(Math.max(r.x0, r.x1), Math.min(r.y0, r.y1), cssW, cssH);
    ctx.save();
    ctx.strokeStyle = cssVar('--text-primary');
    ctx.lineWidth = 2;
    ctx.setLineDash([]);
    ctx.strokeRect(ax, ay, bx - ax, by - ay);
    ctx.restore();
  }
  drawScaleBar(ctx, cssW, cssH);
}

function drawGraticule(ctx, w, h) {
  const [lonA, latA] = PROJ.inverse(...s2w(0, h, w, h));
  const [lonB, latB] = PROJ.inverse(...s2w(w, 0, w, h));
  const span = Math.max(lonB - lonA, latB - latA);
  const step = span > 4 ? 1 : span > 2 ? 0.5 : span > 0.8 ? 0.25 : 0.1;

  ctx.save();
  ctx.strokeStyle = cssVar('--grid');
  ctx.fillStyle = cssVar('--muted');
  ctx.lineWidth = 1;
  ctx.font = '10px system-ui, sans-serif';

  for (let lon = Math.ceil(lonA / step) * step; lon <= lonB; lon += step) {
    ctx.beginPath();
    for (let i = 0; i <= 24; i++) {
      const lat = latA + (latB - latA) * i / 24;
      const [sx, sy] = w2s(...PROJ.forward(lon, lat), w, h);
      i ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy);
    }
    ctx.stroke();
    /* label meridians along the top edge: the bottom carries the hint and scale bar */
    const [lx, ly] = w2s(...PROJ.forward(lon, latB), w, h);
    ctx.fillText(lon.toFixed(step < 1 ? 2 : 0) + '°E', lx + 3, Math.max(ly + 11, 11));
  }
  for (let lat = Math.ceil(latA / step) * step; lat <= latB; lat += step) {
    ctx.beginPath();
    for (let i = 0; i <= 24; i++) {
      const lon = lonA + (lonB - lonA) * i / 24;
      const [sx, sy] = w2s(...PROJ.forward(lon, lat), w, h);
      i ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy);
    }
    ctx.stroke();
    const [lx, ly] = w2s(...PROJ.forward(lonA, lat), w, h);
    ctx.fillText(lat.toFixed(step < 1 ? 2 : 0) + '°N', Math.max(lx + 3, 4), ly - 3);
  }
  ctx.restore();
}

function drawScaleBar(ctx, w, h) {
  const targetPx = Math.min(140, w * 0.3);
  const metres = targetPx / state.cam.k;
  const nice = [1e3,2e3,5e3,1e4,2e4,5e4,1e5,2e5];
  let pick = nice[0];
  for (const n of nice) if (n <= metres) pick = n;
  const px = pick * state.cam.k;
  const x = w - px - 14, y = h - 16;
  ctx.save();
  ctx.strokeStyle = cssVar('--text-primary');
  ctx.fillStyle = cssVar('--text-secondary');
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(x, y - 4); ctx.lineTo(x, y); ctx.lineTo(x + px, y); ctx.lineTo(x + px, y - 4);
  ctx.stroke();
  ctx.font = '11px system-ui, sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText((pick / 1000) + ' km', x + px / 2, y - 7);
  ctx.restore();
}

/* --- map interaction --- */
let drag = null, rectDrag = null;

mapCv.addEventListener('pointerdown', e => {
  mapCv.setPointerCapture(e.pointerId);
  const r = mapCv.getBoundingClientRect();
  const sx = e.clientX - r.left, sy = e.clientY - r.top;
  if (state.mode === 'rect') {
    const [wx, wy] = s2w(sx, sy, r.width, r.height);
    rectDrag = { x0: wx, y0: wy };
  } else {
    drag = { sx, sy, cx: state.cam.cx, cy: state.cam.cy, moved: false };
  }
});

mapCv.addEventListener('pointermove', e => {
  const r = mapCv.getBoundingClientRect();
  const sx = e.clientX - r.left, sy = e.clientY - r.top;

  if (rectDrag) {
    const [wx, wy] = s2w(sx, sy, r.width, r.height);
    state.rect = { x0: rectDrag.x0, y0: rectDrag.y0, x1: wx, y1: wy };
    drawMap();
    return;
  }
  if (drag) {
    if (Math.abs(sx - drag.sx) + Math.abs(sy - drag.sy) > 3) drag.moved = true;
    state.cam.cx = drag.cx - (sx - drag.sx) / state.cam.k;
    state.cam.cy = drag.cy + (sy - drag.sy) / state.cam.k;
    drawMap();
    return;
  }
  updateReadout(sx, sy, r.width, r.height);
});

mapCv.addEventListener('pointerup', e => {
  const r = mapCv.getBoundingClientRect();
  if (rectDrag) {
    const [wx, wy] = s2w(e.clientX - r.left, e.clientY - r.top, r.width, r.height);
    if (Math.abs(wx - rectDrag.x0) > 500 && Math.abs(wy - rectDrag.y0) > 500) {
      state.rect = { x0: rectDrag.x0, y0: rectDrag.y0, x1: wx, y1: wy };
      applyRectSelection();
    } else {
      state.rect = null; clearSelection();
    }
    rectDrag = null;
    return;
  }
  if (drag && !drag.moved && state.mode === 'clus') {
    pickCluster(e.clientX - r.left, e.clientY - r.top, r.width, r.height);
  }
  drag = null;
});

mapCv.addEventListener('pointerleave', () => { readoutEl.innerHTML = readoutBase(); });

mapCv.addEventListener('wheel', e => {
  e.preventDefault();
  const r = mapCv.getBoundingClientRect();
  const sx = e.clientX - r.left, sy = e.clientY - r.top;
  const [wx, wy] = s2w(sx, sy, r.width, r.height);
  const f = Math.exp(-e.deltaY * 0.0015);
  zoomAbout(wx, wy, f, sx, sy, r.width, r.height);
}, { passive: false });

function zoomAbout(wx, wy, f, sx, sy, w, h) {
  const c = state.cam;
  const base = Math.min(w / (WORLD.x1 - WORLD.x0), h / (WORLD.y1 - WORLD.y0)) * 0.94;
  c.k = Math.max(base * 0.8, Math.min(base * 60, c.k * f));
  c.cx = wx - (sx - w / 2) / c.k;
  c.cy = wy + (sy - h / 2) / c.k;
  drawMap();
}

function gridIndexAt(wx, wy) {
  const g = state.grid; if (!g) return -1;
  const m = g.meta;
  const px = Math.floor((wx - m.originX) / m.cell);
  const py = Math.floor((m.originY - wy) / m.cell);
  if (px < 0 || py < 0 || px >= m.width || py >= m.height) return -1;
  return py * m.width + px;
}

const ZONE_NAME = { 0:'no data', 1:'non-viable', 2:'sub-optimal', 3:'optimal' };

function readoutBase() {
  const g = state.grid;
  const a = AQUIFERS.find(x => x.id === state.aquifer);
  const gridArea = g ? g.meta.optimalAreaKm2 : 0;
  let sel;
  if (state.selArea !== null) {
    sel = `<b>${fmtKm2(state.selArea)} km&#178;</b> selected`;
    if (state.mode === 'clus' && state.cluster) sel += ` (cluster ${state.cluster})`;
    if (state.mode === 'rect') sel += ' (in box)';
  } else {
    sel = `<b>${fmtInt(a.area)} km&#178;</b> &mdash; Table 2 optimal-zone area`;
  }
  const alt = (state.selArea === null && gridArea)
    ? ` &middot; grid measures ${fmtInt(gridArea)} km&#178; <a href="#" id="useGrid">use this</a>`
    : '';
  return `<span>${sel}${alt}</span>`;
}

function updateReadout(sx, sy, w, h) {
  const [wx, wy] = s2w(sx, sy, w, h);
  const [lon, lat] = PROJ.inverse(wx, wy);
  const i = gridIndexAt(wx, wy);
  const code = i >= 0 ? state.grid.codes[i] : 0;
  const z = zoneOf(code);
  const cl = clusterOf(code);
  let s = readoutBase();
  s += `<span>${lat.toFixed(3)}°N ${lon.toFixed(3)}°E</span>`;
  s += `<span>${fmtInt(wx)} E, ${fmtInt(wy)} N</span>`;
  s += `<span><b>${ZONE_NAME[z]}</b>${cl ? ' &middot; cluster ' + cl : ''}</span>`;
  readoutEl.innerHTML = s;
}

function pickCluster(sx, sy, w, h) {
  const [wx, wy] = s2w(sx, sy, w, h);
  const i = gridIndexAt(wx, wy);
  if (i < 0) return;
  const code = state.grid.codes[i];
  if (zoneOf(code) !== CLASS.OPTIMAL) { clearSelection(); return; }
  const cid = clusterOf(code);
  const g = state.grid;
  if (cid === 0) {
    /* optimal pixel with no cluster label: select all such pixels */
    state.selLayer = buildSelLayer(g, c => c === 3);
    state.cluster = null;
  } else {
    state.selLayer = buildSelLayer(g, c => c === code);
    state.cluster = cid;
  }
  const cellKm2 = (g.meta.cell * g.meta.cell) / 1e6;
  state.selArea = state.selLayer.count * cellKm2;
  state.rect = null;
  afterSelection();
}

function applyRectSelection() {
  const g = state.grid, m = g.meta, r = state.rect;
  const x0 = Math.min(r.x0, r.x1), x1 = Math.max(r.x0, r.x1);
  const y0 = Math.min(r.y0, r.y1), y1 = Math.max(r.y0, r.y1);
  const px0 = Math.max(0, Math.floor((x0 - m.originX) / m.cell));
  const px1 = Math.min(m.width  - 1, Math.floor((x1 - m.originX) / m.cell));
  const py0 = Math.max(0, Math.floor((m.originY - y1) / m.cell));
  const py1 = Math.min(m.height - 1, Math.floor((m.originY - y0) / m.cell));
  const inBox = (code, i) => {
    if (zoneOf(code) !== CLASS.OPTIMAL) return false;
    const x = i % m.width, y = (i / m.width) | 0;
    return x >= px0 && x <= px1 && y >= py0 && y <= py1;
  };
  state.selLayer = buildSelLayer(g, inBox);
  state.selArea = state.selLayer.count * (m.cell * m.cell) / 1e6;
  state.cluster = null;
  afterSelection();
}

function clearSelection() {
  state.selArea = null; state.selLayer = null; state.cluster = null; state.rect = null;
  afterSelection();
}

function afterSelection() {
  drawMap();
  readoutEl.innerHTML = readoutBase();
  renderParams();
  recompute();
}

/* ======================= 9. PARAM TABLE =========================== */

const PARAM_ROWS = [
  { key:'area',  sd:null,     label:'Area, A',            unit:'km²',  src:'Table 2 / map', step:1,    dp:0 },
  { key:'h',     sd:'hSd',    label:'Thickness, h',       unit:'m',         src:'well petrophysics', step:1, dp:0 },
  { key:'ntg',   sd:'ntgSd',  label:'Net-to-gross, NTG',  unit:'frac.',     src:'well petrophysics', step:0.01, dp:2 },
  { key:'por',   sd:'porSd',  label:'Porosity, φ',   unit:'%',         src:'porosity–depth model', step:0.5, dp:1 },
  { key:'swirr', sd:'swirrSd',label:'Irreducible water, S<sub>wirr</sub>', unit:'frac.', src:'literature (0.27)', step:0.01, dp:2, glob:true },
  { key:'eff',   sd:'effSd',  label:'Efficiency, E',      unit:'frac.',     src:'literature (2 %)', step:0.005, dp:3, glob:true },
  { key:'rho',   sd:'rhoSd',  label:'CO₂ density, ρ', unit:'kg/m³', src:'CoolProp EOS', step:1, dp:0 }
];

function renderParams() {
  const tb = document.getElementById('paramRows');
  const p = effectiveParams();
  const dflt = defaultParams(state.aquifer);
  tb.innerHTML = '';
  for (const r of PARAM_ROWS) {
    const tr = document.createElement('tr');
    const isArea = r.key === 'area';
    const locked = isArea && state.selArea !== null;
    /* Global terms are "modified" relative to the published defaults; per-aquifer
       terms relative to that aquifer's own published row. */
    const refM = r.glob ? GLOBAL_DEFAULT[r.key] : dflt[r.key];
    const refS = r.sd ? (r.glob ? GLOBAL_DEFAULT[r.sd] : dflt[r.sd]) : null;
    const dirtyM = !locked && Math.abs(p[r.key] - refM) > 1e-9;
    const dirtyS = r.sd && Math.abs(p[r.sd] - refS) > 1e-9;
    tr.innerHTML =
      `<td><div class="pname">${r.label}</div><div class="psrc">${locked ? 'from map selection' : r.src}</div></td>` +
      `<td class="n"><input type="number" data-k="${r.key}" step="${r.step}" value="${p[r.key].toFixed(r.dp)}"` +
        `${locked ? ' disabled' : ''} class="${dirtyM ? 'dirty' : ''}"></td>` +
      `<td class="n">${r.sd
        ? `<input type="number" data-k="${r.sd}" step="${r.step}" min="0" value="${p[r.sd].toFixed(r.dp)}" class="${dirtyS ? 'dirty' : ''}">`
        : '<span class="psrc">fixed</span>'}</td>` +
      `<td class="psrc">${r.unit}</td>`;
    tb.appendChild(tr);
  }
  tb.querySelectorAll('input').forEach(inp => {
    inp.addEventListener('input', () => {
      const k = inp.dataset.k;
      const v = parseFloat(inp.value);
      if (!isFinite(v)) return;
      if (k === 'swirr' || k === 'swirrSd' || k === 'eff' || k === 'effSd') state.global[k] = v;
      else state.params[k] = v;
      recompute(true);
    });
  });
}

/* ========================== 10. CHARTS ============================ */

function bindTip(cv, tipEl, hitFn) {
  cv.addEventListener('mousemove', e => {
    const r = cv.getBoundingClientRect();
    const hit = hitFn(e.clientX - r.left, e.clientY - r.top, r.width, r.height);
    if (!hit) { tipEl.classList.remove('on'); return; }
    tipEl.innerHTML = hit.html;
    tipEl.classList.add('on');
    const tw = tipEl.offsetWidth, th = tipEl.offsetHeight;
    let x = hit.x - tw / 2, y = hit.y - th - 10;
    x = Math.max(2, Math.min(r.width - tw - 2, x));
    if (y < 2) y = hit.y + 14;
    tipEl.style.left = x + 'px';
    tipEl.style.top = y + 'px';
  });
  cv.addEventListener('mouseleave', () => tipEl.classList.remove('on'));
}

let histHit = () => null, tornHit = () => null, cmpHit = () => null;

function drawHistogram() {
  const cv = document.getElementById('hist');
  const { ctx, w, h } = fitCanvas(cv, 230);
  const res = state.result;
  if (!res || !res.p10) {
    ctx.fillStyle = cssVar('--muted');
    ctx.font = '13px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('No optimal zones — no capacity to show', w / 2, h / 2);
    histHit = () => null;
    return;
  }

  const padL = 44, padR = 12, padT = 16, padB = 30;
  const plotW = w - padL - padR, plotH = h - padT - padB;
  const hi = Math.max(res.p10 * 1.35, pctOf(res.sorted, 0.985));
  const NB = 42, bw = hi / NB;
  const bins = new Array(NB).fill(0);
  for (let i = 0; i < res.samples.length; i++) {
    const b = Math.floor(res.samples[i] / bw);
    if (b >= 0 && b < NB) bins[b]++;
  }
  const maxC = Math.max(...bins) || 1;
  const X = v => padL + (v / hi) * plotW;
  const Y = c => padT + plotH - (c / maxC) * plotH;

  /* gridlines */
  ctx.strokeStyle = cssVar('--grid'); ctx.lineWidth = 1;
  ctx.fillStyle = cssVar('--muted');
  ctx.font = '10px system-ui, sans-serif';
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (let i = 0; i <= 4; i++) {
    const c = maxC * i / 4, y = Math.round(Y(c)) + 0.5;
    ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(w - padR, y); ctx.stroke();
    ctx.fillText(Math.round(c / state.trials * 1000) / 10 + '%', padL - 7, y);
  }

  /* bars */
  ctx.fillStyle = cssVar('--series-1');
  const pxw = plotW / NB;
  for (let i = 0; i < NB; i++) {
    if (!bins[i]) continue;
    const x = padL + i * pxw, y = Y(bins[i]);
    roundTopRect(ctx, x + 1, y, Math.max(1, pxw - 2), padT + plotH - y, 3);
  }

  /* axis */
  ctx.strokeStyle = cssVar('--axis');
  ctx.beginPath();
  ctx.moveTo(padL, padT + plotH + 0.5); ctx.lineTo(w - padR, padT + plotH + 0.5);
  ctx.stroke();
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  ctx.fillStyle = cssVar('--muted');
  for (let i = 0; i <= 4; i++) {
    const v = hi * i / 4;
    ctx.fillText(fmtGt(v), X(v), padT + plotH + 6);
  }
  ctx.fillStyle = cssVar('--text-secondary');
  ctx.fillText('CO₂ storage capacity (Gt)', padL + plotW / 2, padT + plotH + 18);

  /* percentile markers */
  const marks = [
    { v: res.p90, t: 'P90' }, { v: res.p50, t: 'P50' }, { v: res.p10, t: 'P10' }
  ];
  ctx.save();
  ctx.strokeStyle = cssVar('--text-primary');
  ctx.fillStyle = cssVar('--text-primary');
  ctx.lineWidth = 1.5;
  ctx.font = '600 10px system-ui, sans-serif';
  ctx.textBaseline = 'alphabetic';
  for (const m of marks) {
    if (m.v > hi) continue;
    const x = Math.round(X(m.v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(x, padT - 4); ctx.lineTo(x, padT + plotH); ctx.stroke();
    ctx.textAlign = x > w - padR - 40 ? 'right' : 'left';
    ctx.fillText(m.t, x + (x > w - padR - 40 ? -3 : 3), padT + 4);
  }
  ctx.restore();

  histHit = (mx, my) => {
    if (mx < padL || mx > w - padR || my < padT || my > padT + plotH) return null;
    const i = Math.floor((mx - padL) / pxw);
    if (i < 0 || i >= NB || !bins[i]) return null;
    return {
      x: padL + (i + 0.5) * pxw, y: Y(bins[i]),
      html: `<b>${fmtGt(i*bw)} – ${fmtGt((i+1)*bw)} Gt</b><br>` +
            `${fmtInt(bins[i])} of ${fmtInt(state.trials)} trials (${(bins[i]/state.trials*100).toFixed(1)} %)`
    };
  };
}

function drawTornado() {
  const cv = document.getElementById('tornado');
  const { ctx, w, h } = fitCanvas(cv, 210);
  const p = effectiveParams();
  if (!p.area || !state.result || !state.result.p50) {
    ctx.fillStyle = cssVar('--muted');
    ctx.font = '13px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('No capacity to analyse', w / 2, h / 2);
    tornHit = () => null;
    return;
  }

  const base = p.area * 1e6 * p.h * p.ntg * (p.por/100) * (1 - p.swirr) * p.eff * p.rho / 1e12;
  const vars = [
    { k:'Thickness',   m:p.h,     s:p.hSd,     lo:0, hi:Infinity, f:v=>v/p.h },
    { k:'Net-to-gross',m:p.ntg,   s:p.ntgSd,   lo:0, hi:1,        f:v=>v/p.ntg },
    { k:'Porosity',    m:p.por,   s:p.porSd,   lo:0, hi:100,      f:v=>v/p.por },
    { k:'CO₂ density', m:p.rho, s:p.rhoSd, lo:0, hi:Infinity, f:v=>v/p.rho },
    { k:'Efficiency',  m:p.eff,   s:p.effSd,   lo:0, hi:1,        f:v=>v/p.eff },
    { k:'Irreducible water', m:p.swirr, s:p.swirrSd, lo:0, hi:1,  f:v=>(1-v)/(1-p.swirr) }
  ];
  for (const v of vars) {
    const q10 = truncQuantile(0.10, v.m, v.s, v.lo, v.hi);
    const q90 = truncQuantile(0.90, v.m, v.s, v.lo, v.hi);
    const a = base * v.f(q10), b = base * v.f(q90);
    v.lowV = Math.min(a, b); v.highV = Math.max(a, b);
    v.q10 = q10; v.q90 = q90;
    v.swing = v.highV - v.lowV;
  }
  vars.sort((a, b) => b.swing - a.swing);

  const padL = 108, padR = 14, padT = 8, padB = 26;
  const plotW = w - padL - padR;
  const rowH = (h - padT - padB) / vars.length;
  const lo = Math.min(base, ...vars.map(v => v.lowV));
  const hi = Math.max(base, ...vars.map(v => v.highV));
  const pad = (hi - lo) * 0.06 || 1;
  const X = v => padL + ((v - (lo - pad)) / ((hi + pad) - (lo - pad))) * plotW;

  ctx.font = '11px system-ui, sans-serif';
  const bars = [];
  for (let i = 0; i < vars.length; i++) {
    const v = vars[i];
    const y = padT + i * rowH + rowH * 0.18, bh = rowH * 0.64;
    const xb = X(base);
    const xl = X(v.lowV), xr = X(v.highV);
    ctx.fillStyle = cssVar('--diverge-lo');
    ctx.fillRect(xl, y, Math.max(0, xb - xl - 1), bh);
    ctx.fillStyle = cssVar('--diverge-hi');
    ctx.fillRect(xb + 1, y, Math.max(0, xr - xb - 1), bh);

    ctx.fillStyle = cssVar('--text-secondary');
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillText(v.k, padL - 9, y + bh / 2);
    bars.push({ v, y, bh, xl, xr });
  }

  /* base line */
  const xb = Math.round(X(base)) + 0.5;
  ctx.strokeStyle = cssVar('--text-primary');
  ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(xb, padT); ctx.lineTo(xb, h - padB + 2); ctx.stroke();

  ctx.fillStyle = cssVar('--muted');
  ctx.font = '10px system-ui, sans-serif';
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let i = 0; i <= 4; i++) {
    const v = (lo - pad) + ((hi + pad) - (lo - pad)) * i / 4;
    ctx.fillText(fmtGt(v), X(v), h - padB + 6);
  }
  ctx.fillStyle = cssVar('--text-secondary');
  ctx.fillText('Capacity at mean inputs = ' + fmtGt(base) + ' Gt', padL + plotW / 2, h - padB + 17);

  tornHit = (mx, my) => {
    for (const b of bars) {
      if (my >= b.y && my <= b.y + b.bh && mx >= Math.min(b.xl, b.xr) - 2 && mx <= Math.max(b.xl, b.xr) + 2) {
        return {
          x: mx, y: b.y,
          html: `<b>${b.v.k}</b><br>P10 input ${b.v.q10.toFixed(3)} → ${fmtGt(b.v.lowV)} Gt<br>` +
                `P90 input ${b.v.q90.toFixed(3)} → ${fmtGt(b.v.highV)} Gt<br>` +
                `swing ${fmtGt(b.v.swing)} Gt`
        };
      }
    }
    return null;
  };
}

function drawComparison() {
  const cv = document.getElementById('cmp');
  const { ctx, w, h } = fitCanvas(cv, 260);
  const rows = state.basin;
  if (!rows) return;

  const padL = 46, padR = 12, padT = 26, padB = 42;
  const plotW = w - padL - padR, plotH = h - padT - padB;
  const hi = Math.max(0.1, ...rows.map(r => Math.max(r.res.p10, r.pub ? r.pub[0] : 0))) * 1.05;
  const Y = v => padT + plotH - (v / hi) * plotH;

  ctx.strokeStyle = cssVar('--grid'); ctx.lineWidth = 1;
  ctx.fillStyle = cssVar('--muted');
  ctx.font = '10px system-ui, sans-serif';
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (let i = 0; i <= 4; i++) {
    const v = hi * i / 4, y = Math.round(Y(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(w - padR, y); ctx.stroke();
    ctx.fillText(fmtGt(v), padL - 7, y);
  }

  const gw = plotW / rows.length;
  const bw = Math.min(26, gw * 0.3);
  const bars = [];
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i];
    const cx = padL + gw * (i + 0.5);
    const x1 = cx - bw - 1, x2 = cx + 1;

    ctx.fillStyle = cssVar('--series-1');
    const y1 = Y(r.res.p50);
    roundTopRect(ctx, x1, y1, bw, padT + plotH - y1, 3);
    bars.push({ x:x1, y:y1, w:bw, h:padT+plotH-y1, r, kind:'run' });

    if (r.pub) {
      ctx.fillStyle = cssVar('--series-2');
      const y2 = Y(r.pub[1]);
      roundTopRect(ctx, x2, y2, bw, padT + plotH - y2, 3);
      bars.push({ x:x2, y:y2, w:bw, h:padT+plotH-y2, r, kind:'pub' });
    }

    /* P90-P10 range whisker on the computed bar */
    ctx.strokeStyle = cssVar('--text-primary');
    ctx.lineWidth = 1.5;
    const xm = Math.round(x1 + bw / 2) + 0.5;
    ctx.beginPath();
    ctx.moveTo(xm, Y(r.res.p90)); ctx.lineTo(xm, Y(r.res.p10));
    ctx.moveTo(xm - 3, Y(r.res.p10)); ctx.lineTo(xm + 3, Y(r.res.p10));
    ctx.moveTo(xm - 3, Y(r.res.p90)); ctx.lineTo(xm + 3, Y(r.res.p90));
    ctx.stroke();

    ctx.fillStyle = r.id === state.aquifer ? cssVar('--text-primary') : cssVar('--muted');
    ctx.font = (r.id === state.aquifer ? '600 ' : '') + '11px system-ui, sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.fillText(r.id, cx, padT + plotH + 7);
  }

  ctx.strokeStyle = cssVar('--axis'); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(padL, padT + plotH + 0.5); ctx.lineTo(w - padR, padT + plotH + 0.5); ctx.stroke();

  /* legend */
  ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
  ctx.font = '11px system-ui, sans-serif';
  let lx = padL;
  const leg = [['--series-1', 'This run (P50)'], ['--series-2', 'Published P50 (Table 2)']];
  for (const [c, t] of leg) {
    ctx.fillStyle = cssVar(c);
    ctx.fillRect(lx, padT - 17, 10, 10);
    ctx.fillStyle = cssVar('--text-secondary');
    ctx.fillText(t, lx + 15, padT - 12);
    lx += ctx.measureText(t).width + 34;
  }
  ctx.strokeStyle = cssVar('--text-primary'); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(lx + 4, padT - 18); ctx.lineTo(lx + 4, padT - 6); ctx.stroke();
  ctx.fillStyle = cssVar('--text-secondary');
  ctx.fillText('P90–P10 range', lx + 11, padT - 12);

  ctx.textAlign = 'center';
  ctx.fillStyle = cssVar('--text-secondary');
  ctx.font = '11px system-ui, sans-serif';
  ctx.fillText('Aquifer (youngest → oldest)', padL + plotW / 2, padT + plotH + 26);

  cmpHit = (mx, my) => {
    for (const b of bars) {
      if (mx >= b.x && mx <= b.x + b.w && my >= b.y - 4 && my <= b.y + b.h) {
        const r = b.r;
        return b.kind === 'run'
          ? { x:b.x+b.w/2, y:b.y, html:`<b>Group ${r.id}</b> — this run<br>P10 ${fmtGt(r.res.p10)} &middot; P50 ${fmtGt(r.res.p50)} &middot; P90 ${fmtGt(r.res.p90)} Gt` }
          : { x:b.x+b.w/2, y:b.y, html:`<b>Group ${r.id}</b> — published<br>P10 ${r.pub[0]} &middot; P50 ${r.pub[1]} &middot; P90 ${r.pub[2]} Gt` };
      }
    }
    return null;
  };
}

/* ======================= 11. RENDER / WIRE ======================== */

function renderStats() {
  const res = state.result, p = effectiveParams();
  const a = AQUIFERS.find(x => x.id === state.aquifer);
  const el = document.getElementById('stats');
  if (!p.area) {
    el.innerHTML =
      `<div class="stat" style="grid-column:1/-1">
         <div class="k">Group ${a.id}</div>
         <div class="v" style="font-size:16px">No optimal zones</div>
         <div class="u">Table 2 reports no optimal zones for this aquifer &mdash; it is too shallow for
         CO&#8322; to reach the 300 kg/m&#179; density cut-off. Enter an area above, or draw a box on the map,
         to explore a hypothetical case using its sub-optimal zone properties.</div>
       </div>`;
    return;
  }
  el.innerHTML =
    `<div class="stat hero"><div class="k">P50 &middot; median</div>
       <div class="v">${fmtGt(res.p50)}<span class="u"> Gt</span></div></div>
     <div class="stat"><div class="k">P10 &middot; high</div>
       <div class="v">${fmtGt(res.p10)}<span class="u"> Gt</span></div></div>
     <div class="stat"><div class="k">P90 &middot; low</div>
       <div class="v">${fmtGt(res.p90)}<span class="u"> Gt</span></div></div>
     <div class="stat"><div class="k">Mean</div>
       <div class="v">${fmtGt(res.mean)}<span class="u"> Gt</span></div></div>`;

  const cmp = document.getElementById('cmpPublished');
  if (a.pub && state.selArea === null) {
    const d = (res.p50 - a.pub[1]) / a.pub[1] * 100;
    const sign = d >= 0 ? '+' : '−';
    cmp.innerHTML = `Published Table 2 for Group ${a.id}: <b>P10 ${a.pub[0]} &middot; P50 ${a.pub[1]} &middot; P90 ${a.pub[2]} Gt</b>. ` +
      `This run's P50 differs by ${sign}${Math.abs(d).toFixed(0)} %.`;
  } else if (state.selArea !== null) {
    cmp.innerHTML = `Area taken from your map selection (${fmtKm2(state.selArea)} km&#178;), so this is not directly comparable with Table 2.`;
  } else {
    cmp.innerHTML = '';
  }
}

function renderBasin() {
  const tb = document.querySelector('#basinTable tbody');
  const tf = document.querySelector('#basinTable tfoot');
  const rows = state.basin;
  tb.innerHTML = '';
  for (const r of rows) {
    const tr = document.createElement('tr');
    if (r.id === state.aquifer) tr.className = 'sel';
    tr.style.cursor = 'pointer';
    tr.innerHTML =
      `<td><span class="swatch-cell"><b>${r.id}</b>
         <span class="psrc">${r.age}</span></span></td>` +
      `<td>${r.p.area ? fmtInt(r.p.area) : '—'}</td>` +
      `<td>${fmtInt(r.p.h)} &plusmn; ${fmtInt(r.p.hSd)}</td>` +
      `<td>${r.p.ntg.toFixed(2)} &plusmn; ${r.p.ntgSd.toFixed(2)}</td>` +
      `<td>${r.p.area ? r.p.por.toFixed(0) + ' &plusmn; ' + r.p.porSd.toFixed(0) : '—'}</td>` +
      `<td>${r.p.area ? fmtInt(r.p.rho) + ' &plusmn; ' + fmtInt(r.p.rhoSd) : '—'}</td>` +
      `<td>${r.p.area ? fmtGt(r.res.p10) : '—'}</td>` +
      `<td><b>${r.p.area ? fmtGt(r.res.p50) : '—'}</b></td>` +
      `<td>${r.p.area ? fmtGt(r.res.p90) : '—'}</td>` +
      `<td class="psrc">${r.pub ? r.pub[1].toFixed(2) : 'none'}</td>`;
    tr.addEventListener('click', () => selectAquifer(r.id));
    tb.appendChild(tr);
  }

  const agg = state.aggregate;
  const sum = k => rows.reduce((s, r) => s + (r.p.area ? r.res[k] : 0), 0);
  const sumP10 = sum('p10'), sumP50 = sum('p50'), sumP90 = sum('p90');

  tf.innerHTML =
    `<tr><td>Sum of rows <span class="psrc">(as Table 2)</span></td>` +
    `<td>${fmtInt(rows.reduce((s,r)=>s+r.p.area,0))}</td>` +
    `<td colspan="4" class="psrc" style="text-align:left">column totals</td>` +
    `<td>${fmtGt(sumP10)}</td><td><b>${fmtGt(sumP50)}</b></td><td>${fmtGt(sumP90)}</td>` +
    `<td class="psrc">${PUBLISHED_TOTAL.p50}</td></tr>` +
    `<tr><td>Basin aggregate <span class="psrc">(summed per trial)</span></td>` +
    `<td class="psrc">&mdash;</td>` +
    `<td colspan="4" class="psrc" style="text-align:left">percentiles of the summed distribution</td>` +
    `<td>${fmtGt(agg.p10)}</td><td><b>${fmtGt(agg.p50)}</b></td><td>${fmtGt(agg.p90)}</td>` +
    `<td class="psrc">&mdash;</td></tr>`;

  document.getElementById('aggNote').innerHTML =
    `<b>Why the two total rows differ, and which one to quote.</b>
     <i>Sum of rows</i> adds each percentile column straight down, which is how Table 2 of the paper
     reports the basin figure &mdash; here
     P10 ${fmtGt(sumP10)} &middot; P50 ${fmtGt(sumP50)} &middot; P90 ${fmtGt(sumP90)} Gt against the
     published P10 ${PUBLISHED_TOTAL.p10} &middot; P50 ${PUBLISHED_TOTAL.p50} &middot;
     P90 ${PUBLISHED_TOTAL.p90} Gt. <i>Basin aggregate</i> instead adds the seven aquifers within each
     Monte Carlo trial and then takes percentiles of that total:
     P10 ${fmtGt(agg.p10)} &middot; P50 ${fmtGt(agg.p50)} &middot; P90 ${fmtGt(agg.p90)} Gt.
     The aggregate is the statistically correct one if the aquifers are treated as independent, and it is
     both <b>narrower at the tails</b> and <b>higher at the median</b>: adding percentiles assumes every
     aquifer lands at its extreme at once, while summing trials lets independent errors cancel and pulls
     the total towards the sum of the means. Quote the row that matches your assumption &mdash; the
     column sum for comparability with the paper, the aggregate if you believe the aquifer uncertainties
     are genuinely independent.`;
}

function computeBasin() {
  const rows = [];
  const per = [];
  for (const a of AQUIFERS) {
    const p = {
      area:a.area, h:a.h, hSd:a.hSd, ntg:a.ntg, ntgSd:a.ntgSd,
      por:a.por, porSd:a.porSd, rho:a.rho, rhoSd:a.rhoSd,
      swirr:state.global.swirr, swirrSd:state.global.swirrSd,
      eff:state.global.eff, effSd:state.global.effSd
    };
    /* the aquifer currently on screen uses the edited parameters */
    if (a.id === state.aquifer) Object.assign(p, effectiveParams());
    const res = p.area ? runMC(p, state.trials, state.seed + a.id.charCodeAt(0)) : null;
    rows.push({ id:a.id, age:a.age, p, pub:a.pub, res: res || { p10:0,p50:0,p90:0,mean:0,samples:[],sorted:[] } });
    if (res) per.push(res.samples);
  }
  /* true aggregate: sum realisation by realisation */
  const n = state.trials;
  const tot = new Float64Array(n);
  for (const s of per) for (let i = 0; i < n; i++) tot[i] += s[i];
  const ts = Float64Array.from(tot).sort();
  state.aggregate = { p10: pctOf(ts, 0.90), p50: pctOf(ts, 0.50), p90: pctOf(ts, 0.10) };
  state.basin = rows;
}

let recomputeTimer = null;
function recompute(debounce) {
  clearTimeout(recomputeTimer);
  const go = () => {
    const p = effectiveParams();
    state.result = p.area ? runMC(p, state.trials, state.seed) : { p10:0,p50:0,p90:0,mean:0,det:0,samples:[],sorted:[] };
    computeBasin();
    renderStats();
    drawHistogram();
    drawTornado();
    renderBasin();
    drawComparison();
    document.getElementById('runInfo').textContent =
      fmtInt(state.trials) + ' trials · seed ' + state.seed;
  };
  if (debounce) recomputeTimer = setTimeout(go, 180); else go();
}

function renderTabs() {
  const tabs = document.querySelector('.tabs');
  tabs.querySelectorAll('.tab').forEach(t => t.remove());
  for (const a of AQUIFERS) {
    const b = document.createElement('button');
    b.className = 'tab';
    b.type = 'button';
    b.setAttribute('role', 'tab');
    b.setAttribute('aria-selected', String(a.id === state.aquifer));
    b.innerHTML = `<span class="dot${a.area ? '' : ' none'}"></span>Group ${a.id}`;
    b.title = a.age + (a.area ? '' : ' — no optimal zones');
    b.addEventListener('click', () => selectAquifer(a.id));
    tabs.appendChild(b);
  }
}

function selectAquifer(id) {
  state.aquifer = id;
  state.params = defaultParams(id);
  state.selArea = null; state.selLayer = null; state.cluster = null; state.rect = null;
  renderTabs();
  const a = AQUIFERS.find(x => x.id === id);
  document.getElementById('mapTitle').textContent = `Group ${id} — optimal injection zones`;
  const nCl = (window.MALAY_GRIDS[id] && window.MALAY_GRIDS[id].clusters.length) || 0;
  document.getElementById('mapNote').innerHTML = a.noOptimal
    ? `${a.age}. No optimal zones: modelled CO₂ density stays below the 300 kg/m³ cut-off at these depths.`
    : `${a.age}. ${nCl} connected optimal zone${nCl === 1 ? '' : 's'} mapped at 200 m resolution. ` +
      `Switch to <b>Pick cluster</b> and click one, or <b>Draw box</b> to measure a sub-area.`;

  loadGrid(id).then(g => {
    state.grid = g;
    drawMap();
    readoutEl.innerHTML = readoutBase();
  });
  renderParams();
  recompute();
}

/* --- controls --- */
document.querySelectorAll('.seg [data-mode]').forEach(btn => {
  btn.addEventListener('click', () => {
    state.mode = btn.dataset.mode;
    document.querySelectorAll('.seg [data-mode]').forEach(b =>
      b.setAttribute('aria-pressed', String(b === btn)));
    mapCv.style.cursor = state.mode === 'rect' ? 'crosshair' : state.mode === 'clus' ? 'pointer' : 'grab';
    document.getElementById('mapHint').textContent =
      state.mode === 'rect' ? 'Drag a box over optimal ground'
      : state.mode === 'clus' ? 'Click a green zone to select it'
      : 'Scroll to zoom · drag to pan';
    if (state.mode === 'all') clearSelection();
  });
});

document.getElementById('clearSel').addEventListener('click', clearSelection);
document.getElementById('resetParams').addEventListener('click', () => {
  state.params = defaultParams(state.aquifer);
  state.global = Object.assign({}, GLOBAL_DEFAULT);
  clearSelection();
  renderParams();
  recompute();
});
document.getElementById('runBtn').addEventListener('click', () => {
  state.seed = (state.seed + 1) | 0;
  document.getElementById('seed').value = state.seed;
  recompute();
});
document.getElementById('trials').addEventListener('change', e => {
  state.trials = parseInt(e.target.value, 10); recompute();
});
document.getElementById('seed').addEventListener('input', e => {
  const v = parseInt(e.target.value, 10);
  if (isFinite(v)) { state.seed = v; recompute(true); }
});
readoutEl.addEventListener('click', e => {
  if (e.target.id === 'useGrid') {
    e.preventDefault();
    state.params.area = state.grid.meta.optimalAreaKm2;
    renderParams(); recompute();
  }
});

document.getElementById('zIn').addEventListener('click', () => {
  const r = mapCv.getBoundingClientRect();
  zoomAbout(state.cam.cx, state.cam.cy, 1.5, r.width/2, r.height/2, r.width, r.height);
});
document.getElementById('zOut').addEventListener('click', () => {
  const r = mapCv.getBoundingClientRect();
  zoomAbout(state.cam.cx, state.cam.cy, 1/1.5, r.width/2, r.height/2, r.width, r.height);
});
document.getElementById('zRst').addEventListener('click', () => { resetCam(); drawMap(); });

/* theme */
function applyTheme(t) {
  if (t) document.documentElement.setAttribute('data-theme', t);
  else document.documentElement.removeAttribute('data-theme');
  const dark = t === 'dark' ||
    (!t && window.matchMedia('(prefers-color-scheme: dark)').matches);
  document.getElementById('themeBtn').textContent = dark ? 'Light mode' : 'Dark mode';
  /* Only the base raster is theme-dependent; the selection highlight uses the
     optimal-zone green, which is identical in both modes. */
  for (const k in gridCache) { gridCache[k].base = null; gridCache[k].baseTheme = null; }
  drawMap();
  drawHistogram(); drawTornado(); drawComparison();
}
document.getElementById('themeBtn').addEventListener('click', () => {
  const cur = document.documentElement.getAttribute('data-theme');
  const dark = cur === 'dark' || (!cur && window.matchMedia('(prefers-color-scheme: dark)').matches);
  const next = dark ? 'light' : 'dark';
  try { localStorage.setItem('mb-theme', next); } catch (_) {}
  applyTheme(next);
});

/* CSV export */
document.getElementById('exportBtn').addEventListener('click', () => {
  const L = [];
  L.push('# Malay Basin CO2 storage capacity - de Jonge-Anderson et al. (2025), IJGGC 143, 104347');
  L.push('# Generated ' + new Date().toISOString());
  L.push(`# trials=${state.trials} seed=${state.seed} Swirr=${state.global.swirr}+-${state.global.swirrSd} E=${state.global.eff}+-${state.global.effSd}`);
  L.push('# Percentiles use the exceedance convention: P10 = high estimate, P90 = low estimate.');
  L.push(['Group','Age','Area_km2','Thickness_m','Thickness_sd','NTG','NTG_sd','Porosity_pct','Porosity_sd',
          'CO2_density_kgm3','CO2_density_sd','P10_Gt','P50_Gt','P90_Gt','Mean_Gt','Published_P10','Published_P50','Published_P90'].join(','));
  for (const r of state.basin) {
    L.push([r.id, r.age, r.p.area, r.p.h, r.p.hSd, r.p.ntg, r.p.ntgSd, r.p.por, r.p.porSd, r.p.rho, r.p.rhoSd,
      r.p.area ? r.res.p10.toFixed(4) : '', r.p.area ? r.res.p50.toFixed(4) : '',
      r.p.area ? r.res.p90.toFixed(4) : '', r.p.area ? r.res.mean.toFixed(4) : '',
      r.pub ? r.pub[0] : '', r.pub ? r.pub[1] : '', r.pub ? r.pub[2] : ''].join(','));
  }
  const a = state.aggregate;
  L.push(['BASIN_TOTAL','','','','','','','','','','',
    a.p10.toFixed(4), a.p50.toFixed(4), a.p90.toFixed(4), '',
    PUBLISHED_TOTAL.p10, PUBLISHED_TOTAL.p50, PUBLISHED_TOTAL.p90].join(','));
  const blob = new Blob([L.join('\n')], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url; link.download = 'malay_basin_capacity.csv';
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

/* static prose that depends on the constants */
document.getElementById('calNote').innerHTML =
  `The paper states the <i>means</i> of storage efficiency (E = 2 %) and irreducible water saturation
   (S<sub>wirr</sub> = 27 %) but not their standard deviations, which the Monte Carlo needs.
   This tool defaults to <b>E = 0.02 &plusmn; 0.01</b> and <b>S<sub>wirr</sub> = 0.27 &plusmn; 0.05</b>.
   The efficiency spread was chosen by fitting: it is the value that best reproduces the published
   P10&ndash;P50&ndash;P90 of Table 2 across all seven aquifers. S<sub>wirr</sub> has almost no effect on
   the spread. Both are editable above, and neither should be treated as a published value.`;

document.getElementById('reproNote').innerHTML =
  `<p>With those defaults and the Table 2 inputs, this implementation returns a basin total of roughly
   P10 23.5 &middot; P50 7.8 &middot; P90 1.5 Gt against the published
   P10 ${PUBLISHED_TOTAL.p10} &middot; P50 ${PUBLISHED_TOTAL.p50} &middot; P90 ${PUBLISHED_TOTAL.p90} Gt
   (arithmetic sums, as the paper reports them). Per-aquifer values agree to within a few per cent.
   Residual differences come from the unstated standard deviations, the random seed, and the fact that
   the paper used 1,000 trials per aquifer &mdash; at that sample size the tails move by several per cent
   run to run. Set <b>Trials</b> to 1,000 to see that variability.</p>`;

document.getElementById('areaNote').innerHTML =
  `<b>The archived rasters do not reproduce every Table 2 area.</b> Measuring the optimal-zone pixels
   directly in <code>9 Optimal zone grids/*.bil</code> gives
   D 3,475 &middot; E 9,118 &middot; F 4,467 &middot; H 9,524 &middot; I 24,402 &middot; J 15,465 &middot; K 12,597 km&#178;,
   against Table 2's 3,348 &middot; 13,894 &middot; 18,108 &middot; 22,290 &middot; 24,924 &middot; 12,898 &middot; 10,643 km&#178;.
   Groups D and I match; E, F and H are substantially smaller in the grids, and J and K somewhat larger.
   Porosity and CO&#8322; density averaged inside those same grid zones reproduce Table 2 exactly for all seven
   aquifers, so the grids and this decoding are sound &mdash; the difference is in zone extent alone, most likely
   a consequence of the halved grid resolution noted in the dataset README, or of the rasters being
   regenerated after the table was finalised.
   <b>The calculator therefore defaults to the published Table 2 areas</b>; the grid-derived area is offered
   beneath the map, and any cluster or box you select on the map is measured from the grid.`;

/* boot */
(function init() {
  let saved = null;
  try { saved = localStorage.getItem('mb-theme'); } catch (_) {}
  if (saved) document.documentElement.setAttribute('data-theme', saved);
  const attr = document.documentElement.getAttribute('data-theme');
  const dark = attr === 'dark' ||
    (!attr && window.matchMedia('(prefers-color-scheme: dark)').matches);
  document.getElementById('themeBtn').textContent = dark ? 'Light mode' : 'Dark mode';

  bindTip(document.getElementById('hist'),    document.getElementById('histTip'), (x,y,w,h) => histHit(x,y,w,h));
  bindTip(document.getElementById('tornado'), document.getElementById('tornTip'), (x,y,w,h) => tornHit(x,y,w,h));
  bindTip(document.getElementById('cmp'),     document.getElementById('cmpTip'),  (x,y,w,h) => cmpHit(x,y,w,h));

  resetCam();
  mapCv.style.cursor = 'grab';
  selectAquifer('I');

  let rt = null;
  window.addEventListener('resize', () => {
    clearTimeout(rt);
    rt = setTimeout(() => { drawMap(); drawHistogram(); drawTornado(); drawComparison(); }, 120);
  });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    if (!document.documentElement.getAttribute('data-theme')) applyTheme(null);
  });
})();
