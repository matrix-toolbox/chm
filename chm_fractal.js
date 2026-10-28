// =============================================================================
// chm_fractal.js -- interactive explorer for the finite-time basin sections
// rendered by chm_fractal_fast.cpp. Loaded by chm_fractal.html.
//
// Drag a square to zoom in, right-click to step out. The selection always has
// the aspect ratio of the square canvas, so every view is a square window and
// pixels stay square. Each new view logs the matching chm_fractal command
// line to the browser console.
//
// KERNEL
// chmKernel() transcribes classify_fast<> (plane and kappa=1 line variants)
// and the per-pixel scalar of the C++ program. JavaScript numbers are IEEE
// doubles with correctly rounded + - * /, and every expression keeps the C++
// evaluation order, so samples agree with a C++ build that does not fuse
// multiply-adds (-ffp-contract=off). Reflection reuse is not ported: zoomed
// windows are never symmetric. The reference engine and --fit are not ported.
//
// THREADS AND COLOUR
// The kernel runs in Web Workers created from a Blob of its own source, so the
// page works from file:// with only chm_fractal.{html,css,js}. Without
// workers the same kernel runs on the main thread in small chunks. Workers
// return the scalar per pixel; the page colours it, so a palette change
// needs no new render.
//
// DEPTH LIMIT
// A window may not be narrower than MIN_REL_WIDTH * max(1, |kappa|, |c|),
// c its centre. At 1e-10 about 0.03% of the pixels of a 1300x1300 boundary
// view already differ between FMA and non-FMA builds of the C++ program;
// at 1e-14 the C++ program refuses because samples stop being distinct.
// =============================================================================
'use strict';

const MIN_REL_WIDTH = 1e-10;  // zoom limit, relative to the coordinate scale
const CACHE_BYTES = 160e6;    // finished frames kept for an instant step back
const MIN_DRAG = 4;           // image px; a shorter drag is a click, ignored

// Starting windows of the --figure presets for their kappa (k0 is the main
// figure). Any other kappa starts on [-R,R]^2, R = R_N max(1,|kappa|).
const FIGURES = [
    { kre: 0, kim: 0, window: [-2.8, 2.8, -2.8, 2.8] },
    { kre: 0, kim: 1, window: [-3.2, 3.2, -3.2, 3.2] },
    { kre: 1, kim: 0, window: [-3.4, 3.4, -3.4, 3.4] },
];

// Colour maps of the C++ scalar: basin of e in [0.04,0.34], basin of u in
// [0.68,0.98], unresolved 0.5. Stops are evenly spaced on [0,1]. With the 17
// Viridis stops and the 2 grey stops paint() reproduces the bytes of the C++
// PPM and --bw PGM. Turbo, Twilight and Cubehelix are matplotlib's maps at
// 33 stops.
const PALETTES = [
    { id: 'viridis', name: 'viridis', stops: [
        [ 68,   1,  84], [ 72,  24, 106], [ 71,  45, 123], [ 66,  64, 134],
        [ 59,  82, 139], [ 51,  99, 141], [ 44, 114, 142], [ 38, 130, 142],
        [ 33, 145, 140], [ 31, 160, 136], [ 40, 174, 128], [ 63, 188, 115],
        [ 94, 201,  98], [132, 212,  75], [173, 220,  48], [216, 226,  25],
        [253, 231,  37],
    ] },
    // matplotlib's terrain, sampled at its breakpoints: sea and shore for e, earth to snow for u
    { id: 'terrain', name: 'terrain', stops: [
        [ 51,  51, 153], [ 43,  67, 169], [ 35,  83, 185], [ 26, 102, 204], [ 18, 118, 220], [  8, 136, 238], [  0, 152, 254], [  0, 164, 220],
        [  0, 178, 178], [  0, 190, 143], [  1, 204, 102], [ 25, 209, 107], [ 49, 214, 112], [ 77, 219, 117], [101, 224, 122], [129, 230, 128],
        [153, 235, 133], [177, 239, 137], [205, 245, 143], [229, 250, 148], [254, 254, 152], [242, 238, 146], [230, 223, 139], [216, 205, 132],
        [204, 190, 125], [190, 172, 118], [178, 156, 111], [166, 141, 105], [152, 123,  97], [140, 108,  91], [129,  94,  86], [141, 109, 102],
        [153, 124, 118], [167, 142, 137], [179, 158, 153], [193, 176, 172], [205, 191, 188], [217, 206, 204], [231, 224, 223], [243, 240, 239],
        [255, 255, 255],
    ] },
    // a woodblock print: e on cream paper, u in indigo ink
    { id: 'hokusai', name: 'hokusai', stops: [
        [252, 248, 236], [242, 233, 212], [226, 212, 182], [205, 189, 156], [128, 146, 164],
        [ 52,  92, 140], [ 30,  62, 112], [ 18,  40,  80], [ 10,  22,  48],
    ] },
    // e in emeralds, u in rubies, a dark seam between them
    { id: 'jade', name: 'jade & ruby', stops: [
        [204, 255, 224], [ 72, 206, 138], [ 12, 136,  84], [  6,  64,  42], [ 12,  12,  14],
        [ 72,   6,  28], [164,  16,  52], [232,  72, 102], [255, 204, 214],
    ] },
    { id: 'turbo', name: 'turbo', stops: [
        [ 48,  18,  59], [ 57,  42, 115], [ 64,  64, 162], [ 68,  86, 199], [ 70, 107, 227], [ 70, 128, 246], [ 66, 148, 255], [ 55, 168, 250],
        [ 40, 188, 235], [ 28, 205, 216], [ 24, 221, 194], [ 31, 233, 175], [ 50, 242, 152], [ 78, 249, 125], [109, 254,  98], [139, 255,  75],
        [164, 252,  60], [185, 246,  53], [205, 236,  52], [223, 223,  55], [238, 207,  58], [248, 190,  57], [253, 172,  52], [254, 150,  43],
        [251, 126,  33], [244, 102,  23], [235,  80,  14], [223,  63,   8], [208,  47,   5], [190,  33,   2], [169,  22,   1], [146,  11,   1],
        [122,   4,   3],
    ] },
    { id: 'twilight', name: 'twilight', stops: [
        [226, 217, 226], [215, 215, 221], [196, 206, 212], [172, 194, 204], [149, 181, 199], [129, 166, 195], [114, 151, 193], [104, 135, 190],
        [ 98, 118, 186], [ 95, 100, 181], [ 94,  81, 173], [ 93,  61, 161], [ 89,  42, 143], [ 81,  27, 119], [ 69,  19,  92], [ 56,  17,  69],
        [ 47,  20,  54], [ 58,  17,  58], [ 74,  19,  66], [ 95,  23,  74], [116,  30,  79], [135,  39,  80], [152,  53,  80], [166,  69,  80],
        [178,  86,  82], [187, 105,  88], [194, 124,  99], [200, 144, 115], [204, 163, 137], [209, 182, 163], [216, 199, 190], [223, 212, 214],
        [226, 217, 226],
    ] },
    { id: 'cubehelix', name: 'cubehelix', stops: [
        [  0,   0,   0], [ 12,   5,  14], [ 21,  11,  29], [ 25,  19,  45], [ 27,  30,  59], [ 25,  42,  70], [ 23,  55,  77], [ 21,  70,  78],
        [ 22,  83,  76], [ 26,  96,  70], [ 35, 106,  62], [ 49, 114,  54], [ 68, 119,  49], [ 90, 122,  47], [114, 123,  50], [138, 122,  59],
        [161, 121,  74], [181, 121,  94], [197, 122, 118], [207, 126, 143], [212, 132, 169], [212, 141, 193], [209, 152, 212], [204, 166, 228],
        [198, 180, 238], [194, 195, 242], [193, 209, 243], [196, 221, 242], [203, 232, 240], [214, 240, 239], [227, 246, 240], [242, 251, 246],
        [255, 255, 255],
    ] },
    // e in ice blues, u in embers, a dark seam between them
    { id: 'icefire', name: 'ice & fire', stops: [
        [242, 251, 255], [159, 216, 240], [ 59, 143, 208], [ 31,  63, 143], [ 13,  13,  20],
        [110,  20,  32], [198,  59,  31], [243, 154,  60], [255, 241, 196],
    ] },
    { id: 'gold', name: 'gold leaf', stops: [
        [ 11,  10,  15], [ 28,  20,  22], [ 58,  36,  26], [ 98,  60,  30], [140,  90,  38],
        [184, 128,  52], [219, 168,  72], [242, 210, 122], [255, 247, 222],
    ] },
    // e from night sky to teal, u from magenta to pale rose
    { id: 'aurora', name: 'aurora', stops: [
        [  4,   6,  26], [  8,  28,  62], [ 12,  78, 104], [ 18, 150, 140], [120, 230, 180],
        [150,  60, 170], [205,  80, 180], [240, 150, 210], [253, 233, 246],
    ] },
    { id: 'ocean', name: 'ocean', stops: [
        [  2,  10,  30], [  4,  30,  66], [  8,  54, 100], [ 14,  82, 130], [ 24, 114, 156],
        [ 40, 148, 176], [ 78, 186, 192], [150, 222, 214], [233, 250, 245],
    ] },
    { id: 'grey', name: 'grey', stops: [[0, 0, 0], [255, 255, 255]] },
];

// viridis() of the C++ program for any evenly spaced stop table
function paint(v, stops, out, o) {
    let t = v;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    const x = t * (stops.length - 1);
    const i = Math.min(stops.length - 2, Math.trunc(x));
    const f = x - i;
    const p = stops[i], q = stops[i + 1];
    out[o]     = (p[0] + f * (q[0] - p[0]) + 0.5) | 0;
    out[o + 1] = (p[1] + f * (q[1] - p[1]) + 0.5) | 0;
    out[o + 2] = (p[2] + f * (q[2] - p[2]) + 0.5) | 0;
}

// ------------------------------------------------------------------ kernel ---
// Self-contained: its source is also shipped to the workers.
function chmKernel() {
    // coordinates() without reflection: sample j of count on [lo, hi]
    function coordinates(lo, hi, count, endpoint) {
        const x = new Float64Array(count);
        for (let j = 0; j < count; ++j) {
            const t = endpoint ? (count > 1 ? j / (count - 1) : 0.5) : (j + 0.5) / count;
            x[j] = lo + (hi - lo) * t;
        }
        return x;
    }

    let P = null;   // per-job constants, rebuilt when the job id changes
    function prepare(job) {
        if (P && P.id === job.id) return;
        const m = job.m, n = job.size * job.ss, unit = job.steps === 2 ? 1 : 2;
        const tableSize = Math.min(job.maxit, 65535) + 1;
        const band = new Float64Array(tableSize);
        for (let k = 0; k < tableSize; ++k) band[k] = 0.5 + 0.5 * Math.cos(0.9 * k);
        P = {
            id: job.id,
            m, N: 2 * m + 1, mm: m * m, m2: 2.0 * m, m2m1: 2.0 * m - 1.0, hm1: 0.5 * (m - 1.0),
            tol2: job.tol * job.tol, maxit: job.maxit, unit, limit: Math.floor(job.maxit / unit),
            kre: job.kre, kim: job.kim, line: job.kre === 1 && job.kim === 0,
            xs: coordinates(job.xi, job.xf, n, job.endpoint),
            ys: coordinates(job.yi, job.yf, n, job.endpoint),
            band, tableSize,
        };
    }

    let hit = 0;   // first-arrival step of the last classification

    // classify_fast<false>: [z:1:kappa] under Q. 0 unresolved, 1 basin of e, 2 basin of u.
    function classifyPlane(zr, zi) {
        const m = P.m, mm = P.mm, m2 = P.m2, hm1 = P.hm1, tol2 = P.tol2, N = P.N;
        const unit = P.unit, limit = P.limit;
        let ar = zr, ai = zi, br = 1.0, bi = 0.0, cr = P.kre, ci = P.kim;
        let s = Math.max(Math.abs(ar), Math.abs(ai), Math.abs(br), Math.abs(bi), Math.abs(cr), Math.abs(ci));
        if (!(s > 0.0 && s < Infinity)) { hit = 0; return 0; }
        let inv = 1.0 / s;
        ar *= inv; ai *= inv; br *= inv; bi *= inv; cr *= inv; ci *= inv;
        for (let n = 0; ; ++n) {
            // close_to_attractor
            const transverse = m * ((br * br + bi * bi) + (cr * cr + ci * ci));
            const den = (ar * ar + ai * ai) + transverse;
            if (den > 0.0 && den < Infinity) {
                const threshold = tol2 * den;
                if (transverse < threshold) { hit = n * unit; return 1; }
                let dr = ar - br, di = ai - bi;
                const ab = dr * dr + di * di;
                dr = ar - cr; di = ai - ci;
                const ac = dr * dr + di * di;
                dr = br - cr; di = bi - ci;
                const bc = dr * dr + di * di;
                if (m * (ab + ac) + mm * bc < N * threshold) { hit = n * unit; return 2; }
            }
            if (n === limit) { hit = P.maxit; return 0; }
            if (n) {
                s = Math.max(Math.abs(ar), Math.abs(ai), Math.abs(br), Math.abs(bi), Math.abs(cr), Math.abs(ci));
                if (!(s > 0.0 && s < Infinity)) { hit = n * unit; return 0; }
                inv = 1.0 / s;
                ar *= inv; ai *= inv; br *= inv; bi *= inv; cr *= inv; ci *= inv;
            }
            // quartic: Q = Phi(Phi(x))/N
            const aar = ar * ar - ai * ai, aai = 2.0 * ar * ai;
            const bbr = br * br - bi * bi, bbi = 2.0 * br * bi;
            const ccr = cr * cr - ci * ci, cci = 2.0 * cr * ci;
            const sr = bbr + ccr, si = bbi + cci;
            const comr = hm1 * (sr * sr - si * si), comi = hm1 * (2.0 * sr * si);
            ar = (aar * aar - aai * aai) + m2 * (bbr * ccr - bbi * cci);
            ai = 2.0 * aar * aai + m2 * (bbr * cci + bbi * ccr);
            br = 2.0 * (aar * ccr - aai * cci) + (bbr * bbr - bbi * bbi) + comr;
            bi = 2.0 * (aar * cci + aai * ccr) + 2.0 * bbr * bbi + comi;
            cr = 2.0 * (aar * bbr - aai * bbi) + (ccr * ccr - cci * cci) + comr;
            ci = 2.0 * (aar * bbi + aai * bbr) + 2.0 * ccr * cci + comi;
        }
    }

    // classify_fast<true>: kappa = 1 exactly, Q restricted to [a:b:b].
    function classifyLine(zr, zi) {
        const m2 = P.m2, m2m1 = P.m2m1, tol2 = P.tol2, N = P.N;
        const unit = P.unit, limit = P.limit;
        let ar = zr, ai = zi, br = 1.0, bi = 0.0;
        let s = Math.max(Math.abs(ar), Math.abs(ai), Math.abs(br), Math.abs(bi));
        if (!(s > 0.0 && s < Infinity)) { hit = 0; return 0; }
        let inv = 1.0 / s;
        ar *= inv; ai *= inv; br *= inv; bi *= inv;
        for (let n = 0; ; ++n) {
            // close_on_line
            const transverse = m2 * (br * br + bi * bi);
            const den = (ar * ar + ai * ai) + transverse;
            if (den > 0.0 && den < Infinity) {
                const threshold = tol2 * den;
                if (transverse < threshold) { hit = n * unit; return 1; }
                const dr = ar - br, di = ai - bi;
                if (m2 * (dr * dr + di * di) < N * threshold) { hit = n * unit; return 2; }
            }
            if (n === limit) { hit = P.maxit; return 0; }
            if (n) {
                s = Math.max(Math.abs(ar), Math.abs(ai), Math.abs(br), Math.abs(bi));
                if (!(s > 0.0 && s < Infinity)) { hit = n * unit; return 0; }
                inv = 1.0 / s;
                ar *= inv; ai *= inv; br *= inv; bi *= inv;
            }
            // quartic_line
            const aar = ar * ar - ai * ai, aai = 2.0 * ar * ai;
            const bbr = br * br - bi * bi, bbi = 2.0 * br * bi;
            const b4r = bbr * bbr - bbi * bbi, b4i = 2.0 * bbr * bbi;
            ar = (aar * aar - aai * aai) + m2 * b4r;
            ai = 2.0 * aar * aai + m2 * b4i;
            br = 2.0 * (aar * bbr - aai * bbi) + m2m1 * b4r;
            bi = 2.0 * (aar * bbi + aai * bbr) + m2m1 * b4i;
        }
    }

    // The C++ scalar of canvas rows r0..r1-1 (top down; +Im upwards), with
    // sample counts as in --stats.
    function renderRows(job, r0, r1) {
        prepare(job);
        const W = job.size, ss = job.ss, xs = P.xs, ys = P.ys, band = P.band, tableSize = P.tableSize;
        const classify = P.line ? classifyLine : classifyPlane;
        const invSS2 = 1.0 / (ss * ss);
        const values = new Float64Array((r1 - r0) * W);
        let cnt0 = 0, cnt1 = 0, cnt2 = 0, maxhit = 0, o = 0;
        for (let r = r0; r < r1; ++r) {
            const j = W - 1 - r;
            for (let i = 0; i < W; ++i) {
                let acc = 0.0;
                for (let sy = 0; sy < ss; ++sy) {
                    const y = ys[j * ss + sy];
                    for (let sx = 0; sx < ss; ++sx) {
                        const lab = classify(xs[i * ss + sx], y);
                        if (lab === 0) { ++cnt0; acc += 0.5; continue; }
                        if (lab === 1) ++cnt1; else ++cnt2;
                        if (hit > maxhit) maxhit = hit;
                        const v = hit < tableSize ? band[hit] : 0.5 + 0.5 * Math.cos(0.9 * hit);
                        acc += lab === 1 ? 0.04 + 0.30 * v : 0.68 + 0.30 * v;
                    }
                }
                values[o++] = acc * invSS2;
            }
        }
        return { values, cnt0, cnt1, cnt2, maxhit };
    }

    return { renderRows };
}

// Worker entry point; runs inside the Blob worker next to chmKernel.
function chmWorker() {
    const kernel = chmKernel();
    self.onmessage = (event) => {
        const { job, r0, r1 } = event.data;
        const out = kernel.renderRows(job, r0, r1);
        self.postMessage(Object.assign({ id: job.id, r0, r1 }, out), [out.values.buffer]);
    };
}

// The polynomial R^4-2R^3-(N-2)R-(N-1) has one positive root; as natural_radius().
function naturalRadius(N) {
    const poly = (t) => t * t * t * (t - 2.0) - (N - 2) * t - (N - 1);
    let lo = 0.0, hi = 4.0;
    while (poly(hi) <= 0.0) hi *= 2.0;
    for (let i = 0; i < 100; ++i) {
        const mid = lo + (hi - lo) * 0.5;
        if (mid === lo || mid === hi) break;
        if (poly(mid) > 0.0) hi = mid; else lo = mid;
    }
    return hi;
}

// --------------------------------------------------------------- scheduler ---
// Rows are handed out in chunks to whichever lane is idle, like
// schedule(dynamic). Results of a superseded job are dropped.
function createPool() {
    const pool = { mode: 'workers', size: 0, start, cancel };
    let lanes = [], answered = false, watchdog = 0;
    let job = null, nextRow = 0, rowsDone = 0, chunkRows = 1, onChunk = null, onDone = null;

    function chunkFor(j) {
        const perLane = Math.ceil(j.size / (8 * lanes.length));
        const bySamples = Math.floor(40000 / (j.size * j.ss * j.ss));
        return Math.max(1, Math.min(perLane, bySamples));
    }

    function feed(lane) {
        if (!job || nextRow >= job.size) return;
        const r0 = nextRow, r1 = Math.min(job.size, r0 + chunkRows);
        nextRow = r1;
        lane.busy = true;
        lane.post({ job, r0, r1 });
    }

    function handle(lane, msg) {
        answered = true;
        lane.busy = false;
        if (job && msg.id === job.id) {
            rowsDone += msg.r1 - msg.r0;
            onChunk(msg, rowsDone / job.size);
            if (rowsDone >= job.size) {
                const done = onDone;
                job = null;
                done();
            }
        }
        feed(lane);
    }

    function mainThreadLane() {
        const kernel = chmKernel();
        const lane = {
            busy: false,
            post(msg) {
                setTimeout(() => {
                    const out = kernel.renderRows(msg.job, msg.r0, msg.r1);
                    handle(lane, Object.assign({ id: msg.job.id, r0: msg.r0, r1: msg.r1 }, out));
                }, 0);
            },
        };
        return lane;
    }

    function useMainThread() {
        if (pool.mode === 'main thread') return;
        clearTimeout(watchdog);
        for (const lane of lanes) if (lane.worker) lane.worker.terminate();
        lanes = [mainThreadLane()];
        pool.mode = 'main thread';
        pool.size = 1;
        if (job) {
            nextRow = 0;
            rowsDone = 0;
            chunkRows = chunkFor(job);
            feed(lanes[0]);
        }
    }

    function start(j, chunkCallback, doneCallback) {
        job = j;
        onChunk = chunkCallback;
        onDone = doneCallback;
        nextRow = 0;
        rowsDone = 0;
        chunkRows = chunkFor(j);
        for (const lane of lanes) if (!lane.busy) feed(lane);
        // a worker that is silently blocked never answers
        if (!answered && pool.mode === 'workers' && !watchdog)
            watchdog = setTimeout(() => { if (!answered) useMainThread(); }, 4000);
    }

    function cancel() { job = null; }

    try {
        if (typeof Worker === 'undefined' || typeof Blob === 'undefined') throw new Error('no workers');
        const source = `'use strict';\n${chmKernel.toString()}\n(${chmWorker.toString()})();\n`;
        const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
        const n = Math.max(1, Math.min(16, navigator.hardwareConcurrency || 4));
        for (let k = 0; k < n; ++k) {
            const worker = new Worker(url);
            const lane = { busy: false, worker, post: (msg) => worker.postMessage(msg) };
            worker.onmessage = (event) => handle(lane, event.data);
            worker.onerror = (event) => { event.preventDefault(); useMainThread(); };
            lanes.push(lane);
        }
        pool.size = n;
    } catch (err) {
        useMainThread();
    }
    return pool;
}

// ------------------------------------------------------------ number text ---
// 1-2-5 tick step for a range, aiming at about `target` divisions; as nice_step().
function niceStep(range, target) {
    if (!(range > 0) || target < 1) return 0;
    const raw = range / target;
    if (!(raw > 0) || !Number.isFinite(raw)) return 0;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    if (!(mag > 0)) return raw;
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
}

function decimalsFor(step) {
    return Math.min(20, Math.max(0, Math.ceil(-Math.log10(step) - Math.log10(1.000001))));
}

function fixed(v, dec) {
    const s = v.toFixed(Math.min(20, Math.max(0, dec)));
    return /^-0(\.0*)?$/.test(s) ? s.slice(1) : s;
}

function magnitude(x) {
    if (!(x > 0) || !Number.isFinite(x)) return '—';
    return x < 1e4 ? String(Number(x.toPrecision(3))) : x.toExponential(1).replace('e+', 'e');
}

// C++ coordinates() refuses a spacing that does not change the endpoints.
function checkResolvable(lo, hi, count, endpoint) {
    const dx = (hi - lo) / (endpoint && count > 1 ? count - 1 : count);
    if (count > 1 && (!(dx > 0) || lo + dx === lo || hi - dx === hi))
        throw new Error('sample spacing is not resolvable in double precision');
}

// ---------------------------------------------------------------------- UI ---
function initExplorer() {
    const $ = (id) => document.getElementById(id);
    const stage = $('stage'), image = $('image'), overlay = $('overlay');
    const ctx = image.getContext('2d');
    const octx = overlay.getContext('2d');
    const form = $('controls');
    const field = (name) => form.elements.namedItem(name);
    const scratch = document.createElement('canvas');
    const pool = createPool();
    const PALETTE_KEY = 'chm_fractal.palette';

    const state = {
        params: null,     // parameters of the image on screen
        palette: PALETTES[0],
        grid: false,      // coordinate grid overlay
        home: null,      // starting window for the current kappa
        view: null,       // current window {xi, xf, yi, yf}
        stack: [],        // earlier windows {view, frame}
        frame: null,      // scalar field of the current view {key, size, values, filled, done}
        drag: null,
        pointer: 'mouse',
    };
    let jobId = 0;
    let pixels = null;    // ImageData the frames are painted into

    // ---------------------------------------------------------- parameters --
    function readParams() {
        const text = (name) => String(field(name).value).trim();
        const real = (name, label) => {
            const t = text(name), v = Number(t);
            if (t === '' || !Number.isFinite(v)) throw new Error(`${label}: expected a finite real number`);
            return v;
        };
        const integer = (name, label, lo, hi) => {
            const t = text(name);
            if (!/^[+-]?\d+$/.test(t)) throw new Error(`${label}: expected an integer`);
            const v = Number(t);
            if (v < lo || v > hi) throw new Error(`${label} must lie in [${lo}, ${hi}]`);
            return v;
        };
        const p = {
            m: integer('m', 'm', 1, 1073741823),
            kre: real('kre', 'Re κ'),
            kim: real('kim', 'Im κ'),
            steps: text('map') === 'phi' ? 1 : 2,
            maxit: integer('iter', 'nₘₐₓ', 0, 1000000),
            tol: real('tol', 'ε'),
            size: integer('size', 'size', 64, 4096),
            ss: integer('ss', 'supersampling', 1, 4),
            endpoint: text('sample') === 'endpoint',
        };
        if (!(p.tol > 0 && p.tol < 1) || p.tol * p.tol === 0)
            throw new Error('ε must satisfy 0 < ε < 1 with a representable square');
        return p;
    }

    // everything that changes the scalar field
    const pixelKey = (p) => [p.m, p.kre, p.kim, p.steps, p.maxit, p.tol, p.size, p.ss, p.endpoint].join('|');

    function setParams(p) {
        if (image.width !== p.size || image.height !== p.size) {
            image.width = p.size;
            image.height = p.size;
        }
        state.params = p;
        $('apply').classList.remove('dirty');
    }

    function homeFor(p) {
        const fig = FIGURES.find((f) => f.kre === p.kre && f.kim === p.kim);
        if (fig) {
            const [xi, xf, yi, yf] = fig.window;
            return { xi, xf, yi, yf };
        }
        const R = naturalRadius(2 * p.m + 1) * Math.max(1.0, Math.hypot(p.kre, p.kim));
        return { xi: -R, xf: R, yi: -R, yf: R };
    }

    const sameWindow = (a, b) => !!a && !!b && a.xi === b.xi && a.xf === b.xf && a.yi === b.yi && a.yf === b.yf;

    // ------------------------------------------------------------ geometry --
    // Continuous position along an axis (pixel k centred at k + 0.5) to the
    // fraction of the window, matching coordinates() and the --ss average.
    function fractionAt(pos, p) {
        if (!p.endpoint) return pos / p.size;
        const n = p.size * p.ss;
        return n > 1 ? ((pos - 0.5) * p.ss + (p.ss - 1) / 2) / (n - 1) : 0.5;
    }

    function positionOf(t, p) {
        if (!p.endpoint) return t * p.size;
        return (t * (p.size * p.ss - 1) - (p.ss - 1) / 2) / p.ss + 0.5;
    }

    function toComplex(X, Y, v, p) {
        return [v.xi + (v.xf - v.xi) * fractionAt(X, p),
                v.yi + (v.yf - v.yi) * fractionAt(p.size - Y, p)];
    }

    const centre = (v) => [0.5 * (v.xi + v.xf), 0.5 * (v.yi + v.yf)];
    const widthOf = (v) => Math.min(v.xf - v.xi, v.yf - v.yi);

    function minWidthAt(x, y, p) {
        return MIN_REL_WIDTH * Math.max(1, Math.hypot(p.kre, p.kim), Math.abs(x), Math.abs(y));
    }

    function atDepthLimit() {
        const [x, y] = centre(state.view);
        return widthOf(state.view) <= minWidthAt(x, y, state.params) * (1 + 1e-6);
    }

    // The selection square: anchored at the press point, growing toward the
    // pointer, kept inside the image and never below the depth limit.
    function dragBox(d) {
        const p = state.params, v = state.view, W = p.size;
        const sx = d.x1 >= d.x0 ? 1 : -1, sy = d.y1 >= d.y0 ? 1 : -1;
        const raw = Math.max(Math.abs(d.x1 - d.x0), Math.abs(d.y1 - d.y0));
        const place = (side) => {
            const left = sx > 0 ? d.x0 : d.x0 - side, top = sy > 0 ? d.y0 : d.y0 - side;
            return {
                left: Math.min(Math.max(left, 0), W - side),
                top: Math.min(Math.max(top, 0), W - side),
                side,
            };
        };
        let box = place(Math.min(raw, W));
        const [cx, cy] = toComplex(box.left + box.side / 2, box.top + box.side / 2, v, p);
        const sideMin = W * minWidthAt(cx, cy, p) / widthOf(v);
        box.limited = box.side < sideMin;
        if (box.limited) box = Object.assign(place(Math.min(sideMin, W)), { limited: true });
        box.blocked = sideMin >= 0.99 * W;   // less than 1% left to gain
        box.mag = W / box.side;
        return box;
    }

    function windowFromBox(box) {
        const p = state.params, v = state.view;
        const [xi, yf] = toComplex(box.left, box.top, v, p);
        const [xf, yi] = toComplex(box.left + box.side, box.top + box.side, v, p);
        const w = { xi, xf, yi, yf };
        const [cx, cy] = centre(w), mw = minWidthAt(cx, cy, p);
        // a limited box lands exactly on the limit, so atDepthLimit() sees it
        if (box.limited || widthOf(w) < mw) {
            w.xi = cx - 0.5 * mw; w.xf = cx + 0.5 * mw;
            w.yi = cy - 0.5 * mw; w.yf = cy + 0.5 * mw;
        }
        return w;
    }

    // ------------------------------------------------------------ painting --
    function paintRows(frame, r0, r1) {
        const W = frame.size, stops = state.palette.stops, v = frame.values;
        if (!pixels || pixels.width !== W) pixels = ctx.createImageData(W, W);
        const px = pixels.data;
        for (let k = r0 * W, o = 4 * r0 * W, end = r1 * W; k < end; ++k, o += 4) {
            paint(v[k], stops, px, o);
            px[o + 3] = 255;
        }
        ctx.putImageData(pixels, 0, 0, 0, r0, W, r1 - r0);
    }

    // every finished row, e.g. after a palette change
    function repaint(frame) {
        if (!frame || frame.size !== image.width) return;
        for (let r = 0; r < frame.size;) {
            if (!frame.filled[r]) { ++r; continue; }
            let e = r;
            while (e < frame.size && frame.filled[e]) ++e;
            paintRows(frame, r, e);
            r = e;
        }
    }

    // ----------------------------------------------------------- rendering --
    function render() {
        const p = state.params, v = state.view, n = p.size * p.ss;
        try {
            checkResolvable(v.xi, v.xf, n, p.endpoint);
            checkResolvable(v.yi, v.yf, n, p.endpoint);
        } catch (err) {
            toast(err.message, 'warn');
            return;
        }
        const job = {
            id: ++jobId, size: p.size, ss: p.ss, endpoint: p.endpoint,
            xi: v.xi, xf: v.xf, yi: v.yi, yf: v.yf,
            m: p.m, steps: p.steps, maxit: p.maxit, tol: p.tol, kre: p.kre, kim: p.kim,
        };
        const frame = {
            key: pixelKey(p), size: p.size, done: false,
            values: new Float64Array(p.size * p.size), filled: new Uint8Array(p.size),
        };
        state.frame = frame;
        setBusy(0);
        updateReadout();
        pool.start(job, (msg, done) => {
            frame.values.set(msg.values, msg.r0 * p.size);
            frame.filled.fill(1, msg.r0, msg.r1);
            paintRows(frame, msg.r0, msg.r1);
            setBusy(done);
        }, () => {
            frame.done = true;
            setBusy(null);
        });
    }

    function snapshot() {
        scratch.width = image.width;
        scratch.height = image.height;
        scratch.getContext('2d').drawImage(image, 0, 0);
        return scratch;
    }

    // stretch the selected part of the old image while the new one renders
    function previewZoom(box) {
        const W = state.params.size, src = snapshot();
        ctx.save();
        ctx.imageSmoothingEnabled = false;
        ctx.drawImage(src, box.left, box.top, box.side, box.side, 0, 0, W, W);
        ctx.restore();
    }

    // shrink the current image into its place inside the parent window
    function previewBack(parent) {
        const p = state.params, v = state.view, W = p.size, src = snapshot();
        const X = (x) => positionOf((x - parent.xi) / (parent.xf - parent.xi), p);
        const Y = (y) => W - positionOf((y - parent.yi) / (parent.yf - parent.yi), p);
        const left = X(v.xi), right = X(v.xf), top = Y(v.yf), bottom = Y(v.yi);
        ctx.save();
        ctx.fillStyle = '#1d1e23';
        ctx.fillRect(0, 0, W, W);
        ctx.drawImage(src, 0, 0, W, W, left, top, Math.max(1, right - left), Math.max(1, bottom - top));
        ctx.restore();
    }

    // keep the newest finished frames within CACHE_BYTES
    function trimCache() {
        let budget = CACHE_BYTES;
        for (let k = state.stack.length - 1; k >= 0; --k) {
            const f = state.stack[k].frame;
            if (!f) continue;
            budget -= f.values.byteLength;
            if (budget < 0) state.stack[k].frame = null;
        }
    }

    function zoomTo(box) {
        const view = windowFromBox(box);
        previewZoom(box);
        state.stack.push({ view: state.view, frame: state.frame && state.frame.done ? state.frame : null });
        trimCache();
        state.view = view;
        render();
    }

    function back() {
        const entry = state.stack.pop();
        if (!entry) { toast('Already at the starting view'); return; }
        const f = entry.frame;
        if (f && f.key === pixelKey(state.params)) {
            pool.cancel();
            state.view = entry.view;
            state.frame = f;
            repaint(f);
            setBusy(null);
            updateReadout();
        } else {
            previewBack(entry.view);
            state.view = entry.view;
            render();
        }
    }

    function goHome() {
        state.view = state.home;
        state.stack = [];
        render();
    }

    // ------------------------------------------------------------- palette --
    function gradientCss(stops) {
        const n = stops.length - 1;
        return `linear-gradient(90deg, ${stops.map((s, k) => `rgb(${s[0]} ${s[1]} ${s[2]}) ${(100 * k / n).toFixed(2)}%`).join(', ')})`;
    }

    function setPalette(id) {
        state.palette = PALETTES.find((pal) => pal.id === id) || PALETTES[0];
        for (const b of $('palettes').children)
            b.setAttribute('aria-pressed', String(b.dataset.id === state.palette.id));
        try { localStorage.setItem(PALETTE_KEY, state.palette.id); } catch (err) { /* storage unavailable */ }
        repaint(state.frame);
        if (state.params) updateReadout();
    }

    function buildPalettes() {
        const box = $('palettes');
        for (const pal of PALETTES) {
            const b = document.createElement('button');
            b.type = 'button';
            b.className = 'swatch';
            b.dataset.id = pal.id;
            b.title = pal.name;
            const bar = document.createElement('span');
            bar.className = 'bar';
            bar.style.background = gradientCss(pal.stops);
            const name = document.createElement('span');
            name.className = 'name';
            name.textContent = pal.name;
            b.append(bar, name);
            b.addEventListener('click', () => setPalette(pal.id));
            box.appendChild(b);
        }
    }

    // ------------------------------------------------------------- overlay --
    function ticks(g, lo, hi, S, fontPx, horizontal) {
        let out = [];
        for (const target of [8, 5, 3, 2]) {
            const step = niceStep(hi - lo, target);
            if (!(step > 0)) return [];
            const dec = decimalsFor(step), first = Math.ceil(lo / step);
            const count = Math.ceil((hi - lo) / step) + 2;
            out = [];
            for (let j = 0; j < count; ++j) {
                const val = (first + j) * step;
                if (val >= lo && val <= hi) out.push({ val, label: fixed(val, dec) });
            }
            const need = horizontal
                ? out.reduce((w, t) => Math.max(w, g.measureText(t.label).width), 0) + 2 * fontPx
                : 4 * fontPx;   // room for a label below one line and above the next
            if (need < S * step / (hi - lo)) break;
        }
        return out;
    }

    // k: canvas px per image px; lw: line width in canvas px
    function drawGrid(g, k, fontPx, lw) {
        const p = state.params, v = state.view, W = p.size, S = W * k;
        g.save();
        g.font = `${fontPx}px ui-monospace, "SF Mono", Menlo, Consolas, monospace`;
        g.textBaseline = 'top';
        const xt = ticks(g, v.xi, v.xf, S, fontPx, true);
        const yt = ticks(g, v.yi, v.yf, S, fontPx, false);
        const X = (val) => positionOf((val - v.xi) / (v.xf - v.xi), p) * k;
        const Y = (val) => (W - positionOf((val - v.yi) / (v.yf - v.yi), p)) * k;
        const snap = (z) => Math.round(z) + 0.5;

        g.beginPath();
        for (const t of xt) { const x = snap(X(t.val)); g.moveTo(x, 0); g.lineTo(x, S); }
        for (const t of yt) { const y = snap(Y(t.val)); g.moveTo(0, y); g.lineTo(S, y); }
        g.lineWidth = 3 * lw;
        g.strokeStyle = 'rgba(0, 0, 0, 0.35)';
        g.stroke();
        g.lineWidth = lw;
        g.strokeStyle = 'rgba(255, 255, 255, 0.8)';
        g.stroke();

        const pad = Math.round(fontPx * 0.45), h = fontPx;
        g.lineJoin = 'round';
        g.lineWidth = Math.max(2, fontPx / 4);
        g.strokeStyle = 'rgba(0, 0, 0, 0.8)';
        g.fillStyle = '#ffffff';
        const label = (text, x, y) => { g.strokeText(text, x, y); g.fillText(text, x, y); };
        const bandTop = S - h - 2 * pad;
        for (const t of yt) {
            let y = Y(t.val) + pad;
            if (y + h > bandTop) y = Y(t.val) - pad - h;
            label(t.label, pad, Math.max(0, y));
        }
        for (const t of xt) {
            const w = g.measureText(t.label).width;
            let x = X(t.val) + pad;
            if (x + w > S - 1) x = X(t.val) - pad - w;
            label(t.label, Math.max(pad, x), S - h - pad);
        }
        g.restore();
    }

    function drawBox(g, k, dpr) {
        const b = state.drag.box, S = state.params.size * k;
        const x = b.left * k, y = b.top * k, s = b.side * k;
        g.save();
        g.fillStyle = 'rgba(0, 0, 0, 0.38)';
        g.beginPath();
        g.rect(0, 0, S, S);
        g.rect(x, y, s, s);
        g.fill('evenodd');
        g.lineWidth = 3 * dpr;
        g.strokeStyle = 'rgba(0, 0, 0, 0.6)';
        g.strokeRect(x, y, s, s);
        g.lineWidth = 1.5 * dpr;
        g.strokeStyle = b.limited ? '#f0a640' : '#ffffff';
        if (b.limited) g.setLineDash([6 * dpr, 4 * dpr]);
        g.strokeRect(x, y, s, s);
        g.setLineDash([]);
        const text = b.blocked ? 'maximum depth' : `×${magnitude(b.mag)}${b.limited ? ' · depth limit' : ''}`;
        const fontPx = 12 * dpr;
        g.font = `600 ${fontPx}px system-ui, sans-serif`;
        g.textBaseline = 'bottom';
        const tw = g.measureText(text).width;
        const tx = Math.min(Math.max(x, 4 * dpr), S - tw - 4 * dpr);
        let ty = y - 4 * dpr;
        if (ty - fontPx < 0) ty = y + s + fontPx + 4 * dpr;
        if (ty > S - 2 * dpr) ty = y + fontPx + 6 * dpr;
        g.lineJoin = 'round';
        g.lineWidth = 3 * dpr;
        g.strokeStyle = 'rgba(0, 0, 0, 0.8)';
        g.strokeText(text, tx, ty);
        g.fillStyle = b.limited ? '#f0a640' : '#ffffff';
        g.fillText(text, tx, ty);
        g.restore();
    }

    function drawOverlay() {
        const dpr = window.devicePixelRatio || 1;
        const rect = stage.getBoundingClientRect();
        const w = Math.max(1, Math.round(rect.width * dpr));
        if (overlay.width !== w || overlay.height !== w) {
            overlay.width = w;
            overlay.height = w;
        }
        octx.clearRect(0, 0, w, w);
        if (!state.params) return;
        const k = w / state.params.size;
        if (state.grid) drawGrid(octx, k, Math.round(11 * dpr), dpr);
        if (state.drag && state.drag.box) drawBox(octx, k, dpr);
    }

    // ------------------------------------------------------------- readout --
    let toastTimer = 0;
    function toast(text, kind) {
        const t = $('toast');
        t.textContent = text;
        t.dataset.kind = kind || 'info';
        t.hidden = false;
        clearTimeout(toastTimer);
        toastTimer = setTimeout(() => { t.hidden = true; }, 3500);
    }

    function setBusy(done) {
        const bar = $('progress');
        if (done === null) {
            bar.hidden = true;
            document.body.dataset.state = 'done';
            return;
        }
        bar.hidden = false;
        bar.firstElementChild.style.width = `${(100 * done).toFixed(1)}%`;
        document.body.dataset.state = 'rendering';
    }

    function updateReadout() {
        const p = state.params, v = state.view, w = widthOf(v);
        const dec = Math.min(17, Math.max(2, Math.ceil(-Math.log10(w)) + 3));
        $('cursor').textContent = '—';
        $('window').textContent = `Re [${fixed(v.xi, dec)}, ${fixed(v.xf, dec)}]  Im [${fixed(v.yi, dec)}, ${fixed(v.yf, dec)}]`;
        const [cx, cy] = centre(v), mw = minWidthAt(cx, cy, p), hw = widthOf(state.home);
        $('zoom').textContent = `×${magnitude(hw / w)}`;
        const full = hw > mw ? Math.log(hw / w) / Math.log(hw / mw) : 1;
        const meter = $('meter');
        meter.style.width = `${(100 * Math.min(1, Math.max(0, full))).toFixed(1)}%`;
        meter.parentElement.classList.toggle('full', atDepthLimit());
        $('limit').textContent = `max ×${magnitude(hw / mw)}`;
        logCommand();
    }

    // the same view with the C++ program, for the browser console
    let lastCommand = '';
    function logCommand() {
        const own = state.palette.id !== 'viridis' && state.palette.id !== 'grey';
        const text = commandLine() + (own ? `\n# ${state.palette.name} is browser-only; the C++ program colours with Viridis, or grey with --bw` : '');
        if (text === lastCommand) return;
        lastCommand = text;
        console.info(text);
    }

    function showCursor(X, Y) {
        const p = state.params, v = state.view;
        const [x, y] = toComplex(X, Y, v, p);
        const dec = Math.min(17, Math.max(2, Math.ceil(-Math.log10(widthOf(v) / p.size)) + 1));
        $('cursor').textContent = `${fixed(x, dec)} ${y < 0 ? '−' : '+'} ${fixed(Math.abs(y), dec)}i`;
    }

    function commandLine() {
        const p = state.params, v = state.view, grey = state.palette.id === 'grey';
        const parts = ['./chm_fractal',
            '--xi', v.xi, '--xf', v.xf, '--yi', v.yi, '--yf', v.yf,
            '--kre', p.kre, '--kim', p.kim, '--m', p.m, '--map', p.steps === 2 ? 'q' : 'phi',
            '--iter', p.maxit, '--tol', p.tol, '--width', p.size, '--height', p.size,
            '--ss', p.ss, '--sample', p.endpoint ? 'endpoint' : 'center'];
        if (grey) parts.push('--bw');
        if (state.grid) parts.push('--grid');
        parts.push('--out', grey ? 'view.pgm' : 'view.ppm');
        return parts.map(String).join(' ');
    }

    function markDirty() {
        let dirty = true;
        try { dirty = pixelKey(readParams()) !== pixelKey(state.params); } catch (err) { dirty = true; }
        $('apply').classList.toggle('dirty', dirty);
        const m = Number(field('m').value);
        const nval = $('nval');
        nval.textContent = Number.isInteger(m) && m >= 1 ? `$N{=}${2 * m + 1}$` : '$N{=}?$';
        if (window.MathJax && window.MathJax.Hub) window.MathJax.Hub.Queue(['Typeset', window.MathJax.Hub, nval]);
    }

    // -------------------------------------------------------------- events --
    function imagePoint(e) {
        const rect = image.getBoundingClientRect(), W = state.params.size;
        const x = (e.clientX - rect.left) / rect.width * W;
        const y = (e.clientY - rect.top) / rect.height * W;
        return [Math.min(Math.max(x, 0), W), Math.min(Math.max(y, 0), W)];
    }

    function cancelDrag() {
        if (!state.drag) return;
        try { stage.releasePointerCapture(state.drag.id); } catch (err) { /* already released */ }
        state.drag = null;
        drawOverlay();
    }

    stage.addEventListener('pointerdown', (e) => {
        state.pointer = e.pointerType;
        if (e.button !== 0 || e.ctrlKey || !state.params) return;
        e.preventDefault();
        if (atDepthLimit()) {
            toast('Maximum zoom depth reached. Right-click to step out.', 'warn');
            return;
        }
        const [x, y] = imagePoint(e);
        state.drag = { id: e.pointerId, x0: x, y0: y, x1: x, y1: y, box: null };
        stage.setPointerCapture(e.pointerId);
    });

    stage.addEventListener('pointermove', (e) => {
        if (!state.params) return;
        const [x, y] = imagePoint(e);
        showCursor(x, y);
        const d = state.drag;
        if (!d || e.pointerId !== d.id) return;
        d.x1 = x;
        d.y1 = y;
        d.box = Math.max(Math.abs(x - d.x0), Math.abs(y - d.y0)) >= MIN_DRAG ? dragBox(d) : null;
        drawOverlay();
    });

    stage.addEventListener('pointerup', (e) => {
        const d = state.drag;
        if (!d || e.pointerId !== d.id) return;
        cancelDrag();
        if (!d.box) return;
        if (d.box.blocked) { toast('Maximum zoom depth reached.', 'warn'); return; }
        if (d.box.limited) toast('Zoom limited to the maximum depth.', 'warn');
        zoomTo(d.box);
    });

    stage.addEventListener('pointercancel', cancelDrag);
    stage.addEventListener('pointerleave', () => { if (!state.drag) $('cursor').textContent = '—'; });

    stage.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        if (state.drag) { cancelDrag(); return; }
        if (state.pointer === 'mouse') back();   // a long press on touch is not a step back
    });

    // A new kappa (or a frame that follows m) starts again from its own
    // window; other parameters re-render the current view.
    form.addEventListener('submit', (e) => {
        e.preventDefault();
        let p;
        try { p = readParams(); } catch (err) { toast(err.message, 'warn'); return; }
        setParams(p);
        const home = homeFor(p);
        if (!sameWindow(home, state.home)) {
            state.home = home;
            state.view = home;
            state.stack = [];
        }
        render();
        drawOverlay();
    });

    form.addEventListener('input', markDirty);

    $('grid').addEventListener('click', () => {
        state.grid = !state.grid;
        $('grid').setAttribute('aria-pressed', String(state.grid));
        drawOverlay();
        updateReadout();
    });

    $('reset').addEventListener('click', goHome);

    $('save').addEventListener('click', () => {
        const p = state.params, out = document.createElement('canvas');
        out.width = p.size;
        out.height = p.size;
        const g = out.getContext('2d');
        g.drawImage(image, 0, 0);
        if (state.grid) drawGrid(g, 1, Math.max(10, Math.round(p.size / 60)), Math.max(1, Math.round(p.size / 1000)));
        out.toBlob((blob) => {
            if (!blob) { toast('PNG export failed', 'warn'); return; }
            const a = document.createElement('a');
            a.href = URL.createObjectURL(blob);
            a.download = 'chm_fractal_view.png';
            document.body.appendChild(a);
            a.click();
            a.remove();
            setTimeout(() => URL.revokeObjectURL(a.href), 10000);
        }, 'image/png');
    });

    if (typeof ResizeObserver !== 'undefined') new ResizeObserver(drawOverlay).observe(stage);
    else window.addEventListener('resize', drawOverlay);

    // start on the main figure (--figure k0) unless the form says otherwise
    buildPalettes();
    let saved = null;
    try { saved = localStorage.getItem(PALETTE_KEY); } catch (err) { saved = null; }
    setPalette(saved || 'viridis');
    let p;
    try { p = readParams(); } catch (err) {
        form.reset();
        p = readParams();
    }
    setParams(p);
    state.home = homeFor(p);
    state.view = state.home;
    render();
    drawOverlay();
}

if (typeof document !== 'undefined') {
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initExplorer);
    else initExplorer();
}
