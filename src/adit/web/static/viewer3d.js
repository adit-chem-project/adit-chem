/* 3D structure view for the browser UI. No external libraries.
   window.ADIT_VIEWER(canvas, scene, opts) builds a viewer; the page's #adit3d canvas is wired automatically. */
(function () {
  "use strict";
  var MAX_SELECTED = 4, CLICK_PX = 4, REBOND_ATOMS = 400, BOND_FACTOR = 1.2;
  var LANG_EN = String(window.ADIT_LANG || "").indexOf("en") === 0;
  function T(ja, en) { return LANG_EN ? en : ja; }

  function inv3(m) {
    var a = m[0][0], b = m[0][1], c = m[0][2], d = m[1][0], e = m[1][1], f = m[1][2], g = m[2][0], h = m[2][1], i = m[2][2];
    var det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if (Math.abs(det) < 1e-12) { return null; }
    return [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
            [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
            [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]];
  }
  function mulRow(v, m) {   /* row vector times matrix */
    return [v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0], v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1], v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2]];
  }
  function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }

  /* vector a -> b; the shortest periodic image when a cell (rows = lattice vectors) is given */
  function mic(a, b, cell, inv) {
    var d = sub(b, a);
    if (!cell || !inv) { return d; }
    var f = mulRow(d, inv);
    f = [f[0] - Math.round(f[0]), f[1] - Math.round(f[1]), f[2] - Math.round(f[2])];
    var best = null, bestLen = Infinity;
    for (var i = -1; i <= 1; i++) { for (var j = -1; j <= 1; j++) { for (var k = -1; k <= 1; k++) {
      var v = mulRow([f[0] + i, f[1] + j, f[2] + k], cell);
      var n = dot(v, v);
      if (n < bestLen) { bestLen = n; best = v; }
    } } }
    return best;
  }
  function distance(p, i, j, cell, inv) { return norm(mic(p[i], p[j], cell, inv)); }
  function angle(p, i, j, k, cell, inv) {
    var u = mic(p[j], p[i], cell, inv), v = mic(p[j], p[k], cell, inv);
    var c = dot(u, v) / (norm(u) * norm(v));
    return Math.acos(Math.min(1, Math.max(-1, c))) * 180 / Math.PI;
  }
  function dihedral(p, i, j, k, l, cell, inv) {
    var b1 = mic(p[i], p[j], cell, inv), b2 = mic(p[j], p[k], cell, inv), b3 = mic(p[k], p[l], cell, inv);
    var n1 = cross(b1, b2), n2 = cross(b2, b3), m = cross(n1, scale(b2, 1 / norm(b2)));
    return Math.atan2(dot(m, n2), dot(n1, n2)) * 180 / Math.PI;
  }
  function measure(p, sel, cell, inv) {
    if (sel.length === 2) { return { kind: T("距離", "distance"), value: distance(p, sel[0], sel[1], cell, inv), unit: "Å", digits: 3 }; }
    if (sel.length === 3) { return { kind: T("角度", "angle"), value: angle(p, sel[0], sel[1], sel[2], cell, inv), unit: "°", digits: 1 }; }
    if (sel.length === 4) { return { kind: T("二面角", "dihedral"), value: dihedral(p, sel[0], sel[1], sel[2], sel[3], cell, inv), unit: "°", digits: 1 }; }
    return null;
  }
  function bondsOf(pos, radii) {
    var out = [];
    for (var i = 0; i < pos.length - 1; i++) {
      for (var j = i + 1; j < pos.length; j++) {
        var d = norm(sub(pos[j], pos[i]));
        if (d > 0.1 && d <= (radii[i] + radii[j]) * BOND_FACTOR) { out.push([i, j]); }
      }
    }
    return out;
  }

  function makeViewer(canvas, scene, opts) {
    opts = opts || {};
    var ctx = canvas.getContext("2d");
    var state = { rx: -1.05, ry: 0.52, zoom: 1, panx: 0, pany: 0 };
    var selected = [];
    var cell = scene.cell || null, inv = cell ? inv3(cell) : null;
    var projected = [];

    function rotate(p) {
      var cy = Math.cos(state.ry), sy = Math.sin(state.ry);
      var x = p[0] * cy + p[2] * sy, z = -p[0] * sy + p[2] * cy;
      var cx = Math.cos(state.rx), sx = Math.sin(state.rx);
      return [x, p[1] * cx - z * sx, p[1] * sx + z * cx];
    }
    function project(p) {
      var r = rotate(p);
      var w = canvas.width, h = canvas.height;
      var s = (Math.min(w, h) * 0.42 / scene.scale) * state.zoom;
      return [w / 2 + r[0] * s + state.panx, h / 2 - r[1] * s + state.pany, r[2], s];
    }
    function measureText() {
      var m = measure(scene.positions, selected, cell, inv);
      if (!m) { return ""; }
      var names = selected.map(function (i) { return scene.symbols[i] + (i + 1); });
      return m.kind + " " + names.join("-") + ": " + m.value.toFixed(m.digits) + " " + m.unit;
    }
    function label(x, y, text, color, ratio, dark) {
      ctx.font = "bold " + (11 * ratio) + "px system-ui, sans-serif";
      var w = ctx.measureText(text).width + 8 * ratio, h = 16 * ratio;
      x = Math.min(Math.max(2, x), Math.max(2, canvas.width - w - 2)); y = Math.min(Math.max(2, y), Math.max(2, canvas.height - h - 2));
      ctx.fillStyle = dark ? "rgba(30,30,30,0.8)" : "rgba(255,255,255,0.86)";
      ctx.fillRect(x, y, w, h);
      ctx.fillStyle = color; ctx.textBaseline = "middle"; ctx.textAlign = "left";
      ctx.fillText(text, x + 4 * ratio, y + h / 2);
      ctx.textBaseline = "alphabetic";
    }

    function draw() {
      var w = canvas.width, h = canvas.height;
      var dark = window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
      var ratio = window.devicePixelRatio || 1;
      ctx.fillStyle = dark ? "#1e1e1e" : "#ffffff";
      ctx.fillRect(0, 0, w, h);
      var i, a, b;
      ctx.strokeStyle = dark ? "#7a7a7a" : "#9a9a9a";
      ctx.lineWidth = 1;
      for (i = 0; i < scene.cell_lines.length; i++) {
        a = project(scene.cell_lines[i][0]); b = project(scene.cell_lines[i][1]);
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
      }
      for (i = 0; i < scene.bonds.length; i++) {
        var ia = scene.bonds[i][0], ib = scene.bonds[i][1];
        a = project(scene.positions[ia]); b = project(scene.positions[ib]);
        var mx = (a[0] + b[0]) / 2, my = (a[1] + b[1]) / 2;
        var width = Math.max(1.5, 0.16 * a[3]);
        ctx.lineWidth = width + 2;
        ctx.strokeStyle = dark ? "rgba(255,255,255,0.30)" : "rgba(0,0,0,0.35)";
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
        ctx.lineWidth = width;
        ctx.strokeStyle = scene.colors[ia];
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(mx, my); ctx.stroke();
        ctx.strokeStyle = scene.colors[ib];
        ctx.beginPath(); ctx.moveTo(mx, my); ctx.lineTo(b[0], b[1]); ctx.stroke();
      }
      projected = scene.positions.map(function (p, k) { var q = project(p); return { x: q[0], y: q[1], z: q[2], r: Math.max(2, scene.radii[k] * 0.55 * q[3]) }; });
      var order = projected.map(function (q, k) { return k; });
      order.sort(function (u, v) { return projected[u].z - projected[v].z; });
      var accent = dark ? "#3A95FF" : "#0A7AFF";
      for (i = 0; i < order.length; i++) {
        var k = order[i], q = projected[k], rad = q.r;
        var g = ctx.createRadialGradient(q.x - rad * 0.35, q.y - rad * 0.35, rad * 0.1, q.x, q.y, rad);
        g.addColorStop(0, "#ffffff"); g.addColorStop(0.35, scene.colors[k]); g.addColorStop(1, scene.colors[k]);
        ctx.fillStyle = g;
        ctx.beginPath(); ctx.arc(q.x, q.y, rad, 0, 2 * Math.PI); ctx.fill();
        ctx.strokeStyle = dark ? "rgba(0,0,0,0.55)" : "rgba(0,0,0,0.35)";
        ctx.lineWidth = 1; ctx.stroke();
        if (selected.indexOf(k) >= 0) {
          ctx.strokeStyle = accent; ctx.lineWidth = 2.5 * ratio;
          ctx.beginPath(); ctx.arc(q.x, q.y, rad + 3 * ratio, 0, 2 * Math.PI); ctx.stroke();
        }
      }
      if (selected.length) {
        ctx.strokeStyle = accent; ctx.lineWidth = 2 * ratio; ctx.setLineDash([6 * ratio, 4 * ratio]);
        for (i = 0; i + 1 < selected.length; i++) {
          a = projected[selected[i]]; b = projected[selected[i + 1]];
          ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
        }
        ctx.setLineDash([]);
        for (i = 0; i < selected.length; i++) {
          a = projected[selected[i]];
          label(a.x + a.r + 4 * ratio, a.y - a.r - 16 * ratio, String(i + 1), accent, ratio, dark);
        }
        var text = measureText();
        if (text) {
          var at;
          if (selected.length === 2) { at = [(projected[selected[0]].x + projected[selected[1]].x) / 2, (projected[selected[0]].y + projected[selected[1]].y) / 2]; }
          else if (selected.length === 3) { at = [projected[selected[1]].x, projected[selected[1]].y + 18 * ratio]; }
          else { at = [(projected[selected[1]].x + projected[selected[2]].x) / 2, (projected[selected[1]].y + projected[selected[2]].y) / 2]; }
          label(at[0] + 6 * ratio, at[1] + 6 * ratio, text, accent, ratio, dark);
        }
      }
      drawAxes(dark);
      if (opts.onSelection) { opts.onSelection(selected.slice(), measureText()); }
    }

    function drawAxes(dark) {
      var ratio = window.devicePixelRatio || 1;
      var arm = 22 * ratio, ox = 14 * ratio + arm, oy = canvas.height - 14 * ratio - arm;
      var axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      var colors = dark ? ["#E8756F", "#7BD07B", "#6FA8E8"] : ["#D9534F", "#5CB85C", "#4A90D9"];
      var names = ["x", "y", "z"];
      ctx.font = (11 * ratio) + "px system-ui, sans-serif";
      ctx.textAlign = "center"; ctx.textBaseline = "middle";
      for (var k = 0; k < 3; k++) {
        var v = rotate(axes[k]);
        var ex = ox + v[0] * arm, ey = oy - v[1] * arm;
        ctx.strokeStyle = colors[k]; ctx.fillStyle = colors[k];
        ctx.lineWidth = 1.6 * ratio; ctx.lineCap = "round";
        ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ex, ey); ctx.stroke();
        ctx.beginPath(); ctx.arc(ex, ey, 2.2 * ratio, 0, 2 * Math.PI); ctx.fill();
        ctx.fillText(names[k], ox + v[0] * (arm + 9 * ratio), oy - v[1] * (arm + 9 * ratio));
      }
      ctx.textAlign = "start"; ctx.textBaseline = "alphabetic";
    }

    function resize() {
      var ratio = window.devicePixelRatio || 1;
      var width = canvas.clientWidth || 420, height = canvas.clientHeight || 320;
      canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
      draw();
    }

    function canvasPoint(e) {
      var rect = canvas.getBoundingClientRect();
      return [(e.clientX - rect.left) * canvas.width / Math.max(1, rect.width), (e.clientY - rect.top) * canvas.height / Math.max(1, rect.height)];
    }
    /* front-most atom drawn under the point, or -1 */
    function pick(x, y) {
      var best = -1, bestZ = -Infinity;
      for (var k = 0; k < projected.length; k++) {
        var q = projected[k];
        if (Math.hypot(q.x - x, q.y - y) <= q.r + 1.5 && q.z > bestZ) { bestZ = q.z; best = k; }
      }
      return best;
    }
    function setSelection(list) {
      selected = [];
      for (var i = 0; i < list.length; i++) {
        var k = list[i];
        if (k >= 0 && k < scene.positions.length && selected.indexOf(k) < 0) { selected.push(k); }
      }
      selected = selected.slice(-MAX_SELECTED);
      draw();
    }
    function click(x, y, additive) {
      if (!opts.measure) { return; }
      var k = pick(x, y);
      if (k < 0) { if (!additive) { setSelection([]); } return; }
      if (additive) {
        var at = selected.indexOf(k);
        if (at >= 0) { var copy = selected.slice(); copy.splice(at, 1); setSelection(copy); }
        else { setSelection(selected.concat([k])); }
      } else { setSelection([k]); }
    }

    var last = null, button = 0, press = null, dragged = false;
    canvas.addEventListener("mousedown", function (e) { last = [e.clientX, e.clientY]; press = last; dragged = false; button = e.button; e.preventDefault(); });
    window.addEventListener("mouseup", function () { last = null; });
    canvas.addEventListener("mouseup", function (e) {
      if (press && !dragged && e.button === 0) { var p = canvasPoint(e); click(p[0], p[1], e.shiftKey); }
      press = null;
    });
    window.addEventListener("mousemove", function (e) {
      if (!last) { return; }
      var dx = e.clientX - last[0], dy = e.clientY - last[1];
      last = [e.clientX, e.clientY];
      if (press && Math.hypot(e.clientX - press[0], e.clientY - press[1]) > CLICK_PX) { dragged = true; }
      if (button === 2) { state.panx += dx; state.pany += dy; }
      else { state.ry += dx * 0.01; state.rx += dy * 0.01; }
      draw();
    });
    canvas.addEventListener("contextmenu", function (e) { e.preventDefault(); });
    canvas.addEventListener("wheel", function (e) {
      e.preventDefault();
      state.zoom *= (e.deltaY < 0 ? 1.1 : 1 / 1.1);
      state.zoom = Math.min(20, Math.max(0.2, state.zoom));
      draw();
    }, { passive: false });
    canvas.addEventListener("dblclick", function () {
      state.rx = -1.05; state.ry = 0.52; state.zoom = 1; state.panx = 0; state.pany = 0; draw();
    });
    window.addEventListener("resize", resize);
    resize();

    return {
      draw: draw, resize: resize, state: state, scene: scene,
      setView: function (rx, ry) { state.rx = rx; state.ry = ry; state.zoom = 1; state.panx = 0; state.pany = 0; draw(); },
      setPositions: function (pos) {
        scene.positions = pos;
        if (pos.length <= REBOND_ATOMS) { scene.bonds = bondsOf(pos, scene.radii); }
        draw();
      },
      selected: function () { return selected.slice(); },
      setSelection: setSelection,
      clearSelection: function () { setSelection([]); },
      measureText: measureText,
      measure: function () { return measure(scene.positions, selected, cell, inv); },
      pick: pick
    };
  }
  window.ADIT_VIEWER = makeViewer;
  window.ADIT_GEOMETRY = { mic: mic, distance: distance, angle: angle, dihedral: dihedral, inv3: inv3 };

  var canvas = document.getElementById("adit3d");
  if (!canvas || !window.ADIT_SCENE) { return; }
  var measureDiv = document.getElementById("adit3d-measure"), clearBtn = document.getElementById("adit3d-clear");
  var sendBtn = document.getElementById("adit3d-send-fixed"), indicesEl = document.getElementById("adit3d-indices");
  var hintText = measureDiv ? measureDiv.textContent : "";
  var viewer = makeViewer(canvas, window.ADIT_SCENE, {
    measure: true,
    onSelection: function (sel, text) {
      var ones = sel.map(function (i) { return i + 1; });
      if (measureDiv) {
        measureDiv.textContent = text || (sel.length === 1 ? T("選択: " + window.ADIT_SCENE.symbols[sel[0]] + ones[0] + " (もう 1 個を Shift+クリックすると距離)",
                                                                  "selected: " + window.ADIT_SCENE.symbols[sel[0]] + ones[0] + " (Shift+click another atom for the distance)") : hintText);
      }
      if (clearBtn) { clearBtn.disabled = sel.length === 0; }
      if (sendBtn) { sendBtn.disabled = sel.length === 0; }
      if (indicesEl) { indicesEl.textContent = sel.length ? ones.join(",") : "-"; }
    }
  });
  var views = document.querySelectorAll("[data-view]");
  for (var v = 0; v < views.length; v++) {
    views[v].addEventListener("click", function (e) {
      var parts = e.currentTarget.getAttribute("data-view").split(",");
      viewer.setView(parseFloat(parts[0]), parseFloat(parts[1]));
    });
  }
  if (clearBtn) { clearBtn.addEventListener("click", function () { viewer.clearSelection(); }); }
  if (sendBtn) {
    sendBtn.addEventListener("click", function () {
      var field = document.getElementById("fixed");
      var sel = viewer.selected();
      if (!field || field.disabled || !sel.length) { return; }
      var have = field.value.replace(/\s+/g, "").split(",").filter(function (s) { return s; });
      var plain = {}, axes = [];
      have.forEach(function (part) {
        if (part.indexOf(":") >= 0) { axes.push(part); return; }
        var m = part.match(/^(\d+)(?:-(\d+))?$/);
        if (!m) { axes.push(part); return; }
        var lo = parseInt(m[1], 10), hi = m[2] ? parseInt(m[2], 10) : lo;
        for (var k = lo; k <= hi; k++) { plain[k] = true; }
      });
      sel.forEach(function (i) { plain[i + 1] = true; });
      var idx = Object.keys(plain).map(Number).sort(function (a, b) { return a - b; });
      var runs = [];
      idx.forEach(function (i) { if (runs.length && runs[runs.length - 1][1] === i - 1) { runs[runs.length - 1][1] = i; } else { runs.push([i, i]); } });
      var text = runs.map(function (r) { return r[0] === r[1] ? String(r[0]) : r[0] + "-" + r[1]; });
      field.value = text.concat(axes).join(",");
      if (measureDiv) { measureDiv.textContent = T("原子 " + sel.map(function (i) { return i + 1; }).join(",") + " を「固定原子」の欄に入れました (プレビューで反映)",
                                                   "atoms " + sel.map(function (i) { return i + 1; }).join(",") + " were put into the fixed-atoms field (preview to apply)"); }
    });
  }
})();
