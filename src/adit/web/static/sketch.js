/* 2D molecule sketcher for the browser UI. No external libraries; SMILES conversion happens on the server. */
(function () {
  "use strict";
  var canvas = document.getElementById("sketchpad");
  if (!canvas) { return; }
  var ctx = canvas.getContext("2d");
  var field = document.getElementById("sketch_data");
  var BOND = 44, HIT = 14, CLICK = 10;
  var state = { atoms: [], bonds: [], elem: "C", tool: "1", selected: -1 };
  var history = [], future = [];
  var colors = window.ADIT_SKETCH_COLORS || {};

  try { state.atoms = (JSON.parse(field.value || "{}").atoms) || []; state.bonds = (JSON.parse(field.value || "{}").bonds) || []; } catch (e) {}

  function snapshot() { history.push(JSON.stringify({ atoms: state.atoms, bonds: state.bonds })); if (history.length > 200) history.shift(); future = []; }
  function restore(text) { var d = JSON.parse(text); state.atoms = d.atoms; state.bonds = d.bonds; state.selected = -1; }
  function save() { field.value = JSON.stringify({ atoms: state.atoms, bonds: state.bonds }); }

  function atomAt(x, y) {
    for (var i = state.atoms.length - 1; i >= 0; i--) {
      var a = state.atoms[i];
      if ((a.x - x) * (a.x - x) + (a.y - y) * (a.y - y) <= HIT * HIT) { return i; }
    }
    return -1;
  }
  function bondAt(x, y) {
    for (var i = 0; i < state.bonds.length; i++) {
      var a = state.atoms[state.bonds[i].a], b = state.atoms[state.bonds[i].b];
      var dx = b.x - a.x, dy = b.y - a.y, len2 = dx * dx + dy * dy;
      if (!len2) { continue; }
      var t = Math.max(0, Math.min(1, ((x - a.x) * dx + (y - a.y) * dy) / len2));
      var px = a.x + t * dx - x, py = a.y + t * dy - y;
      if (px * px + py * py <= (HIT * 0.7) * (HIT * 0.7)) { return i; }
    }
    return -1;
  }
  function bondBetween(i, j) {
    for (var k = 0; k < state.bonds.length; k++) {
      var b = state.bonds[k];
      if ((b.a === i && b.b === j) || (b.a === j && b.b === i)) { return k; }
    }
    return -1;
  }
  function removeAtom(i) {
    state.bonds = state.bonds.filter(function (b) { return b.a !== i && b.b !== i; });
    state.bonds.forEach(function (b) { if (b.a > i) b.a--; if (b.b > i) b.b--; });
    state.atoms.splice(i, 1);
    state.selected = -1;
  }
  function placeRing(name, cx, cy) {
    var n = name === "cyclopentane" ? 5 : 6, aromatic = name === "benzene";
    var r = BOND / (2 * Math.sin(Math.PI / n)), start = state.atoms.length, k;
    for (k = 0; k < n; k++) {
      var t = 2 * Math.PI * k / n - Math.PI / 2;
      state.atoms.push({ x: cx + r * Math.cos(t), y: cy + r * Math.sin(t), elem: "C", charge: 0 });
    }
    for (k = 0; k < n; k++) {
      state.bonds.push({ a: start + k, b: start + (k + 1) % n, order: (aromatic && k % 2 === 0) ? 2 : 1 });
    }
  }

  function draw() {
    var dark = window.ADIT_DARK();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = dark ? "#2b2d35" : "#ffffff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    var ratio = canvas.width / canvas.clientWidth;
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    var fg = dark ? "#e6e6e6" : "#202020";
    var i;
    for (i = 0; i < state.bonds.length; i++) {
      var bo = state.bonds[i], a = state.atoms[bo.a], b = state.atoms[bo.b];
      var dx = b.x - a.x, dy = b.y - a.y, len = Math.hypot(dx, dy) || 1;
      var nx = -dy / len * 3.5, ny = dx / len * 3.5;
      ctx.strokeStyle = fg; ctx.lineWidth = 2;
      var offs = bo.order === 1 ? [0] : bo.order === 2 ? [-1, 1] : [-1, 0, 1];
      for (var k = 0; k < offs.length; k++) {
        ctx.beginPath();
        ctx.moveTo(a.x + nx * offs[k], a.y + ny * offs[k]);
        ctx.lineTo(b.x + nx * offs[k], b.y + ny * offs[k]);
        ctx.stroke();
      }
    }
    ctx.font = "13px system-ui, sans-serif";
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    for (i = 0; i < state.atoms.length; i++) {
      var at = state.atoms[i];
      if (i === state.selected) {
        ctx.strokeStyle = "#2f7ae5"; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(at.x, at.y, 12, 0, 2 * Math.PI); ctx.stroke();
      }
      if (at.elem === "C" && at.charge === 0) { continue; }
      ctx.fillStyle = dark ? "#2b2d35" : "#ffffff";
      ctx.beginPath(); ctx.arc(at.x, at.y, 10, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = colors[at.elem] || fg;
      var text = at.elem + (at.charge > 0 ? "+".repeat(at.charge) : at.charge < 0 ? "−".repeat(-at.charge) : "");
      ctx.fillText(text, at.x, at.y);
    }
    save();
  }

  function pos(e) {
    var r = canvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  var down = null, downAtom = -1;
  canvas.addEventListener("mousedown", function (e) {
    var p = pos(e);
    if (e.button === 2) { return; }
    down = p; downAtom = atomAt(p[0], p[1]);
  });
  canvas.addEventListener("contextmenu", function (e) {
    e.preventDefault();
    var p = pos(e), i = atomAt(p[0], p[1]);
    snapshot();
    if (i >= 0) { removeAtom(i); } else { var k = bondAt(p[0], p[1]); if (k >= 0) { state.bonds.splice(k, 1); } }
    draw();
  });
  canvas.addEventListener("mouseup", function (e) {
    if (!down) { return; }
    var p = pos(e), moved = Math.hypot(p[0] - down[0], p[1] - down[1]) > CLICK;
    snapshot();
    if (state.tool === "eraser") {
      var i = atomAt(p[0], p[1]);
      if (i >= 0) { removeAtom(i); } else { var k = bondAt(p[0], p[1]); if (k >= 0) { state.bonds.splice(k, 1); } }
    } else if (state.tool.indexOf("ring:") === 0) {
      placeRing(state.tool.slice(5), p[0], p[1]);
    } else if (downAtom >= 0 && moved) {
      var target = atomAt(p[0], p[1]);
      var order = parseInt(state.tool, 10);
      if (target >= 0 && target !== downAtom) {
        var e0 = bondBetween(downAtom, target);
        if (e0 >= 0) { state.bonds[e0].order = order; } else { state.bonds.push({ a: downAtom, b: target, order: order }); }
      } else {
        var a0 = state.atoms[downAtom];
        var ang = Math.atan2(p[1] - a0.y, p[0] - a0.x);
        ang = Math.round(ang / (Math.PI / 6)) * (Math.PI / 6);
        state.atoms.push({ x: a0.x + BOND * Math.cos(ang), y: a0.y + BOND * Math.sin(ang), elem: state.elem, charge: 0 });
        state.bonds.push({ a: downAtom, b: state.atoms.length - 1, order: order });
        state.selected = state.atoms.length - 1;
      }
    } else if (downAtom >= 0) {
      state.atoms[downAtom].elem = state.elem;
      state.selected = downAtom;
    } else {
      var kb = bondAt(p[0], p[1]);
      if (kb >= 0) {
        var want = parseInt(state.tool, 10);
        state.bonds[kb].order = state.bonds[kb].order === want ? (want % 3) + 1 : want;
      } else {
        state.atoms.push({ x: p[0], y: p[1], elem: state.elem, charge: 0 });
        state.selected = state.atoms.length - 1;
      }
    }
    down = null; downAtom = -1;
    draw();
  });

  function pick(group, value, key) {
    var nodes = document.querySelectorAll("[data-" + group + "]");
    for (var i = 0; i < nodes.length; i++) {
      nodes[i].classList.toggle("active", nodes[i].getAttribute("data-" + group) === value);
    }
    state[key] = value;
  }
  document.querySelectorAll("[data-elem]").forEach(function (b) {
    b.addEventListener("click", function () { pick("elem", b.getAttribute("data-elem"), "elem"); });
  });
  document.querySelectorAll("[data-tool]").forEach(function (b) {
    b.addEventListener("click", function () { pick("tool", b.getAttribute("data-tool"), "tool"); });
  });
  document.querySelectorAll("[data-charge]").forEach(function (b) {
    b.addEventListener("click", function () {
      if (state.selected < 0) { return; }
      snapshot();
      state.atoms[state.selected].charge += parseInt(b.getAttribute("data-charge"), 10);
      draw();
    });
  });
  var clear = document.getElementById("sketch_clear");
  if (clear) { clear.addEventListener("click", function () { snapshot(); state.atoms = []; state.bonds = []; state.selected = -1; draw(); }); }
  var undo = document.getElementById("sketch_undo");
  if (undo) {
    undo.addEventListener("click", function () {
      if (!history.length) { return; }
      future.push(JSON.stringify({ atoms: state.atoms, bonds: state.bonds }));
      restore(history.pop()); draw();
    });
  }
  var redo = document.getElementById("sketch_redo");
  if (redo) {
    redo.addEventListener("click", function () {
      if (!future.length) { return; }
      history.push(JSON.stringify({ atoms: state.atoms, bonds: state.bonds }));
      restore(future.pop()); draw();
    });
  }
  window.addEventListener("keydown", function (e) {
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z") { e.preventDefault(); (e.shiftKey ? redo : undo).click(); }
    if (e.key === "Delete" && state.selected >= 0) { snapshot(); removeAtom(state.selected); draw(); }
  });

  function resize() {
    var ratio = window.devicePixelRatio || 1;
    canvas.width = Math.round(canvas.clientWidth * ratio);
    canvas.height = Math.round(canvas.clientHeight * ratio);
    draw();
  }
  window.addEventListener("resize", resize);
  pick("elem", "C", "elem"); pick("tool", "1", "tool");
  resize();
})();
