/* Isosurface box of the analysis page: fetches stats / triangles from /isosurface.json and draws them with ADIT_VIEWER. */
(function () {
  "use strict";
  var root = document.getElementById("adit-iso");
  if (!root || !window.ADIT_ISO_URL || !window.ADIT_VIEWER) { return; }
  var LANG_EN = String(window.ADIT_LANG || "").indexOf("en") === 0;
  function T(ja, en) { return LANG_EN ? en : ja; }
  var file = document.getElementById("adit-iso-file"), level = document.getElementById("adit-iso-level");
  var range = document.getElementById("adit-iso-range"), note = document.getElementById("adit-iso-note");
  var pos = document.getElementById("adit-iso-pos"), neg = document.getElementById("adit-iso-neg");
  var opacity = document.getElementById("adit-iso-opacity"), stride = document.getElementById("adit-iso-stride");
  var show = document.getElementById("adit-iso-show"), clear = document.getElementById("adit-iso-clear"), useStride = document.getElementById("adit-iso-use-stride");
  var canvas = document.getElementById("adit3d-iso");
  var viewer = null, data = null, suggested = null;

  function rgb(hex) { return [parseInt(hex.substr(1, 2), 16), parseInt(hex.substr(3, 2), 16), parseInt(hex.substr(5, 2), 16)]; }
  function fmt(x) { return Number(x).toPrecision(4); }
  function url(withLevel) {
    var u = window.ADIT_ISO_URL + "&file=" + encodeURIComponent(file.value);
    if (withLevel) { u += "&level=" + encodeURIComponent(level.value.trim()) + "&stride=" + encodeURIComponent(stride.value || "1"); }
    return u;
  }
  function shifted(tri, center) {
    var out = new Array(tri.length);
    for (var i = 0; i < tri.length; i++) { out[i] = tri[i] - center[i % 3]; }
    return out;
  }
  function apply() {
    if (!viewer || !data) { return; }
    var list = [];
    for (var i = 0; i < data.surfaces.length; i++) {
      var s = data.surfaces[i];
      list.push({ tri: shifted(s.tri, data.center), normal: s.normal, rgb: rgb(s.level >= 0 ? pos.value : neg.value), opacity: Number(opacity.value) / 100 });
    }
    viewer.setSurfaces(list);
  }
  function showStats(st) {
    var unit = st.unit ? " " + st.unit : "";
    var shape = st.shape.join("×");
    range.textContent = T("値の範囲: 最小 " + fmt(st.min) + unit + "、最大 " + fmt(st.max) + unit + " (格子 " + shape + "、" + st.n_points.toLocaleString() + " 点" + (st.periodic ? "、周期的" : "") + ")。参考: |最大| の 10 % = " + fmt(st.tenth_of_abs_max) + unit + " (この 10 % は単に |最大| の 1/10 で、化学的な根拠はありません。等値は利用者が決めます)",
                          "value range: min " + fmt(st.min) + unit + ", max " + fmt(st.max) + unit + " (grid " + shape + ", " + st.n_points.toLocaleString() + " points" + (st.periodic ? ", periodic" : "") + "). For reference: 10 % of |max| = " + fmt(st.tenth_of_abs_max) + unit + " (that 10 % is just one tenth of |max|, with no chemical basis; the level is your choice)");
  }
  function fetchJSON(u, done) {
    var x = new XMLHttpRequest();
    x.open("GET", u, true);
    x.onreadystatechange = function () {
      if (x.readyState !== 4) { return; }
      if (x.status !== 200) { note.textContent = T("読めません: ", "cannot read: ") + x.responseText; return; }
      try { done(JSON.parse(x.responseText)); } catch (e) { note.textContent = T("応答を読めません", "cannot parse the response"); }
    };
    x.send();
  }
  function build(d) {
    data = d;
    showStats(d.stats);
    if (!viewer || viewer.scene !== d.scene) { viewer = window.ADIT_VIEWER(canvas, d.scene, { measure: true }); }
  }
  function loadStats() {
    note.textContent = ""; useStride.hidden = true; suggested = null;
    fetchJSON(url(false), function (d) { d.surfaces = []; build(d); viewer.clearSurfaces(); });
  }
  function draw() {
    if (!level.value.trim()) { note.textContent = T("等値を入れてください (上の範囲を見て決めます。既定値はありません)", "enter a level (choose it from the range above; there is no default)"); return; }
    note.textContent = T("計算しています…", "computing…"); useStride.hidden = true; suggested = null;
    fetchJSON(url(true), function (d) {
      build(d);
      apply();
      var parts = d.surfaces.map(function (s) { return (s.level >= 0 ? "+" : "") + fmt(s.level) + ": " + s.n_triangles.toLocaleString() + T(" 枚", " triangles"); });
      var text = d.surfaces.length ? T("三角形 ", "Triangles ") + parts.join(" / ") + "." : "";
      if (d.notes && d.notes.length) { text += (text ? " " : "") + d.notes.join(" "); }
      note.textContent = text;
      if (d.suggested_stride) {
        suggested = d.suggested_stride;
        useStride.textContent = T("間引きを " + suggested + " にする", "Set the stride to " + suggested); useStride.hidden = false;
      }
    });
  }
  file.addEventListener("change", loadStats);
  show.addEventListener("click", draw);
  clear.addEventListener("click", function () { if (viewer) { viewer.clearSurfaces(); } if (data) { data.surfaces = []; } note.textContent = ""; });
  useStride.addEventListener("click", function () { if (suggested) { stride.value = suggested; draw(); } });
  pos.addEventListener("input", apply); neg.addEventListener("input", apply); opacity.addEventListener("input", apply);
  loadStats();
})();
