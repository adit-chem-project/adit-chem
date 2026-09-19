/* Trajectory player for the analysis page: fetches the frames, drives an ADIT_VIEWER and a small series plot. */
(function () {
  "use strict";
  var root = document.getElementById("adit-play");
  if (!root || !window.ADIT_FRAMES_URL || !window.ADIT_VIEWER) { return; }
  var LANG_EN = String(window.ADIT_LANG || "").indexOf("en") === 0;
  function T(ja, en) { return LANG_EN ? en : ja; }
  var canvas = document.getElementById("adit3d-play"), series = document.getElementById("adit-play-series");
  var slider = document.getElementById("adit-play-slider"), pos = document.getElementById("adit-play-pos");
  var note = document.getElementById("adit-play-note"), toggle = document.getElementById("adit-play-toggle");
  var btn = { first: document.getElementById("adit-play-first"), prev: document.getElementById("adit-play-prev"),
              next: document.getElementById("adit-play-next"), last: document.getElementById("adit-play-last") };
  var data = null, viewer = null, frame = 0, timer = null;

  function count() { return data ? data.n_frames : 0; }
  function label(i) {
    var t = data.times_fs ? T("、t = " + data.times_fs[i].toFixed(1) + " fs", ", t = " + data.times_fs[i].toFixed(1) + " fs") : "";
    return T("フレーム " + (i + 1) + " / " + count() + t, "frame " + (i + 1) + " / " + count() + t);
  }
  function xs() { return data.times_fs || data.frame_index; }
  function drawSeries() {
    if (!series || !data) { return; }
    var ratio = window.devicePixelRatio || 1;
    var w = series.clientWidth || 420, h = series.clientHeight || 180;
    series.width = Math.round(w * ratio); series.height = Math.round(h * ratio);
    var ctx = series.getContext("2d");
    var dark = window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
    ctx.fillStyle = dark ? "#1e1e1e" : "#ffffff"; ctx.fillRect(0, 0, series.width, series.height);
    var rows = [];
    if (data.energies_ev) { rows.push([T("エネルギー [eV]", "energy [eV]"), data.energies_ev]); }
    if (data.temperatures_k) { rows.push([T("温度 [K]", "temperature [K]"), data.temperatures_k]); }
    if (!rows.length) { series.style.display = "none"; return; }
    var x = xs(), xmin = Math.min.apply(null, x), xmax = Math.max.apply(null, x), xr = Math.max(xmax - xmin, 1e-9);
    var left = 64 * ratio, right = series.width - 10 * ratio, top = 8 * ratio, bottom = series.height - 22 * ratio;
    var rowH = (bottom - top) / rows.length;
    ctx.font = (10 * ratio) + "px system-ui, sans-serif";
    rows.forEach(function (row, r) {
      var y = row[1], ymin = Math.min.apply(null, y), ymax = Math.max.apply(null, y), yr = Math.max(ymax - ymin, 1e-9);
      var y0 = top + r * rowH, y1 = y0 + rowH - 6 * ratio;
      ctx.strokeStyle = dark ? "#555" : "#ccc"; ctx.lineWidth = 1; ctx.strokeRect(left, y0, right - left, y1 - y0);
      ctx.strokeStyle = dark ? "#ddd" : "#404048"; ctx.lineWidth = 1.2 * ratio; ctx.beginPath();
      for (var i = 0; i < y.length; i++) {
        var px = left + (x[i] - xmin) / xr * (right - left), py = y1 - (y[i] - ymin) / yr * (y1 - y0);
        if (i === 0) { ctx.moveTo(px, py); } else { ctx.lineTo(px, py); }
      }
      ctx.stroke();
      ctx.fillStyle = dark ? "#ccc" : "#404048"; ctx.textAlign = "right"; ctx.textBaseline = "middle";
      ctx.fillText(ymax.toPrecision(5), left - 4 * ratio, y0 + 6 * ratio);
      ctx.fillText(ymin.toPrecision(5), left - 4 * ratio, y1 - 6 * ratio);
      ctx.save(); ctx.translate(10 * ratio, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText(row[0], 0, 0); ctx.restore();
    });
    ctx.fillStyle = dark ? "#ccc" : "#404048"; ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(data.times_fs ? T("時刻 [fs]", "time [fs]") : (data.is_md ? T("フレーム", "frame") : T("ステップ", "step")), (left + right) / 2, bottom + 6 * ratio);
    ctx.textAlign = "left"; ctx.fillText(String(xmin), left, bottom + 2 * ratio); ctx.textAlign = "right"; ctx.fillText(String(xmax), right, bottom + 2 * ratio);
    var mx = left + (x[frame] - xmin) / xr * (right - left);
    ctx.strokeStyle = "#0A7AFF"; ctx.lineWidth = 1.6 * ratio; ctx.beginPath(); ctx.moveTo(mx, top); ctx.lineTo(mx, bottom); ctx.stroke();
    series._plot = { left: left, right: right, xmin: xmin, xr: xr };
  }
  function go(i) {
    var n = count();
    if (!n) { return; }
    frame = ((i % n) + n) % n;
    slider.value = frame;
    viewer.setPositions(data.positions[frame]);
    pos.textContent = label(frame);
    drawSeries();
  }
  function play() { if (timer || !count()) { return; } timer = setInterval(function () { go(frame + 1); }, 60); toggle.textContent = T("停止", "Pause"); }
  function pause() { if (timer) { clearInterval(timer); timer = null; } toggle.textContent = T("再生", "Play"); }

  fetch(window.ADIT_FRAMES_URL).then(function (r) { return r.json(); }).then(function (d) {
    data = d;
    note.textContent = (d.notes || []).join(" ");
    if (!d.n_frames) { canvas.style.display = "none"; if (series) { series.style.display = "none"; } slider.disabled = true; return; }
    viewer = window.ADIT_VIEWER(canvas, d.scene, { measure: false });
    slider.min = 0; slider.max = d.n_frames - 1; slider.value = 0;
    slider.addEventListener("input", function () { pause(); go(parseInt(slider.value, 10)); });
    toggle.addEventListener("click", function () { if (timer) { pause(); } else { play(); } });
    btn.first.addEventListener("click", function () { pause(); go(0); });
    btn.prev.addEventListener("click", function () { pause(); go(frame - 1); });
    btn.next.addEventListener("click", function () { pause(); go(frame + 1); });
    btn.last.addEventListener("click", function () { pause(); go(count() - 1); });
    if (series) {
      series.addEventListener("click", function (e) {
        var p = series._plot; if (!p) { return; }
        var rect = series.getBoundingClientRect();
        var xv = p.xmin + ((e.clientX - rect.left) * series.width / Math.max(1, rect.width) - p.left) / (p.right - p.left) * p.xr;
        var x = xs(), best = 0;
        for (var i = 1; i < x.length; i++) { if (Math.abs(x[i] - xv) < Math.abs(x[best] - xv)) { best = i; } }
        pause(); go(best);
      });
      window.addEventListener("resize", drawSeries);
    }
    go(0);
  }).catch(function (err) { note.textContent = T("フレームを読み込めません: ", "cannot load the frames: ") + err; });
})();
