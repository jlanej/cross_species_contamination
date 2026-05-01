"""Inline JavaScript + CSS that powers the cohort report's interactivity.

The report deliberately ships as a *single* HTML file with no
external assets, so this module just exposes two string constants
(``COHORT_CSS`` and ``COHORT_JS``) that the renderer concatenates
into ``<style>`` / ``<script>`` blocks.

The vanilla-JS module here is intentionally tiny (~3 KB minified)
and degrades gracefully to a vertically-scrolling table when JS is
disabled or fails to evaluate.  Tables opt in by adding the class
``paginated`` plus optional ``data-page-size`` and ``data-table-id``
attributes.

Features:

* Pagination with first/prev/next/last + page indicator.
* Click-to-sort column headers (numeric vs string detected from the
  cell content).
* A free-text filter input that hides non-matching rows.
* Domain / detection-band filter dropdowns when the relevant
  ``<th data-filter="…">`` markers are present.
"""

from __future__ import annotations


COHORT_CSS = r"""
.legend { font-size: 11.5px; color: #444; margin: 0.3em 0 0.6em;
          display: flex; flex-wrap: wrap; gap: 0.4em 0.9em; }
.legend .lg-item { display: inline-flex; align-items: center; gap: 4px; }
.legend .lg-sw { display: inline-block; width: 10px; height: 10px;
                  border: 1px solid #444; }
table.paginated { font-variant-numeric: tabular-nums; }
table.paginated thead th { cursor: pointer; user-select: none;
                            position: sticky; top: 0; }
table.paginated thead th[data-sort]::after { content: " \2195"; opacity: 0.4; }
table.paginated thead th[data-sort="asc"]::after { content: " \25B2"; opacity: 1; }
table.paginated thead th[data-sort="desc"]::after { content: " \25BC"; opacity: 1; }
.tbl-controls { display: flex; flex-wrap: wrap; gap: 0.7em;
                 align-items: center; margin: 0.3em 0; font-size: 12.5px; }
.tbl-controls input[type=search],
.tbl-controls select { font-size: 12px; padding: 2px 6px; }
.tbl-controls .pager button { font-size: 12px; padding: 2px 8px;
                                margin-right: 2px; cursor: pointer; }
.tbl-controls .pager button[disabled] { opacity: 0.4; cursor: default; }
.tbl-controls .pager-info { color: #555; }
.quadrant-grid { display: grid; grid-template-columns: 1fr 1fr;
                  gap: 0.6em; margin: 0.6em 0; }
.quadrant-grid .quadrant { background: #fafafa; padding: 0.6em;
                            border: 1px solid #ddd; border-radius: 4px; }
.quadrant-grid .quadrant h4 { margin: 0 0 0.3em; font-size: 13px;
                               text-transform: uppercase; letter-spacing: 0.04em; }
.partition-grid { display: grid; grid-template-columns: repeat(3, 1fr);
                   gap: 0.6em; }
.partition-grid > div { border: 1px solid #ddd; border-radius: 4px;
                         padding: 0.5em 0.7em; background: #fafafa; }
.partition-grid h4 { margin: 0 0 0.4em; font-size: 13px; }
.partition-grid table { font-size: 12px; }
details.species-details { margin: 0.4em 0; padding: 0.3em 0.6em;
                           background: #fafafa; border-left: 3px solid #56B4E9; }
details.species-details summary { cursor: pointer; font-weight: 600; }
.dl-link { font-size: 12px; }
"""


COHORT_JS = r"""
(function () {
  function parseNum(s) {
    if (s == null) return NaN;
    s = String(s).trim().replace(/,/g, '').replace(/%$/, '');
    if (s === '' || s === 'NA' || s === 'N/A') return NaN;
    var n = parseFloat(s);
    return isNaN(n) ? NaN : n;
  }
  function cellSortValue(td) {
    var v = td.dataset.sort != null ? td.dataset.sort : td.textContent;
    var n = parseNum(v);
    return isNaN(n) ? String(v).toLowerCase() : n;
  }
  function compare(a, b) {
    var an = typeof a === 'number';
    var bn = typeof b === 'number';
    if (an && bn) return a - b;
    if (an) return -1;
    if (bn) return 1;
    return a < b ? -1 : a > b ? 1 : 0;
  }
  function setupTable(tbl) {
    var pageSize = parseInt(tbl.dataset.pageSize, 10) || 25;
    var thead = tbl.querySelector('thead');
    var tbody = tbl.querySelector('tbody');
    if (!thead || !tbody) return;
    var rows = Array.prototype.slice.call(tbody.querySelectorAll('tr'));
    if (rows.length === 0) return;

    var visibleRows = rows.slice();
    var page = 0;

    var controls = document.createElement('div');
    controls.className = 'tbl-controls';

    var search = document.createElement('input');
    search.type = 'search';
    search.placeholder = 'filter…';
    search.style.minWidth = '160px';
    controls.appendChild(search);

    // Per-column dropdown filters (attribute data-filter on <th>).
    var headers = thead.querySelectorAll('th');
    var filterSelects = [];
    headers.forEach(function (th, idx) {
      if (!th.dataset.filter) return;
      var sel = document.createElement('select');
      var options = new Set();
      rows.forEach(function (r) {
        var c = r.children[idx];
        if (c) options.add(c.textContent.trim());
      });
      var opt0 = document.createElement('option');
      opt0.value = '';
      opt0.textContent = 'all ' + th.dataset.filter;
      sel.appendChild(opt0);
      Array.from(options).sort().forEach(function (v) {
        var o = document.createElement('option');
        o.value = v;
        o.textContent = v;
        sel.appendChild(o);
      });
      sel.addEventListener('change', function () {
        page = 0;
        applyFilters();
      });
      controls.appendChild(sel);
      filterSelects.push({ idx: idx, sel: sel });
    });

    var pager = document.createElement('span');
    pager.className = 'pager';
    var btnFirst = document.createElement('button'); btnFirst.textContent = '«';
    var btnPrev = document.createElement('button'); btnPrev.textContent = '‹';
    var btnNext = document.createElement('button'); btnNext.textContent = '›';
    var btnLast = document.createElement('button'); btnLast.textContent = '»';
    var info = document.createElement('span');
    info.className = 'pager-info';
    pager.appendChild(btnFirst);
    pager.appendChild(btnPrev);
    pager.appendChild(info);
    pager.appendChild(btnNext);
    pager.appendChild(btnLast);
    controls.appendChild(pager);

    tbl.parentNode.insertBefore(controls, tbl);

    function render() {
      var totalPages = Math.max(1, Math.ceil(visibleRows.length / pageSize));
      if (page >= totalPages) page = totalPages - 1;
      if (page < 0) page = 0;
      rows.forEach(function (r) { r.style.display = 'none'; });
      var start = page * pageSize;
      var end = Math.min(start + pageSize, visibleRows.length);
      for (var i = start; i < end; i++) {
        visibleRows[i].style.display = '';
      }
      info.textContent = ' page ' + (page + 1) + ' / ' + totalPages +
        ' · ' + visibleRows.length + ' rows ';
      btnFirst.disabled = btnPrev.disabled = page === 0;
      btnLast.disabled = btnNext.disabled = page >= totalPages - 1;
    }

    function applyFilters() {
      var q = search.value.trim().toLowerCase();
      visibleRows = rows.filter(function (r) {
        if (q && r.textContent.toLowerCase().indexOf(q) === -1) return false;
        for (var k = 0; k < filterSelects.length; k++) {
          var f = filterSelects[k];
          if (!f.sel.value) continue;
          var c = r.children[f.idx];
          if (!c || c.textContent.trim() !== f.sel.value) return false;
        }
        return true;
      });
      render();
    }

    search.addEventListener('input', function () { page = 0; applyFilters(); });
    btnFirst.addEventListener('click', function () { page = 0; render(); });
    btnPrev.addEventListener('click', function () { page--; render(); });
    btnNext.addEventListener('click', function () { page++; render(); });
    btnLast.addEventListener('click', function () { page = Infinity; render(); });

    headers.forEach(function (th, idx) {
      if (th.dataset.nosort != null) return;
      th.dataset.sort = '';
      th.addEventListener('click', function () {
        var dir = th.dataset.sort === 'asc' ? 'desc' : 'asc';
        headers.forEach(function (h) { if (h !== th) h.dataset.sort = ''; });
        th.dataset.sort = dir;
        var sign = dir === 'asc' ? 1 : -1;
        rows.sort(function (a, b) {
          var av = cellSortValue(a.children[idx]);
          var bv = cellSortValue(b.children[idx]);
          return sign * compare(av, bv);
        });
        rows.forEach(function (r) { tbody.appendChild(r); });
        applyFilters();
      });
    });

    // Initial default sort marker (data-default-sort="N|desc")
    var def = tbl.dataset.defaultSort;
    if (def) {
      var bits = def.split('|');
      var col = parseInt(bits[0], 10);
      var dir = bits[1] === 'asc' ? 'asc' : 'desc';
      var th = headers[col];
      if (th) {
        th.dataset.sort = dir;
        var sign = dir === 'asc' ? 1 : -1;
        rows.sort(function (a, b) {
          var av = cellSortValue(a.children[col]);
          var bv = cellSortValue(b.children[col]);
          return sign * compare(av, bv);
        });
        rows.forEach(function (r) { tbody.appendChild(r); });
      }
    }

    applyFilters();
  }
  function init() {
    document.querySelectorAll('table.paginated').forEach(setupTable);
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
"""


__all__ = ["COHORT_CSS", "COHORT_JS", "TIER_PICKER_CSS", "TIER_PICKER_JS"]


# ---------------------------------------------------------------------------
# Confidence-tier picker (inline radio/select toggle)
# ---------------------------------------------------------------------------
# Renders a small selector at the top of the cohort body that swaps the
# visible <section.csc-tier> block on change.  Both sections (sensitive
# + high-confidence) remain in the DOM so the user can flip back and
# forth without losing scroll position; we just toggle their `display`
# style.  Falls back gracefully when JS is disabled (the first section
# is visible by default; subsequent ones are hidden via inline style
# applied at first-render).

TIER_PICKER_CSS = r"""
.tier-picker { background: #f1f6ff; border: 1px solid #c2d4f0;
               border-radius: 4px; padding: 0.5em 0.8em;
               margin: 0.6em 0 0.8em; font-size: 13px;
               position: sticky; top: 0; z-index: 5; }
.tier-picker label { font-size: 13px; }
.tier-picker select { margin-left: 0.4em; font-size: 13px; padding: 2px 6px; }
.tier-picker .tier-picker-help { color: #444; font-size: 12px;
                                  margin-left: 0.6em; }
section.csc-tier { display: none; }
section.csc-tier.csc-tier-active { display: block; }
.tier-banner { background: #fffbe6; border-left: 3px solid #d9b300;
               padding: 0.3em 0.6em; margin: 0 0 0.6em;
               font-size: 12.5px; color: #5a4500; }
h2.confidence-concordance { border-top: 2px solid #888;
                             padding-top: 0.4em; margin-top: 1.2em; }
"""

TIER_PICKER_JS = r"""
(function () {
  var sel = document.getElementById('csc-tier-select');
  if (!sel) return;
  var sections = document.querySelectorAll('section.csc-tier');
  function activate(value) {
    sections.forEach(function (s) {
      var match = s.dataset.tier === value;
      s.classList.toggle('csc-tier-active', match);
    });
  }
  // Activate the first option on load (sensitive tier).
  if (sections.length > 0) {
    activate(sections[0].dataset.tier);
    sel.value = sections[0].dataset.tier;
  }
  sel.addEventListener('change', function () { activate(sel.value); });
})();
"""
