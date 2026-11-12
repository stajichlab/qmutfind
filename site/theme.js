// Theme buttons: set the theme, save the choice, and mark the active button.
(function () {
  var root = document.documentElement;
  function mark() {
    document.querySelectorAll('[data-set-theme]').forEach(function (b) {
      b.setAttribute('aria-pressed', String(b.getAttribute('data-set-theme') === root.getAttribute('data-theme')));
    });
  }
  document.querySelectorAll('[data-set-theme]').forEach(function (button) {
    button.addEventListener('click', function () {
      var theme = button.getAttribute('data-set-theme');
      root.setAttribute('data-theme', theme);
      try { localStorage.setItem('qmutfind-theme', theme); } catch (e) {}
      mark();
      document.dispatchEvent(new Event('qmutfind-theme'));
    });
  });
  mark();
})();
