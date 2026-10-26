// The README figures come in a light and a dark version inside <picture>, and on GitHub the
// reader's system theme picks one. On this site the theme is set by the light switch, which
// writes data-bs-theme on <html>, so the figures follow that attribute instead.
(function () {
  function sync() {
    var dark = document.documentElement.getAttribute("data-bs-theme") === "dark";
    document.querySelectorAll("picture").forEach(function (p) {
      var img = p.querySelector("img");
      if (!img) return;
      if (!img.dataset.light) {
        // the first pass keeps both paths and drops the <source>, so that the system theme
        // no longer competes with the switch
        var source = p.querySelector('source[media*="prefers-color-scheme"]');
        if (!source) return;
        img.dataset.light = img.getAttribute("src");
        img.dataset.dark = source.getAttribute("srcset");
        source.remove();
      }
      img.src = dark ? img.dataset.dark : img.dataset.light;
    });
  }
  document.addEventListener("DOMContentLoaded", function () {
    sync();
    new MutationObserver(sync).observe(document.documentElement,
      { attributes: true, attributeFilter: ["data-bs-theme"] });
  });
})();
