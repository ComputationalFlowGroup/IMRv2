// MathJax setup for pymdownx.arithmatex (generic mode): equations written as
// Subscribe only after the CDN script has initialized MathJax. Material's
// document$ emits immediately; subscribing at config load can otherwise
// access startup.output before it exists and leave equations unrendered.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*",
    processHtmlClass: "arithmatex"
  },
  startup: {
    typeset: false,
    ready: () => {
      MathJax.startup.defaultReady();
      document$.subscribe(() => {
        // Serialize page changes with initialization and previous rendering.
        MathJax.startup.promise = MathJax.startup.promise.then(() => {
          MathJax.typesetClear();
          MathJax.texReset();
          return MathJax.typesetPromise();
        }).catch(error => console.error("MathJax typesetting failed:", error));
      });
    }
  }
};
