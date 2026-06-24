window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  // Do not restrict to .arithmatex elements — docstring math rendered by
  // mkdocstrings uses \(...\) / \[...\] delimiters but does not carry the
  // arithmatex class, so MathJax must scan the full page.
  options: {
    ignoreHtmlClass: "md-code|highlight",
  },
};
