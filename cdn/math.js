// math.js -- MathJax for any page in the Catalog, at any depth.
//
// Loaded on its own by the appendix pages (<script src="../cdn/math.js">),
// and by cdn.js for the pages at the root.  Keeping it in one file means the
// delimiters are configured in exactly one place: the Catalog writes inline
// math as $...$, which is NOT a MathJax default and has to be asked for.
document.write(
'<script type="text/x-mathjax-config">'+
'                MathJax.Hub.Config({"HTML-CSS": { preferredFont: "TeX", availableFonts: ["STIX","TeX"], linebreaks: { automatic:true }, EqnChunk: (MathJax.Hub.Browser.isMobile ? 10 : 50) },'+
'                    tex2jax: { inlineMath: [ ["$", "$"], ["\\\\(","\\\\)"] ], displayMath: [ ["$$","$$"], ["\\[", "\\]"] ], processEscapes: true, ignoreClass: "tex2jax_ignore|dno" },'+
'                    TeX: {'+
'                        extensions: ["begingroup.js"],'+
'                        noUndefined: { attributes: { mathcolor: "red", mathbackground: "#FFEEEE", mathsize: "90%" } },'+
'                        Macros: { href: "{}" }'+
'                    },'+
'                    displayAlign: "left",'+
'                    messageStyle: "none",'+
'                    styles: { ".MathJax_Display, .MathJax_Preview, .MathJax_Preview > *": { "background": "inherit" } },'+
'                    SEEditor: "mathjaxEditing"'+
'            });'+
'            </script>'+
'<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS_HTML-full"></script>'
);
