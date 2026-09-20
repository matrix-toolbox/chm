// cdn.js -- everything the pages at the root of the Catalog load.
//
// The paths are resolved against this file rather than against the including
// page, so the script also works from a subdirectory (an appendix page needs
// only <script src="../cdn.js">).  Appendix pages that want maths alone, and
// not the Bootstrap layer, load cdn/math.js instead.
(function () {
    var here = document.currentScript.src.replace(/cdn\.js(\?.*)?$/, '');
    document.write(
    '<script src="' + here + 'cdn/math.js"><\/script>'+
    '<script src="' + here + 'cdn/jquery.min.js"><\/script>'+
    '<link rel="stylesheet" href="' + here + 'cdn/bootstrap_CHM.min.css">'+
    '<script src="' + here + 'cdn/bootstrap.min.js"><\/script>'+
    '<link rel="stylesheet" href="' + here + 'cdn/chm.css">'+
    '<link rel="icon" type="image/png" href="' + here + 'cdn/icon.png">'
    );
})();
