
function ct(ct_data, ct_scale,container) {
	var ct_data = ct_data; 
    //var hits = hits;
    //var hit_colors = hit_colors;
    //var hmms = hmms;
    //var hmm_colors = hmm_colors;
    
  $("#" + container).html("");
  var ui = $("#" + container);
    ui.append('<div id="ctHits">')
    var ctHits = $("#ctHits");
    ctHits.append("<h3>Cluster Comparison</h3>")
  for (var i in ct_data) {
    ctHits.append('<div id="hit-' + i + '">')
    var hitDiv = $("#hit-"+i);
    if (ct_data[i].hasOwnProperty("dist")) {
        hitDiv.append("<h4>" + ct_data[i]["id"] + " - Distance to Top Cluster: " + ct_data[i]["dist"] + " - (Location: " + ct_data[i]["start"] + "-" + ct_data[i]["end"] + ", Size: " + ct_data[i]["size"] + " nts)" + " </h4>");}
        else{
  	hitDiv.append("<h4>" + ct_data[i]["id"] + " (Location: " + ct_data[i]["start"] + "-" + ct_data[i]["end"] + ", Size: " + ct_data[i]["size"] + " nts)" 
                  + "</h4>");}
    var svg = $(Arrower.drawClusterToolsSVG(ct_data[i],ct_scale));
    svg.css("clear", "both");
    svg.addClass("arrower-svg");
    hitDiv.append(svg);
    hitDiv.find("svg.arrower-svg").each(
    	function(){
    		$(this).attr("width",$(this).find("g")[0].getBoundingClientRect().width);
    		$(this).attr("height",$(this).find("g")[0].getBoundingClientRect().height*1.3);
    	});
  }
    
/*ui.append('<div style="text-align:right" class="desc-container" id="legend_container">');
var legendCont = $('#legend_container');
legendCont.append("<p style='display:inline;font-size:x-small'><a title='Show Legend' id='toggle' href='##'align='right'>Toggle Legend</a></p>");

legendCont.append('<div style="text-align:left" class="desc-ui hidden" style="margin-top: 2px;" id="legend">');

var legendDiv = $('#legend');

legendDiv.append("<h4 style='display:inline'>Gene Hits</h4>");
legendDiv.append('<div id="geneQueries">');
    var geneQueries = $("#geneQueries");
    for (var i in hits){
        geneQueries.append('<div id="geneHit-' + i + '">');
        var geneQuery = $("#geneHit-"+i);
        geneQuery.append('<svg width="20px" height="20px"><rect width="20px" height="20px" style="fill:' + hit_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>');
        geneQuery.append('\t'+hits[i]+'\t');
    }
    legendDiv.append("<h4 style='display:inline'>HMM Hits</h4>");
    legendDiv.append('<div id="hmmQueries">');
    var hmmQueries = $("#hmmQueries");
    for (var i in hmms){
        hmmQueries.append('<div  id="hmmHit-' + i + '">');
        var hmmQuery = $("#hmmHit-"+i);
        hmmQueries.append('<svg width="20px" height="20px"><rect width="20px" height="20px" style="fill:' + hmm_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>');
        hmmQueries.append('\t'+hmms[i]+'\t');
    }
   
  var toggle_btn = $("#toggle");
  toggle_btn.click(function(handler){
    if (legendDiv.hasClass("hidden")) {
      legendDiv.removeClass("hidden");
    } else {
      legendDiv.addClass("hidden")
      }
    handler.stopPropagation();
  });*/
}