function renderSegmentPlot(){
    if(selectedProtein === ""){
        d3.select("#graphDiv3").select("svg").remove();
        document.getElementById('no_protein_div_3').innerHTML = '<div style="width: 640px; height: 440px;"><h3>Select a protein to get started.</h3></div>';
        return;
    } else {
        document.getElementById('no_protein_div_3').innerHTML = "<div></div>";
    }

    // remove old svg
    d3.select("#graphDiv3").select("svg").remove();

    d3.json("http://127.0.0.1:5000/get-segment-data?protein="+selectedProtein+"&samples="+selectedSample, function(error, data) {
        if (error) throw error;
        var rect_patches = data.peptide_patches;
        var mod_patches = data.mod_patches;
        var plot_height = data.height;
        var peptide_seq = data.seqq;
        var modification_color_map = data.modification_color_map;
        var min_peptide = data.min_peptide;
        var max_peptide = data.max_peptide;
        var min_intensity = min_peptide[0];
        var max_intensity = max_peptide[0];
        //round min_intensity and max_intensity to 10 decimals
        min_intensity = Math.round(min_intensity * 100000000) / 100000000;
        max_intensity = Math.round(max_intensity * 100000000) / 100000000;

        var min_color = min_peptide[1];
        var max_color = max_peptide[1];
        var color = d3.scaleLinear().range([min_color, max_color]).domain([1, 2]);
        
        modification_color_map_keys = Object.keys(modification_color_map);
        modification_color_map_values = Object.values(modification_color_map);
        // create map from values to keys
        var colors_to_mod_map = new Map();
        for (var i = 0; i < modification_color_map_keys.length; i++) {
            colors_to_mod_map.set(modification_color_map_values[i], modification_color_map_keys[i]);
        }

    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 30, bottom: 30, left: 30},
        width = screen.width - margin.left - margin.right,
        height = 700 - margin.top - margin.bottom;
        
    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv3")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    // make scrolable svg
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    
    var scale_factor = width/peptide_seq.length;
    // var scale_factor = (width - margin.left - margin.right)/peptide_length;

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;
    

    var xScale = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([0, width]);
    // based on width of window calculate tick_values_distance
    // var tick_values_distance = Math.floor(width/peptide_length);
    var tick_values_distance = 3;
    var xAxis = d3.axisTop(xScale)
        // .tickValues(d3.range(0, peptide_length, tick_values_distance))
        .tickValues(d3.range(0, peptide_length))
    
        .tickFormat(function(d,i){ return peptide_chars[d];
            // if (i%tick_values_distance == 0){
            //     return peptide_chars[d];
            // }
            // else{
            //     // set tick size
            //     return "";
            // }
        });        
        xAxis.tickSizeOuter(12);
        
        // TODO: adjust when zooming
    svg.append("g")
        .attr("transform", "translate(" + 0 + ")")
        .call(xAxis);    
        
    var rx = 2
    var ry = 2
    var stroke_width = 0.2

    // sample data for rectangles
    var rects = svg.selectAll("foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*scale_factor)
        .attr("y", d=> d[1]*scale_factor)
        .attr("width", d=> d[2]*scale_factor)
        .attr("height", d=> d[3]*scale_factor)
        .attr("fill", d=> d[4])
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", 0.8)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)
        .on("mouseover",mouseover)
        .on("mousemove", mousemove_segments)
        .on("mouseout", mouseleave);
        

    var tooltip = d3.select("#graphDiv3")
        .append("div")
        .style("opacity", 0)
        .attr("class", "tooltip")
        .style("background-color", "white")
        .style("border", "solid")
        .style("border-width", "2px")
        .style("border-radius", "5px")
        .style("padding", "5px")

    // Three function that change the tooltip when user hover / move / leave a cell
    function mouseover(d) {
        tooltip.style("opacity", 1)
    }
    function mousemove_modification(d) {
        // get modification type from colors_to_mod_map map
        var mod_type = colors_to_mod_map.get(d[4]);
        var mod_position = d[0];
        var mod_char = peptide_chars[mod_position];
        tooltip.html("<p>Intensity: " + expo(d[5], 3) + "</p><p>Modification Type: " + mod_type + "</p><p> Position: " + mod_position + "</p><p>Modification Character: " + mod_char + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
    }
    function mousemove_segments(d) {
        tooltip.html("<p>Peptide: " + peptide_seq.substring(d[0],  d[0] + d[2]) + "</p><p>Intensity: " + expo(d[5], 3) + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
    }
    
    function mouseleave(d) {
        tooltip.style("opacity", 0)
    }
    function expo(x, f) {
        return Number.parseFloat(x).toExponential(f);
      }


    var mod_rects = svg.selectAll("boo")
        .data(mod_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*scale_factor)
        .attr("y", d=> (d[1])*scale_factor)
        .attr("width", d=> d[2]*scale_factor)
        .attr("height", d=> d[3]*scale_factor)
        .attr("fill", d=> d[4])
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)
        .on("mouseover",mouseover)
        .on("mousemove", mousemove_modification)
        .on("mouseout", mouseleave);

    // Inspired by: https://www.d3-graph-gallery.com/graph/custom_legend.html
    //   modification labels
    var size = 10
    var legend = svg.selectAll("legend")
        .data(modification_color_map_keys)
        .enter()
        .append("rect")
        .attr("x", 0)
        .attr("y", function(d,i){ return 100 + i*(size+5)}) // 100 is where the first dot appears. 25 is the distance between dots
        .attr("width", size)
        .attr("height", size)
        .style("fill", function(d, i){ return modification_color_map[d]})
    
    // Add one dot in the legend for each name.
    legend_text = svg.selectAll("myLabels")
        .data(modification_color_map_keys)
        .enter()
        .append("text")
        .attr("x", 0 + size*1.2)
        .attr("y", function(d,i){ return 100 + i*(size+5) + (size/2)}) // 100 is where the first dot appears. 25 is the distance between dots
        .text(function(d){ return d})
        .attr("text-anchor", "left")
        .style("fill", "black")
        .style("alignment-baseline", "middle")
        .on("mouseover", function(d,i){
            d3.select(this).style("font-weight", "bold");
        })
        .on("mouseout", function(d,i){
            d3.select(this).style("font-weight", "normal");
        })
        .on("click", function(d,i){
            // switch between color in color map and standard color 
            var current_color = d3.select(this).style("fill");  
            if (current_color == "black"){ 
                d3.select(this).style("fill", modification_color_map[d]);
                d3.select(this).style("font-weight", "bold");
            }
            else{
                d3.select(this).style("fill", "black");
                d3.select(this).style("font-weight", "normal");
            }
            })


    // move legend and legend_text  to Center LEFT 
    legend.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");
    legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");

    color_bar_width = 20;
    color_bar_height = 180;

    // find difrence in magnintude and use as steps
    var max_exp = expo(max_intensity,1);
    var min_exp = expo(min_intensity,1);
    var max_exp_last_char = max_exp.substring(max_exp.length-1);
    var min_exp_last_char = min_exp.substring(min_exp.length-1);
    var steps = Math.abs(max_exp_last_char - min_exp_last_char)+1;

    var step = (Math.log(max_intensity) - Math.log(min_intensity)) / (steps - 1);
    log_steps = [];
    for (var i = 0; i < steps; i++) {
        log_steps.push(Math.exp(Math.log(min_intensity) + i * step));
    }

    //parse log_steps to string
    var log_steps_string = [];
    for (var i = 0; i < log_steps.length; i++) {
        log_steps_string.push(expo(parseFloat(log_steps[i]).toPrecision(3), 3));
    }

    //add max_intensity to log_steps_floor
    // make firsttt item equal expo(max_intensity,3)
    // delete first and last item in log_steps_string
    log_steps_string = log_steps_string.slice(1, log_steps_string.length - 1);
    
    log_steps_string.push(expo(max_intensity,3));
    log_steps_string.unshift(expo(min_intensity,3));
    
    var defs = svg.append("defs");
    var linearGradient = defs.append("linearGradient")
        .attr("id", "linear-gradient")
        .attr("x1", "0%")
        .attr("y1", "0%")
        .attr("x2", "0%")
        .attr("y2", "100%");

    linearGradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", color(2))

    linearGradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", color(1))

    var rect = svg.append("rect")
        .attr("x", 200)
        .attr("y", 0)
        .attr("width", color_bar_width)
        .attr("height", color_bar_height)
        .style("fill", "url(#linear-gradient)");

    color_legend_text = svg.selectAll("colorLabels")
        .data(log_steps_string)
        .enter()
        .append("text")
        .attr("x", 220)
        .attr("y", function(d,i){ return color_bar_height - i*(color_bar_height/(steps-1))})
        .text(function(d){ return "-" +d})
        .attr("text-anchor", "right")
        .style("alignment-baseline", "middle")

        //move color_legend_text next to rect
    color_legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    rect.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    });
}

$('document').ready(function(){
    renderSegmentPlot();
});
