const { zoom } = require("d3-zoom");

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
    
    // get <p> paragraph field inside the graphDiv3 div to change Peptide Segments Plot text
    var p = document.getElementById("graphDiv3").getElementsByTagName("p")[0];
    p.innerHTML = "Peptide Segments Plot - Protein: " + selectedProtein;
    p.style.fontWeight = "bold";
    p.style.fontStyle = "italic";

    
    
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

    
    var segment_height = 6;
    var values_distance = 10;

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;

    var isScrollDisplayed = values_distance * peptide_length > width;
    console.log("isScrollDisplayed", isScrollDisplayed);

    var xScale_ticks = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([0, peptide_length*values_distance]);
    
    var xScale_labels = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([values_distance/2, peptide_length*values_distance + values_distance/2]);

    var xAxis_ticks = d3.axisTop()
        .scale(xScale_ticks).tickPadding(10)
        .tickValues(d3.range(0, peptide_length))
        // .tickFormat(function(d,i){ return peptide_chars[d];
        .tickFormat(function(d,i){ return null
        });        
        xAxis_ticks.tickSizeOuter(12);
    
    var xAxis_labels = d3.axisTop()
        .scale(xScale_labels).tickPadding(8)
        .tickSize(0)
        .tickValues(d3.range(0, peptide_length))
        .tickFormat(function(d,i){ 
            return peptide_chars[d];
        });
        
    svg.append("g")
        .attr("class", "xAxis_ticks")
        .attr("transform", "translate(0, -1)")
        .call(xAxis_ticks);
    
    svg.append("g")
        .attr("class", "xAxis_labels")
        .attr("transform", "translate(0, -1)")
        .call(xAxis_labels)
        .selectAll("text")
        .style("font-weight", function(d,i){
            if(i%10 === 0){
                return "bold";
            } else {
                return "normal";
            }
        })
        .style("font-size", function(d,i){
            if(i%10 === 0){
                return "1.2em";
            } else {
                return "1em";
            }
        });

    var rx = 3
    var ry = 3
    var stroke_width = 0.1
    var opacity = 0.7
    var opacity_mod = 0.8
    opacity_mod_grey = 0.5
    // sample data for rectangles
    var rects = svg.selectAll("foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*values_distance)
        .attr("y", d=> d[1]*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        .attr("fill", d=> d[4])
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity)
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
        if (d[0]-3  >= 0) {
            var three_chars_before = peptide_chars.slice(d[0]-3, d[0]); 
        } else {
            var three_chars_before = peptide_chars.slice(0, d[0]);
        }
        if (d[0]+ d[2]+3 <= peptide_length) {
            var three_chars_after = peptide_chars.slice(d[0] + d[2], d[0] + d[2]+3);
        } else {
            var three_chars_after = peptide_chars.slice(d[0] + d[2], peptide_length);
        }
        // make text smaller in html
        three_chars_before_italic = '<i>' +three_chars_before.join("")+'</i>';
        three_chars_after_italic = '<i>' +three_chars_after.join("") +'</i>';
        
        // get all modification on this segments
        mod_types_and_positions = [];
        mod_positions = [];
        for (var i = 0; i < mod_patches.length; i++) {
            if (mod_patches[i][0] >= d[0] && mod_patches[i][0] <= d[0] + d[2] && mod_patches[i][1] == d[1]) {
                pos = mod_patches[i][0];
                mod_positions.push(pos);
                _type =colors_to_mod_map.get(mod_patches[i][4]);
                mod_types_and_positions.push(_type + "(" + pos+")");
            }
        }
        
        mod_types_and_positions_str = mod_types_and_positions.join(", "); 
        if (mod_types_and_positions_str.length > 0) {
        tooltip.html("<p>Peptide: " + three_chars_before_italic+ '.' + peptide_seq.substring(d[0],  d[0] + d[2])+ '.' + three_chars_after_italic + "</p><p>Intensity: " + expo(d[5], 3) + "</p><p>Modifications: " + mod_types_and_positions_str + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        } else {
            tooltip.html("<p>Peptide: " + three_chars_before_italic+ '.' + peptide_seq.substring(d[0],  d[0] + d[2])+ '.' + three_chars_after_italic + "</p><p>Intensity: " + expo(d[5], 3) + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        }

        // highlight the peptide sequence in the x-axis 
        var x_axis_highlight = svg.selectAll(".xAxis_labels")
            .selectAll("text")
            .style("fill", function(dd,i){
                if(i >= d[0] && i < d[0] + d[2]){
                    // check if modification at position i
                    if (mod_positions.includes(i)) {
                        return "red";
                    } else {
                        return d[4];
                    }
                } else {
                    return "black";
                }
            })
            // increase size of the text
            .style("font-size", function(dd,i){
                if(i >= d[0] && i < d[0] + d[2]){
                    return "1.5em";
                } else {
                    return "1em";
                }
            })
            // make bold
            .style("font-weight", function(dd,i){
                if(i >= d[0] && i < d[0] + d[2]){
                    return "bold";
                } else {
                    return "normal";
                }
            });

    }
    
    function mouseleave(d) {
        tooltip.style("opacity", 0)
        var x_axis_highlight = svg.selectAll(".xAxis_labels")
        .selectAll("text")
        .style("fill", 'black');
        // TODO: make it back to normal
        // ##########
        // ##########
        // ##########
        // ##########
    }
    function expo(x, f) {
        return Number.parseFloat(x).toExponential(f);
      }


    var mod_rects = svg.selectAll("boo")
        .data(mod_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*values_distance)
        .attr("y", d=> (d[1])*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        // .attr("fill", d=> d[4])
        .attr("fill", 'grey')
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity_mod)
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
        .on("mouseover", function(d,i){
            d3.select(this)
            .style("cursor", "pointer")
            .style("stroke", "black")
            .style("stroke-width", "2px");
        })
        .on("mouseout", function(d,i){
            d3.select(this)
            .style("cursor", "none")
            .style("stroke", "none");
        })
        // implement ONCLICK
        .on("click", function(d,i){
            // iterate over all modifications where the color is the same as the one clicked
            // and change the opacity to 0
            mod_rects.filter(function(dd){
                color_val = modification_color_map[d];
                if (color_val == dd[4]) {
                    // check color
                    col = d3.select(this).style("fill");
                    console.log(col);
                    if (col == "rgb(128, 128, 128)") {// if grey
                        console.log("color is grey");
                        d3.select(this).style("fill", color_val);
                        d3.select(this).attr("opacity", opacity_mod);
                    } else {
                        d3.select(this).style("fill", "grey");
                        d3.select(this).attr("opacity", opacity_mod_grey);
                    }
                    // opacity_ = d3.select(this).attr("opacity");
                    // // change opacity
                    // if (opacity_ == 0) {
                    //     d3.select(this).attr("opacity", opacity_mod);
                    // } else {
                    //     d3.select(this).attr("opacity", 0);
                    // }
                }
            })
        });
    
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
            d3.select(this).style("font-weight", "bold").style("cursor", "pointer");
        })
        .on("mouseout", function(d,i){
            d3.select(this).style("font-weight", "normal").style("cursor", "none");
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
        log_steps_string.push(expo(parseFloat(log_steps[i]).toPrecision(2), 2));
    }

    log_steps_string = log_steps_string.slice(1, log_steps_string.length - 1);
    
    log_steps_string.push(expo(max_intensity,2));
    log_steps_string.unshift(expo(min_intensity,2));
    
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


    var color_legend_x = 220;
    color_legend_text = svg.selectAll("colorLabels")
        .data(log_steps_string)
        .enter()
        .append("text")
        .attr("x", color_legend_x)
        .attr("y", function(d,i){ return color_bar_height - i*(color_bar_height/(steps-1))})
        .text(function(d){ return "-" +d})
        .attr("text-anchor", "right")
        .style("alignment-baseline", "middle")

    // add right line to next to color_legend_text
    var legendLine = svg.append("line")
    .attr("x1", color_legend_x)
    .attr("y1", 0)
    .attr("x2", color_legend_x)
    .attr("y2", color_bar_height)
    .attr("stroke", "black")
    .attr("stroke-width", 1);
    
    //move color_legend_text next to rect
    color_legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    // move legendLine next to color_legend_text
    legendLine.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");

    
    rect.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    });
}

$('document').ready(function(){
    renderSegmentPlot();
});
