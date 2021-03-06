function renderSegmentPlot(){
    if(getSelectedProtein() === ""){
        d3.select("#graphDiv3").select("svg").remove();
        document.getElementById('no_protein_div_3').innerHTML = '<div style="width: 640px; height: 440px;"><h3>Select a protein to get started.</h3></div>';
        return;
    } else {
        document.getElementById('no_protein_div_3').innerHTML = "<div></div>";
    }
    
    // remove old svg
    d3.select("#graphDiv3").select("svg").remove();
    
    // get <p> paragraph field inside the graphDiv3 div to change Peptide Segments Plot text
    var selectedSamplesKey = "all";
    if(getSelectedSamples().length > 0){
        selectedSamplesKey = getSelectedSamples();
    }
    var p = document.getElementById("graphDiv3").getElementsByTagName("p")[0];
    p.innerHTML = "<h5>Peptide Segments Plot</h5><h6>Protein: " + getSelectedProtein()+"; Samples: "+getSelectedSamples()+"</h6>";
    
    
    d3.json("http://127.0.0.1:5000/get-segment-data?protein="+getSelectedProtein()+"&samples="+getSelectedSamples(), function(error, data) {
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
    var histogram_data = data.histogram_data;
    //round min_intensity and max_intensity to 10 decimals
    var min_intensity = Math.round(min_intensity * 100000000) / 100000000;
    var max_intensity = Math.round(max_intensity * 100000000) / 100000000;

    var min_color = min_peptide[1];
    var max_color = max_peptide[1];
    var color = d3.scaleLinear().range([min_color, max_color]).domain([1, 2]);

    var modification_color_map_keys = Object.keys(modification_color_map);
    var modification_color_map_values = Object.values(modification_color_map);
    // create map from values to keys
    var colors_to_mod_map = new Map();
    for (var i = 0; i < modification_color_map_keys.length; i++) {
        colors_to_mod_map.set(modification_color_map_values[i], modification_color_map_keys[i]);
    }

    var selector_height = 100
    var histogram_height = selector_height / 4;

    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 20, bottom: 30, left: 20};
    var margin_overview = {top: 30, right: 15, bottom: 10, left: 15};
    var width = screen.width * 0.9;
    var mod_label_size = 10
    var color_bar_width = 20;
    var color_bar_height = 180;
    var height = plot_height*5 + margin_overview.bottom + 12 + color_bar_height;
        
    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv3")
    .append("svg")
    .attr("width", width)
    .attr("height", height + selector_height + histogram_height + 50)
    .append("g")
    
    var segment_height = 5;
    var values_distance = 10;

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;

    var xScale_ticks = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([0, peptide_length*values_distance]);
    
    var xScale_labels = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([values_distance/2, peptide_length*values_distance + values_distance/2]);

    var xAxis_ticks = d3.axisTop()
        .scale(xScale_ticks).tickPadding(10)
        .tickValues(d3.range(0, peptide_length))
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
        
    var x_ticks = svg.append("g")
        .attr("class", "xAxis_ticks")
        .attr("transform", "translate(" + 0 + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
        .call(xAxis_ticks);
    
    var x_labels = svg.append("g")
        .attr("class", "xAxis_labels")
        .attr("transform", "translate(" + 0 + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
        .call(xAxis_labels);

    // make every other label smaller
    var adjust_xaxis = function (scale) {
        // make every second label have font size by scale smaller
        x_labels.attr("font-size", function(d,i){
            if (i % 2 == 0) {
                return scale * mod_label_size;
            } else {
                return mod_label_size;
            }
        });
    }

    //make zoom fuctionality for x axis
    
    // TODO
    // adjust_xaxis(2);
    
    
    var rx = 2
    var ry = 2
    var stroke_width = 0.1
    var opacity = 0.7
    var opacity_mod = 0.8
    var opacity_mod_grey = 0.5
    var opacity_highlight_legend = 0.85
    // sample data for rectangles


    var rects = svg.selectAll("foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*values_distance)
        .attr("y", d => (histogram_height + selector_height + margin_overview.bottom +12) + d[1]*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        .attr("fill", d=> d[4])
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)
        .on("mouseover",mouseover)
        .on("mouseenter", mousemove_segments)
        .on("click", mouseclick)
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
    function mouseover(d){
        tooltip.style("opacity", 1);
    }

    function mouseover_mod(d) {
        makeModificationTooltip(colors_to_mod_map, d, tooltip);
        
        highlightSegmentAndModInXaxis(rect_patches, d, mod_patches, colors_to_mod_map, peptide_seq, svg);

        // COLORING of modificaiton in legend and in segment
        d3.select(this)  // Highlight color of mod in segment
           .style("fill", d[4])
           .attr("opacity", opacity_mod);

        var mod_type = colors_to_mod_map.get(d[4]);
        // Highlight color in  mod legend
        var legend_rects_ = d3.selectAll(".legend_rects");
        var legend_rects_filtered = legend_rects_.filter(function(rect) {
            return (rect == mod_type)
        });
        var legend_rects_filtered_fill = legend_rects_filtered.style("fill");
        var legend_rects_filtered_fill_hex = d3.rgb(legend_rects_filtered_fill).toString();
        var mod_rect_hex_color = d3.rgb(d[4]).toString();
        // if color of mod rect is NOT same as color of mod legend rect change legend rect color
        if (legend_rects_filtered_fill_hex != mod_rect_hex_color){ 
            legend_rects_filtered.style("opacity", opacity_highlight_legend);
            legend_rects_filtered.style("fill", d[4]);
        }
    }

    function mouseclick (d){
        try{
            removeProteinSelectionPlot();
        }catch(e){}
        var modPositionsAndTypes = getModificationPositions(mod_patches, d, colors_to_mod_map, peptide_seq);
        var mod_types_and_positions = modPositionsAndTypes[0];
        var mod_positions = modPositionsAndTypes[1];
        var mod_types_and_positions_str = mod_types_and_positions.join(", "); 
        proteinStartPos = d[0]+1;
        proteinEndPos = proteinStartPos+d[2] - 1;
        selectedSequence=data.seqq.substring(proteinStartPos-1, proteinEndPos);
        toolTipWidth = selectedSequence.length*15;
        if(toolTipWidth > 600){
            toolTipWidth = 600;
        }
        if (toolTipWidth < 150) {
            toolTipWidth = 150;
        } 
        tooltip.html('<div class="container" style="width: '+toolTipWidth+'px;"><p>Intensity: ' + expo(d[5], 3) + '</p><p>Modification: ' + mod_types_and_positions_str + '</p><div id="peptide_selection_view"></div></div>')
        .style("left", (d3.event.pageX + 10) + "px")
        .style("top", (d3.event.pageY - 10) + "px");
        if (mod_types_and_positions_str.length > 0) {
        tooltip.html('<div class="container" style="width: '+toolTipWidth+'px;"><p>Intensity: ' + expo(d[5], 3) + '</p><p>Modification: ' + mod_types_and_positions_str + '</p><div id="peptide_selection_view"></div></div>')
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        } else {
            tooltip.html('<div class="container" style="width: '+toolTipWidth+'px"><p>Intensity: ' + expo(d[5], 3) + '</p><div id="peptide_selection_view"></div></div>')
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        }
        renderProteinSelectionPlot();
    }

    function mousemove_segments(d) {
        // get all modification on this segments - retuns [mod_types_and_positions, mod_positions]
        var modPositionsAndTypes = getModificationPositions(mod_patches, d, colors_to_mod_map, peptide_seq);
        var mod_types_and_positions = modPositionsAndTypes[0];
        var mod_positions = modPositionsAndTypes[1];
        var mod_types_and_positions_str = mod_types_and_positions.join(", "); 
        var tempProteinStartPos = d[0]+1;
        var tempProteinEndPos = proteinStartPos+d[2];
        var tempSelectedSequence=data.seqq.substring(tempProteinStartPos, tempProteinEndPos);
        var tempToolTipWidth = tempSelectedSequence.length*15;
        if(tempToolTipWidth > 600){
            tempToolTipWidth = 600;
        }
        if (tempToolTipWidth < 150) {
            tempToolTipWidth = 150;
        } 
        if (mod_types_and_positions_str.length > 0) {
        tooltip.html('<div class="container" style="width: '+tempToolTipWidth+'px"><p>Intensity:' + expo(d[5], 3) + "</p><p>Modifications: " + mod_types_and_positions_str + "</p></div>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        } else {
            tooltip.html('<div class="container" style="width: '+tempToolTipWidth+'px"><p>Intensity:'  + expo(d[5], 3) + "</p></div>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        }
        // highlight the peptide sequence in the x-axis 
        var x_axis_highlight = highlightSeqInXAxis(svg, d, mod_positions, 0, false);
    }
    
    function mouseleave(d) {
        tooltip.style("opacity", 0)
        setXaxisNormal(svg);
    }

    function mouseleave_mod(d) {
        tooltip.style("opacity", 0)
        var x_axis_highlight = setXaxisNormal(svg);
        var mod_fill = d3.select(this).style("fill");
        if (mod_fill != "grey" && mod_fill != "rgb(128, 128, 128)") {
            var legend_rects_ = d3.selectAll(".legend_rects");
            // check if mod_fill is in legend_rects fill
            var legend_rects_filtered = legend_rects_.filter(function(rect, i) {
                var opacity_rect = d3.select(legend_rects_.nodes()[i]).style("opacity");
                var rect_fill = d3.select(legend_rects_.nodes()[i]).style("fill");
                return (rect_fill == mod_fill) && (opacity_rect != opacity_highlight_legend); //if selected and not temporarily selected (determined by highligt legend opacity)
            });
            // if legend_rects_filtered is empty, then the mod type is not selected in legend
            if (legend_rects_filtered.empty()) {
                d3.select(this)
                    .style("fill", "grey")
                    .attr("opacity", opacity);
            }
        }
        // make all legend rects grey if they have opacity = opacity_highlight_legend
        var legend_rects_ = d3.selectAll(".legend_rects");
        var legend_rects_length = legend_rects_.size();
        for (var i = 0; i < legend_rects_length; i++) {
            var legend_rects_opacity = legend_rects_.nodes()[i].style.opacity;
            if (legend_rects_opacity == opacity_highlight_legend) {
                legend_rects_.nodes()[i].style.fill = "grey";
                legend_rects_.nodes()[i].style.opacity = opacity_mod;
                // set attr fill grey
                d3.select(legend_rects_.nodes()[i])
                    .style("fill", "grey")
                    .attr("opacity", opacity_mod);
            }           
        }
    }
           

    var mod_rects = svg.selectAll("boo")
        .data(mod_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*values_distance)
        .attr("y", d=> (histogram_height + selector_height + margin_overview.bottom+12) + d[1]*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        // .attr("fill", d=> d[4])
        .attr("fill", 'grey')
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity_mod)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)
        .on("mouseenter",mouseover_mod)
        .on("mouseleave", mouseleave_mod);


    // Inspired by: https://www.d3-graph-gallery.com/graph/custom_legend.html
    //   modification labels
    var mod_label_size = 10
    var legend = svg.selectAll("legend")
        .data(modification_color_map_keys)
        .enter()
        .append("rect")
        .attr("class", "legend_rects")
        .attr("x", 0)
        .attr("y", function(d,i){ return 100 + i*(mod_label_size+5)}) // 100 is where the first dot appears. 25 is the distance between dots
        .attr("width", mod_label_size)
        .attr("height", mod_label_size)
        // .style("fill", function(d, i){ return modification_color_map[d]})
        .style("fill",'grey')
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
            color_val = modification_color_map[d];
            if(selectedMods.has(d)){
                selectedMods.delete(d);
            } else {
                selectedMods.set(d, color_val);
            }
            toggleColorModRects(mod_rects, color_val, opacity_mod);
            // fill box with color of text next to it (if it is grey)
            var box_fill = d3.select(this).style("fill");
            if (box_fill == "rgb(128, 128, 128)" || box_fill == "grey") {
                d3.select(this).style("fill", modification_color_map[d]);
            } else {
                d3.select(this).style("fill", "grey");
            }
        });

    
    // Add one dot in the legend for each name.
    var legend_text = svg.selectAll("myLabels")
        .data(modification_color_map_keys)
        .enter()
        .append("text")
        .attr("x", 0 + mod_label_size*1.2)
        .attr("y", function(d,i){ return 100 + i*(mod_label_size+5) + (mod_label_size/2)}) // 100 is where the first dot appears. 25 is the distance between dots
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
            legend.filter(function(dd, ii){// call legend.click on i'th element
                if (ii == i) {
                    d3.select(this).dispatch("click");
                }
            })
        })

    // move legend and legend_text  to Center LEFT 
    legend.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");
    legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");

    if(rect_patches.length > 1) {

    // find difrence in magnintude and use as steps
    var max_exp = expo(max_intensity,1);
    var min_exp = expo(min_intensity,1);
    var max_exp_last_char = max_exp.substring(max_exp.length-1);
    var min_exp_last_char = min_exp.substring(min_exp.length-1);
    var steps = Math.abs(max_exp_last_char - min_exp_last_char)+1;

    var step = (Math.log(max_intensity) - Math.log(min_intensity)) / (steps - 1);
    var log_steps = [];
    for (var i = 0; i < steps; i++) {
        log_steps.push(Math.exp(Math.log(min_intensity) + i * step));
    }

    //parse log_steps to string
    var log_steps_string = [];
    for (var i = 0; i < log_steps.length; i++) {
        log_steps_string.push(expo(parseFloat(log_steps[i]).toPrecision(2), 2));
    }

    var log_steps_string = log_steps_string.slice(1, log_steps_string.length - 1);
    
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
    var color_legend_text = svg.selectAll("colorLabels")
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

    // make color_legend text label on top of color bar
    var color_legend_text_label = svg.append("text")
        .attr("x", color_legend_x - color_bar_width)
        .attr("y", -20)
        .text("Intensity")
        .attr("text-anchor", "right")
        .style("alignment-baseline", "middle")
        // .style("font-size", "12px")
        .attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");      
    }

    var legend_title = svg.append("text")
    .attr("x", 0)
    // same y as color_legend_text_label
    .attr("y", (height/2 - margin.top + 100) - 15)
    .text("Modification Types")
    .attr("text-anchor", "left")
    .style("text-anchor", "left");

    var sub_segment_height = selector_height/height * segment_height;
    var sub_values_distance = Math.min(values_distance, width/peptide_length);
    var selector_width = Math.min(values_distance*peptide_length, Math.round(parseFloat((width/values_distance*width)/peptide_length)));
    var selector_range = Math.min(width, values_distance*peptide_length);
     // add title to legend and legend_text

    var sub_segments = svg.selectAll("sub_foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*sub_values_distance)
        .attr("y", d=>  histogram_height + d[1]*sub_segment_height)
        .attr("width", d=> d[2]*sub_values_distance)
        .attr("height", d=> d[3]*sub_segment_height)
        .attr("fill", d=> d[4])
        .attr("opacity", 0.5);

    var selector = svg.append("rect")
        .attr("class", "mover")
        .attr("x", 0)
        .attr("y", histogram_height)
        .attr("height", selector_height)
        .attr("width", selector_width)
        .attr("fill", "grey")
        .attr("opacity", 0.5)
        .attr("pointer-events", "all")
        .attr("cursor", "ew-resize")
        .call(d3.drag().on("drag", drag_plt));

    var left_selector_border = svg.append("rect")
        .attr("class", "left_selector_border")
        .attr("x", 0)
        .attr("y", histogram_height)
        .attr("height", selector_height)
        .attr("width", 4)
        .attr("fill", "grey")
        .attr("opacity", 1)
        .attr("pointer-events", "all")
        .attr("cursor", "ew-resize")
        .call(d3.drag().on("drag", zoom_plt_left));
    
    var right_selector_border = svg.append("rect")
        .attr("class", "right_selector_border")
        .attr("x", selector_width)
        .attr("y", histogram_height)
        .attr("height", selector_height)
        .attr("width", 4)
        .attr("fill", "grey")
        .attr("opacity", 1)
        .attr("pointer-events", "all")
        .attr("cursor", "ew-resize")
        .call(d3.drag().on("drag", zoom_plt_right));

    function drag_plt() {
        var selector_width = parseInt(d3.select('.mover').attr('width'));

        var left_selector_border_x = parseInt(d3.select('.left_selector_border').attr('x'));

        var overview_scale = selector_range / selector_width;

        var x = parseInt(d3.select(this).attr("x"));
        var nx = x + d3.event.dx,
            w = parseInt(d3.select(this).attr("width"));

        if(nx < 0 - left_selector_border_x || nx + w > selector_range - left_selector_border_x) return;

        d3.select(this).attr("x", nx);

        rects.attr("transform", "translate(" + (-selector_range/selector_width * nx - overview_scale *left_selector_border_x) + "," + 0 + ")");
        mod_rects.attr("transform", "translate(" + (-selector_range/selector_width * nx - overview_scale *left_selector_border_x) + "," + 0 + ")");
        x_ticks.attr("transform", "translate(" + (-selector_range/selector_width * nx - overview_scale *left_selector_border_x) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")");
        x_labels.attr("transform", "translate(" + (-selector_range/selector_width * nx - overview_scale *left_selector_border_x) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")");
        left_selector_border.attr("transform", "translate(" + nx + "," + 0 + ")");
        right_selector_border.attr("transform", "translate(" + nx + "," + 0 + ")");
    }

    function zoom_plt_right() {
        var selector_x = parseFloat(d3.select('.mover').attr('x')),
            selector_width = parseFloat(d3.select('.mover').attr('width'));
        var left_selector_border_x = parseFloat(d3.select('.left_selector_border').attr('x'));

        var x = parseFloat(d3.select(this).attr("x"));
        var nx = parseFloat(d3.event.x),
            ndx = parseFloat(x + d3.event.dx);

        if(nx < selector_x + left_selector_border_x + 40 || nx > selector_range) return;
        var scaleFactor = parseFloat(selector_width/(ndx - left_selector_border_x));

        d3.select(this).attr("x", ndx)
        selector.attr("width", ndx - left_selector_border_x);

        rects.attr("width", function() {
            return this.getAttribute('width') * scaleFactor;
        });
        rects.attr("x", function() {
            return this.getAttribute('x') * scaleFactor;
        });
           
        mod_rects.attr("width", function() {
            return this.getAttribute('width') * scaleFactor;
        });
        mod_rects.attr("x", function() {
            return this.getAttribute('x') * scaleFactor;
        });

        var new_selector_width = parseFloat(ndx - left_selector_border_x);
        var overview_scale = parseFloat(selector_range / new_selector_width);

        mod_rects.attr("transform", function(){
            return "translate(" +  (scaleFactor * (selector_range/ndx) - (selector_x + left_selector_border_x) *overview_scale)  + "," + 0 + ")";
        });
        rects.attr("transform", function() {
            return "translate(" + (scaleFactor * (selector_range/ndx) - (selector_x + left_selector_border_x) *overview_scale) + "," + 0 + ")";
        });

        var tick_r = xScale_ticks.range();
        xScale_ticks.range([0, tick_r[1] * scaleFactor]);

        var label_r = xScale_labels.range();
        xScale_labels.range([label_r[0]*scaleFactor, label_r[1] * scaleFactor]);

        x_ticks
            .attr("transform", "translate(" + (scaleFactor * (selector_range/ndx) - (selector_x + left_selector_border_x) *overview_scale) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
            .call(xAxis_ticks);

        x_labels
            .attr("transform", "translate(" + (scaleFactor * (selector_range/ndx) - (selector_x + left_selector_border_x) *overview_scale) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
            .call(xAxis_labels);
    }

    function zoom_plt_left() {
        var selector_width = parseFloat(d3.select('.mover').attr('width'));
        var selector_x = parseFloat(d3.select('.mover').attr('x'));
        var right_selector_border_x = parseFloat(d3.select('.right_selector_border').attr('x'));
        
        var x = parseFloat(d3.select(this).attr("x"));
        var nx = parseFloat(d3.event.x),
            ndx = parseFloat(x + d3.event.dx);
        
        if(nx < 0 || ndx > right_selector_border_x - 40) return;
        var scaleFactor = parseFloat(selector_width/(right_selector_border_x - ndx));

        d3.select(this).attr("x", ndx)
        selector.attr("transform", "translate(" + ndx + "," + 0 + ")");
        selector.attr('width', right_selector_border_x - ndx);

        rects.attr("width", function() {
            return this.getAttribute('width') * scaleFactor;
        });

        rects.attr("x", function() {
            return this.getAttribute('x') * scaleFactor;
        });
           
        mod_rects.attr("width", function() {
            return this.getAttribute('width') * scaleFactor;
        });
        mod_rects.attr("x", function() {
            return this.getAttribute('x') * scaleFactor;
        });

        var new_selector_width = parseFloat(right_selector_border_x - ndx);
        var overview_scale = parseFloat(selector_range / new_selector_width);

        mod_rects.attr("transform", function(){
            return "translate(" +  (scaleFactor * (selector_range/right_selector_border_x) - (ndx + selector_x) * overview_scale)  + "," + 0 + ")";
        });
        rects.attr("transform", function() {
            return "translate(" + (scaleFactor * (selector_range/right_selector_border_x)  - (ndx + selector_x) * overview_scale) + "," + 0 + ")";
        });

        var tick_r = xScale_ticks.range();
        xScale_ticks.range([0, tick_r[1] * scaleFactor]);

        var label_r = xScale_labels.range();
        xScale_labels.range([label_r[0]*scaleFactor, label_r[1] * scaleFactor]);

        x_ticks
            .attr("transform", "translate(" + (scaleFactor * (selector_range/right_selector_border_x)  - (ndx + selector_x) * overview_scale) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
            .call(xAxis_ticks);

        x_labels
            .attr("transform", "translate(" + (scaleFactor * (selector_range/right_selector_border_x)  - (ndx + selector_x) * overview_scale) + "," + (histogram_height + selector_height + margin_overview.bottom - 1+12) + ")")
            .call(xAxis_labels);
    }

    // get modification positions
    makeHistogramPlot(histogram_data, selector_height, peptide_length, width, svg);

    });
}

(function() {
    renderSegmentPlot();
    // check if pressed
    window.setInterval(function(){
        if (getTabPressed() == true){
            renderSegmentPlot();
            setDocumentLabels();
        }
      }, 500);     
 })();


function highlightSegmentAndModInXaxis(rect_patches, d, mod_patches, colors_to_mod_map, peptide_seq, svg) {
    // get segment in rect_patches where d[0] and d[1] are located and call highlightSeqInXAxis
    var segment = rect_patches.filter(function (rect) {
        return (rect[0] <= d[0] && d[0] <= rect[0] + rect[2] && rect[1] <= d[1] && d[1] < rect[1] + rect[3]);
    });

    // get all modification on this segments - retuns [mod_types_and_positions, mod_positions]
    var modPositionsAndTypes = getModificationPositions(mod_patches, segment[0], colors_to_mod_map, peptide_seq);
    var mod_positions = modPositionsAndTypes[1];
    var specific_mod_pos = d[0];
    var x_axis_highlight = highlightSeqInXAxis(svg, segment[0], mod_positions, specific_mod_pos, true);
}

function makeModificationTooltip(colors_to_mod_map, d, tooltip) {
    tooltip.style("opacity", 1);
    var mod_type = colors_to_mod_map.get(d[4]);
    tooltip.html("<p>Intensity: " + expo(d[5], 3) + "</p><p>Modification Type: " + mod_type + "</p>")
        .style("left", (d3.event.pageX + 10) + "px")
        .style("top", (d3.event.pageY - 10) + "px");
}

function setXaxisNormal(svg) {
    return svg.selectAll(".xAxis_labels")
        .selectAll("text")
        .style("fill", 'black')
        .style("font-size", "1em")
        .style("font-weight", "normal")
        .style("opacity", 1);
}

function toggleColorModRects(mod_rects, color_val, opacity_mod) {
    selectionColorVal = color_val;
    selectionOpacityVal = opacity_mod;
    mod_rects.filter(function (dd) {
        if (color_val == dd[4]) {
            // check color
            col = d3.select(this).style("fill");
            if (col == "rgb(128, 128, 128)" || col == 'grey') { // if grey
                d3.select(this).style("fill", color_val);
                d3.select(this).attr("opacity", opacity_mod);
            } else {
                d3.select(this).style("fill", "grey");
                d3.select(this).attr("opacity", opacity_mod);
            }
        }
    }); 
}

function makeHistogramPlot(histogram_data, selector_height, peptide_length, width, svg) {
    var modification_positions = Object.values(histogram_data.x);
    var unique_values = Array.from(new Set(modification_positions)).length; // get number of unique modification positions
    var num_bins = peptide_length ; // determine bin size based on number of unique modification positions
    var histogram_height = selector_height / 4;

    // make histogram from values
    var histogram = d3.histogram()
        .domain(d3.extent(modification_positions)) // then the domain of the graphic
        .thresholds(num_bins); // then the numbers of bins
    var bins = histogram(modification_positions);

    // make histogram x and y axis
    var x = d3.scaleLinear()
        .domain([1, peptide_length+1])
        .range([0, width]);
    var y = d3.scaleLinear()
        .domain([0, d3.max(bins, d => d.length)])
        .range([histogram_height, 0]);
    var xAxis = d3.axisBottom(x);

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + 0 + "," + histogram_height + ")")
        .call(xAxis);

    // make histogram bars
    var bar = svg.selectAll(".bar")
        .data(bins)
        .enter().append("g")
        .attr("class", "bar")
        .attr("transform", function (d) { return "translate(" + x(d.x0) + "," + y(d.length) + ")"; });
    bar.append("rect")
        .attr("x", 1)
        .attr('fill', 'grey')
        .style("opacity", 0.7)
        .attr("width", x(bins[0].x1) - x(bins[0].x0) - 1)        
        .attr("height", function (d) { return histogram_height - y(d.length); });

    // delete x-axis ticks and labels
    svg.selectAll(".x.axis .tick").remove();
    svg.selectAll(".x.axis .domain").remove();
}

function getModificationPositions(mod_patches, d, colors_to_mod_map, peptide_sequence) {
    var mod_types_and_positions = [];
    var mod_positions = [];
    for (var i = 0; i < mod_patches.length; i++) {
        if (mod_patches[i][0] >= d[0] && mod_patches[i][0] <= d[0] + d[2] && mod_patches[i][1] == d[1]) {
            var pos = mod_patches[i][0];
            var letter = peptide_sequence[pos]; //letter at position
            mod_positions.push(pos);
            var _type = colors_to_mod_map.get(mod_patches[i][4]);
            mod_types_and_positions.push(_type + "(" + letter + ")");
        }
    }
    // make array of arrays mod_types_and_positions and mod_positions
    return [mod_types_and_positions, mod_positions];
}

function highlightSeqInXAxis(svg, d, mod_positions, specific_mod_pos, has_specific) {
    return svg.selectAll(".xAxis_labels")
        .selectAll("text")
        // make bold
        .style('font-weight', function (dd, i) {
            if (i >= d[0] && i < d[0] + d[2] && mod_positions.includes(i)) {
                return 'bold';
            } else {
                return 'normal';
            }
        })
        // increase size of the text
        .style("font-size", function (dd, i) {
            if (has_specific && i == specific_mod_pos) {
                return "1.5em";
            // }else if 
            } else if (i >= d[0] && i < d[0] + d[2] && mod_positions.includes(i)) {
                return "1.2em";

            } else if (i >= d[0] && i < d[0] + d[2]) {
                return "1.1em";    
            } else {
                return "1em";
            }
        })
        // opacity
        .style("opacity", function (dd, i) {   
            if (i >= d[0] && i < d[0] + d[2]) { // and within the range
                return 1;
            }
            else {
                return 0.4;
            }
        })};


    function expo(x, f) {
        return Number.parseFloat(x).toExponential(f);
    }
