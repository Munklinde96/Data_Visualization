//start file


$('document').ready(function(){
    d3.json("http://127.0.0.1:5000/get-segment-data", function(error, data) {
        if (error) throw error;
        console.log(data);
        var rect_patches = data.peptide_patches;
        var mod_patches = data.mod_patches;
        var plot_height = data.height;
        var peptide_seq = data.seqq;
        var modification_color_map = data.modification_color_map;
        var min_peptide = data.min_peptide;
        var max_peptide = data.max_peptide;
        var min_intensity = min_peptide[0];
        var max_intensity = max_peptide[0];
        var min_color = min_peptide[1];
        var max_color = max_peptide[1];
        var color = d3.scaleLinear().range([min_color, max_color]).domain([1, 2, 3, 4, 5]);


    
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 130, bottom: 30, left: 130},
        width = screen.width - margin.left - margin.right,
        height = 700 - margin.top - margin.bottom;
    
        
    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv3")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
    // // add box around everything without margin
    // svg.append("rect")
    //     .attr("width", width)
    //     .attr("height", height)
    //     .attr("fill", "none")
    //     .attr("stroke", "black")
    //     .attr("stroke-width", 1);
    // Add X axis with labels from top
    //make random list of chars of lenth 100
    // var random_charlist= [];
    // for (var i = 0; i < 100; i++) {
    //     random_charlist.push(String.fromCharCode(Math.floor(Math.random() * 26) + 97));
    // }
    // make scale factor based on window size
    var scale_factor = width/peptide_seq.length;

    // var scale_factor = (width - margin.left - margin.right)/peptide_length;

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;
    var peptide_chars_with_numbers = peptide_chars.map(function(d,i){return d+i;});
    
    var xScale = d3.scaleLinear()
        .domain([0, peptide_length])
        .range([0, width]);

    // set ticks to be chars of peptide seqq
    // based on width of window calculate tick_values_distance
    var tick_values_distance = Math.floor(width/peptide_length);
    var xAxis = d3.axisTop(xScale)
        .tickValues(d3.range(0, peptide_length, tick_values_distance))
        .tickFormat(function(d,i){console.log(d); return peptide_chars[d];});
        // TODO: adjust when zooming

    svg.append("g")
        .attr("transform", "translate(" + 0 + ")")
        .call(xAxis);
        
        // var xAxis = d3.axisTop(xScale)
        //     .tickFormat(function(d) {return d[0];});
    // var xAxis = svg.append("g")
    //     .attr("transform", "translate(" + margin.left + ")")
    //     .call(d3.axisTop(xScale).tickSize(0))
    //     .selectAll("text")
        
    //calculate the scale factor to multiple size of each rectangle
    var scale_factor = (width - margin.left - margin.right)/peptide_length;

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
        .attr("stroke", "black")
        .attr("stroke-width", 0.2 *scale_factor)
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
        tooltip.html("<p>Peptide: " + d[0] + "</p><p>Intensity: " + d[1] + "</p><p>Modtype: " + d[3] + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
    }
    function mousemove_segments(d) {
        tooltip.html("<p>Peptide: " + d[0] + "</p><p>Intensity: " + d[1] + "</p>")
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
    }
    
    function mouseleave(d) {
        tooltip.style("opacity", 0)
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
        .attr("stroke", "black")
        .attr("stroke-width", 0.2 *scale_factor)
        .on("mouseover",mouseover)
        .on("mousemove", mousemove_modification)
        .on("mouseout", mouseleave);

    // Inspired by: https://www.d3-graph-gallery.com/graph/custom_legend.html
    keys = Object.keys(modification_color_map);
    var size = 10
    var legend = svg.selectAll("legend")
        .data(keys)
        .enter()
        .append("rect")
        .attr("x", 100)
        .attr("y", function(d,i){ return 100 + i*(size+5)}) // 100 is where the first dot appears. 25 is the distance between dots
        .attr("width", size)
        .attr("height", size)
        .style("fill", function(d, i){ return modification_color_map[d]})
    
    // Add one dot in the legend for each name.
    legend_text = svg.selectAll("myLabels")
        .data(keys)
        .enter()
        .append("text")
        .attr("x", 100 + size*1.2)
        .attr("y", function(d,i){ return 100 + i*(size+5) + (size/2)}) // 100 is where the first dot appears. 25 is the distance between dots
        .text(function(d){ return d})
        .attr("text-anchor", "left")
        .style("alignment-baseline", "middle")
        
       


    // move legend and legend_text  to Center LEFT 
    legend.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");
    legend_text.attr("transform", "translate(" + 0+ "," + (height/2 - margin.top) + ")");
    
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
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", 20)
        .attr("height", 170)
        .style("fill", "url(#linear-gradient)");

    rect.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    //add text to top and bottom of rect
    svg.append("text")
        .attr("x", 0)
        .attr("y", 0)
        .attr("text-anchor", "middle")
        .style("alignment-baseline", "ideographic")
        .text(String(max_intensity))
        .attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");

    svg.append("text")
        .attr("x", 0)
        .attr("y", 0)
        .attr("text-anchor", "middle")
        .style("alignment-baseline", "hanging")
        .text(String(min_intensity))
        .attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100 + 170) + ")");
    });
});
