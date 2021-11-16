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
        //round min_intensity and max_intensity to 10 decimals
        min_intensity = Math.round(min_intensity * 100000000) / 100000000;
        max_intensity = Math.round(max_intensity * 100000000) / 100000000;

        var min_color = min_peptide[1];
        var max_color = max_peptide[1];
        var color = d3.scaleLinear().range([min_color, max_color]).domain([1, 2]);
        
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
    var scale_factor = (width - margin.left - margin.right)/peptide_length;

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;
    var peptide_chars_with_numbers = peptide_chars.map(function(d,i){return d+i;});
    
    var xScale = d3.scaleBand()
        .domain(peptide_chars_with_numbers.map(function(d) {return [d];}))
        .range([0, peptide_length*scale_factor]); // find matching length

    var xAxis = d3.axisTop()
        .scale(xScale)
        .tickFormat(function (d) {
            return peptide_chars_with_numbers[d];
        });

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
        .on("mouseover", mouseover)
        .on("mousemove", mousemove)
        .on("mouseleave", mouseleave);

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
    var mouseover = function(d) {
        console.log("mouse over")
        tooltip.style("opacity", 1)
    }
    var mousemove = function(d) {
        tooltip
        .html("The exact value of<br>this cell is: " + d[3])
        .style("left", (d3.mouse(d)[0]*scale_factor) + "px")
        .style("top", (d3.mouse(d)[1]*scale_factor) + "px")
    }
    var mouseleave = function(d) {
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
        .attr("stroke-width", 0.2 *scale_factor);

    // Inspired by: https://www.d3-graph-gallery.com/graph/custom_legend.html
    keys = Object.keys(modification_color_map);
    var size = 10
    var legend = svg.selectAll("legend")
        .data(keys)
        .enter()
        .append("rect")
        .attr("x", 0)
        .attr("y", function(d,i){ return 100 + i*(size+5)}) // 100 is where the first dot appears. 25 is the distance between dots
        .attr("width", size)
        .attr("height", size)
        .style("fill", function(d, i){ return modification_color_map[d]})
    
    // Add one dot in the legend for each name.
    legend_text = svg.selectAll("myLabels")
        .data(keys)
        .enter()
        .append("text")
        .attr("x", 0 + size*1.2)
        .attr("y", function(d,i){ return 100 + i*(size+5) + (size/2)}) // 100 is where the first dot appears. 25 is the distance between dots
        .text(function(d){ return d})
        .attr("text-anchor", "left")
        .style("alignment-baseline", "middle")

    // move legend and legend_text  to Center LEFT 
    legend.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");
    legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top) + ")");

    color_bar_width = 20;
    color_bar_height = 180;

    var steps = 5;
    var step = (Math.log(max_intensity) - Math.log(min_intensity)) / (steps - 1);
    log_steps = [];
    for (var i = 0; i < steps; i++) {
        log_steps.push(Math.exp(Math.log(min_intensity) + i * step));
    }
    
    //floor log_steps to nearest power of 10
    var log_steps_floor = [];
    for (var i = 0; i < log_steps.length; i++) {
        log_steps_floor.push(Math.pow(10, Math.floor(Math.log(log_steps[i]) / Math.log(10))));
    }

    //parse log_steps to string
    var log_steps_string = [];
    for (var i = 0; i < log_steps.length; i++) {
        log_steps_string.push(parseFloat(log_steps_floor[i]).toPrecision(1));
    }

    //add max_intensity to log_steps_floor
    log_steps_string.push(max_intensity);
    log_steps_string.unshift(min_intensity);
    


    
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
        .attr("y", function(d,i){ return color_bar_height - i*(color_bar_height/6)})
        .text(function(d){ return "-" +d})
        .attr("text-anchor", "right")
        .style("alignment-baseline", "middle")

        //move color_legend_text next to rect
    color_legend_text.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    rect.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");

    });
});
