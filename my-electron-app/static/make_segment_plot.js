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
        .attr("stroke-width", 0.2 *scale_factor);

    var mod_rects = svg.selectAll("boo")
        .data(mod_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*scale_factor)
        .attr("y", d=> d[1]*scale_factor)
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
        .attr("x", 100)
        .attr("y", function(d,i){ console.log(d); return 100 + i*(size+5)}) // 100 is where the first dot appears. 25 is the distance between dots
        .attr("width", size)
        .attr("height", size)
        .style("fill", function(d, i){ console.log(d) ;return modification_color_map[d]})
    
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
    


    });
});
