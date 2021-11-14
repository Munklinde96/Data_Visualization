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

    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;
    var x = d3.scaleBand()
        .domain(peptide_chars)
        .range([0, peptide_length]); // find matching length
    svg.append("g")
        .attr("transform", "translate(0," + margin.top + ")")
        .call(d3.axisTop(x));
        
    // sample data for rectangles

    var rects = svg.selectAll("foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0])
        .attr("y", d=> d[1])
        .attr("width", d=> d[2])
        .attr("height", d=> d[3])
        .attr("fill", d=> d[4])
        .attr("stroke", "black")
        .attr("stroke-width", 2)
        .attr("fill", "teal");

    // var mod_rects = svg.selectAll("boo")
    //     .data(mod_patches)
    //     .enter()
    //     .append("rect")
    //     .attr("x", d => d[0])
    //     .attr("y", d=> d[1])
    //     .attr("width", d=> d[2])
    //     .attr("height", d=> d[3])
    //     .attr("fill", d=> d[4])
    //     .attr("stroke", "black")
    //     .attr("stroke-width", 2)
    //     .attr("fill", "red");
});
});
