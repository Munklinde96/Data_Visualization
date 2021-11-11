//start file


function init() {
    var rect_patches = null;
    var mod_patches = null;
    var plot_height = null;
    var peptide_seq = null;


    d3.json("http://127.0.0.1:5000/get-segment-data", function(error, data) {
        if (error) throw error;
        rect_patches = data.peptide_patches;
        mod_patches = data.mod_patches;
        plot_height = data.height;
        peptide_seq = data.peptide_seq;
    });
    
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 130, bottom: 30, left: 130},
        width = 700 - margin.left - margin.right,
        height = plot_height - margin.top - margin.bottom;
    
    // add box around everything without margin
    svg.append("rect")
        .attr("width", width)
        .attr("height", height)
        .attr("fill", "none")
        .attr("stroke", "black")
        .attr("stroke-width", 1);

    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv2")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Add X axis with labels from top
    var peptide_chars = peptide_seq.split("");
    var peptide_length = peptide_seq.length;
    var x = d3.scaleBand()
        .domain(peptide_chars)
        .range([0, peptide_length]); // find matching length
    svg.append("g")
        .attr("transform", "translate(0," + margin.top + ")")
        .call(d3.axisTop(x));
        
    // sample data for rectangles
    var data = [ 
                {x: 20, y: 60, width: 30, hight: 50},
                {x: 50, y: 80, width: 100, hight: 150},
                {x: 200, y: 400, width: 10, hight: 100}];

    var rects = svg.selectAll("foo")
        .data(data)
        .enter()
        .append("rect")
        .attr("x", d=> d.x)
        .attr("y", d=> d.y)
        .attr("width", d=> d.width)
        .attr("height", d=> d.hight)
        .attr("stroke", "black")
        .attr("stroke-width", 2)
        .attr("fill", "teal");
}
