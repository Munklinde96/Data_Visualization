function renderSampleModPlot(){
    if(selectedProtein === ""){
        d3.select("#graphDiv2").select("svg").remove();
        document.getElementById('no_protein_div_2').innerHTML = '<div style="width: 640px; height: 440px;"><h3>Select a protein to get started.</h3></div>';
        return;
    } else {
        document.getElementById('no_protein_div_2').innerHTML = "<div></div>";
    }
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 130, bottom: 30, left: 130},
        width = 900 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;
    
    // remove old svg
    d3.select("#graphDiv2").select("svg").remove();
    
    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv2")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");
    //Read the data
    d3.json("http://127.0.0.1:5000/get-sample-mod-data?protein="+selectedProtein, function(error, data) {
        // validate request
        if (error) throw error;
        // build modification and sample categories
        var maxValue = null;
        var minValue = null;
        var modStructData = null;
        var differentSamples = [];
        for (const [sample, mods] of Object.entries(data)) {
            for (const [mod, value] of Object.entries(mods)){
                if(maxValue == null || value > maxValue){
                    maxValue = value;
                }
                if(minValue == null || value < minValue) {
                    minValue =  value
                }
                if(differentSamples.includes(sample) == false){
                    differentSamples.push(sample);
                }
                if(mod && sample && value){
                    if(modStructData == null){
                        modStructData = [{
                        sample: sample,
                        mod: mod,
                        value: value,
                    }];
                    } else {
                        modStructData.push({
                        sample: sample,
                        mod: mod,
                        value: value,
                    });
                    }
                }
            }
        }
        // Build color scale
        var myColor = d3.scaleSequential()
        .domain([minValue,maxValue])
        .interpolator(d3.interpolateBlues);

        // Labels of row and columns
        var myGroups = d3.map(modStructData, function(d){return d.sample;}).keys()
        var myVars = d3.map(modStructData, function(d){return d.mod;}).keys();
                
        // Build X scales and axis:
        var x = d3.scaleBand()
        .range([ 0, width ])
        .domain(myGroups)
        .padding(0.01);
        svg.append("g")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x))

        // Build Y scales and axis:
        var y = d3.scaleBand()
        .range([ height, 0 ])
        .domain(myVars)
        .padding(0.01);
        svg.append("g").call(d3.axisLeft(y));

        // create a tooltip
        var tooltip = d3.select("#graphDiv2")
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
            tooltip.style("opacity", 1)
        }
        var mousemove = function(d) {
        tooltip
            .html("The exact value of<br>this cell is: " + d.value)
            .style("left", (d3.mouse(this)[0]) + "px")
            .style("top", (d3.mouse(this)[1]) + "px")
        }
        var mouseleave = function(d) {
        tooltip.style("opacity", 0)
        }

        var mouseclick = function(d, i){
            for(let n = 0; n < modStructData.length; n++){
                // remove old if clicked on one it already contains
                if(differentSamples.length === selectedSample.length || (modStructData[n].sample === d.sample && selectedSample.includes(Number(d.sample)))){
                    let target = d3.select("#sample_id_heatmap_"+n);
                    target.style('stroke', 'none');
                } else if(Number(d.sample) === Number(modStructData[n].sample)){
                    let target = d3.select("#sample_id_heatmap_"+n);
                    target.style('stroke', 'red');
                    target.style('stroke-width', 2);
                }
            }
            if(selectedSample.includes(Number(d.sample)) === false){
                selectedSample.push(Number(d.sample));
            } else {
                selectedSample = selectedSample.filter(function(val, index, err){
                    return val != d.sample;
                });
            }
            renderSegmentPlot();
        }

        svg.selectAll()
            .data(modStructData, function(d){ return d.sample+":"+d.mod })
            .enter()
            .append("rect")
            .attr("x", function(d) {return x(d.sample) })
            .attr("y", function(d) { return y(d.mod) })
            .attr("id", function(d, i){return "sample_id_heatmap_"+i })
            .attr("width", x.bandwidth() )
            .attr("height", y.bandwidth() )
            .style("fill", function(d) { return myColor(d.value)} )
            .on("mouseover", mouseover)
            .on("mousemove", mousemove)
            .on("mouseleave", mouseleave)
            .on("click", mouseclick)
    });
}

$('document').ready(function(){
    renderSampleModPlot();
});
