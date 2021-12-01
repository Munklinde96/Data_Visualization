function renderSampleModPlot(){
    if(getSelectedProtein === ""){
        d3.select("#graphDiv2").select("svg").remove();
        document.getElementById('no_protein_div_2').innerHTML = '<div style="width: 640px; height: 440px;"><h3>Select a protein to get started.</h3></div>';
        return;
    } else {
        document.getElementById('no_protein_div_2').innerHTML = "<div></div>";
    }
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 130, bottom: 30, left: 130},
        width = screen.width*0.35,
        height = 500 - margin.top - margin.bottom;
    
    // remove old svg
    d3.select("#graphDiv2").select("svg").remove();
    
    // append the svg object to the body of the page
    var svg = d3.select("#graphDiv2")
        .append("svg")
        .attr("width", width + margin.left + 2*margin.right)
        .attr("height", height + margin.top + 2* margin.bottom)
        .append("g")
        .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");
    //Read the data
    d3.json("http://127.0.0.1:5000/get-sample-mod-data?protein="+getSelectedProtein(), function(error, data) {
        // validate request
        if (error) throw error;
        // build modification and sample categories
        var maxValue = 0;
        var minValue = 0;
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
        .interpolator(d3.interpolateOranges);

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
        .call(d3.axisBottom(x));

        // x-axis label in middle of x-axis
        svg.append("text")
        .attr("text-anchor", "middle")
        .attr("transform", "translate("+ (width/2) +","+(height+margin.top+20)+")")
        .text("Sample");

        // Build Y scales and axis:
        var y = d3.scaleBand()
        .range([ height, 0 ])
        .domain(myVars)
        .padding(0.01);
        svg.append("g").call(d3.axisLeft(y));

        // add label to y axis
        svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0 - margin.left)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Modification Type");


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
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        }
        var mouseleave = function(d) {
        tooltip.style("opacity", 0)
        }

        var mouseclick = function(d, i){
            for(let n = 0; n < modStructData.length; n++){
                // remove old if clicked on one it already contains
                if(modStructData[n].sample === d.sample && getSelectedSamples().includes(Number(d.sample))){
                    let target = d3.select("#sample_id_heatmap_"+n);
                    target.style('stroke', 'none');
                } else if(Number(d.sample) === Number(modStructData[n].sample)){
                    let target = d3.select("#sample_id_heatmap_"+n);
                    target.style('stroke', 'red');
                    target.style('stroke-width', 2);
                }
            }
            if(getSelectedSamples().includes(Number(d.sample)) === false){
                var tempSamples = getSelectedSamples();
                tempSamples.push(Number(d.sample));
                setSelectedSamples(tempSamples);
            } else {
                setSelectedSamples(getSelectedSamples().filter(function(val, index, err){
                    return val != d.sample;
                }));
            }
            setDocumentLabels();
            //renderSegmentPlot();
        }

        svg.selectAll()
            .data(modStructData, function(d){ return d.sample+":"+d.mod })
            .enter()
            .append("rect")
            .attr("x", function(d) {return x(d.sample) })
            .attr("y", function(d) { return y(d.mod) })
            .attr("rx",3)
            .attr("ry",3)
            .attr("id", function(d, i){return "sample_id_heatmap_"+i })
            .attr("width", x.bandwidth() )
            .attr("height", y.bandwidth() )
            .style("fill", function(d) { return myColor(d.value)} )
            .on("mouseover", mouseover)
            .on("mousemove", mousemove)
            .on("mouseleave", mouseleave)
            .on("click", mouseclick)

        // COLORBAR LEGEND
        color_bar_width = 30;
        color_bar_height = height;
    
        var defs = svg.append("defs");
        var linearGradient = defs.append("linearGradient")
        .attr("id", "linear-gradient_sample_mod_heat")
        .attr("x1", "0%")
        .attr("y1", "100%")
        .attr("x2", "0%")
        .attr("y2", "0%")
        .attr("spreadMethod", "pad");


        // get 100 points from my color and add to linearGradient
        var gradinet_steps = 100;
        for(let i = 0; i < gradinet_steps ; i++){
            value = minValue + (maxValue - minValue) * i / gradinet_steps ;
            color = myColor(value);
            offset = (i / gradinet_steps) * 100 + "%";
            linearGradient.append("stop")
            .attr("offset", offset)
            .attr("stop-color", color)
        }

        var rect = svg.append("rect")
        .attr("x", width +5)
        .attr("y", 0)
        .attr("width", color_bar_width)
        .attr("height", color_bar_height)
        .style("fill", "url(#linear-gradient_sample_mod_heat)");

        // legend axis and scale
        // create a scale and axis for the legend
        var legendScale = d3.scaleLinear()
        .domain([maxValue, minValue])
        .range([0,color_bar_height]);
        
        var legendAxis = d3.axisRight(legendScale)
        .ticks(5)
        .tickSize(5)
        .tickFormat(d3.format(".2f"));
        // set tick for last value also
        legendAxis.tickValues(legendAxis.scale().ticks(5).concat(legendAxis.scale().domain()));

        // add the legend axis to the svg
        svg.append("g")
        .attr("transform", "translate(" + (width + color_bar_width + 5) + ",0)")
        .call(legendAxis);
 
        // add text label on top of colorbar
        svg.append("text")
        .attr("text-anchor", "middle")
        .attr("transform", "translate("+ (width + color_bar_width + 35 + color_bar_width/2) +","+(color_bar_height/2)+")rotate(-90)")
        .text("Modification Frequency");

        
    });
}

(function() {
    renderSampleModPlot();
 })();