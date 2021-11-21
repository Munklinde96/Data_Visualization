function renderProteinModPlot(){
    console.log("drawing mod heatap")
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 130, bottom: 30, left: 130},
        width = 900 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;
    // remove old svg
    d3.select("#mod_heatmap").select("svg").remove();
    
    // append the svg object to the body of the page
    var svg = d3.select("#mod_heatmap")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    //Read the data
    d3.json("http://127.0.0.1:5000/get-protein-mod-data?min_mod_count="+minModificationCount+"&normalization_type="+normalizationType+"&samples="+selectedSample, function(error, data) {
        // validate request
        if (error) throw error;
        // build modification and protein categories
        var maxValue = null;
        var minValue = 0;
        var proteinStructData = null;
        for (const [protein, mods] of Object.entries(data)) {
            for (const [mod, value] of Object.entries(mods)){
                if(maxValue == null || value > maxValue){
                    maxValue = value;
                }
                if(minValue == null || value < minValue) {
                    minValue =  value
                }
                if(mod && protein && value){
                    if(proteinStructData == null){
                        proteinStructData = [{
                        protein: protein,
                        mod: mod,
                        value: value,
                    }];
                    } else {
                        proteinStructData.push({
                        protein: protein,
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
        .interpolator(d3.interpolateViridis);

        // Labels of row and columns
        var myGroups = d3.map(proteinStructData, function(d){return d.protein;}).keys()
        var myVars = d3.map(proteinStructData, function(d){return d.mod;}).keys();
                
        // Build X scales and axis:
        var x = d3.scaleBand()
        .range([ 0, width ])
        .domain(myGroups)
        svg.append("g")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x))
        .selectAll("text")  
        .style("text-anchor", "end")
        .attr("dx", "-.8em")
        .attr("dy", ".15em")
        .attr("transform", "rotate(-45)");    

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
            for(let n = 0; n < proteinStructData.length; n++){
                // remove old
                if(selectedProtein == proteinStructData[n].protein){
                    let target = d3.select("#mod_prot_id_heatmap_"+n);
                    target.style('stroke', 'none');
                }
                if(d.protein === proteinStructData[n].protein && d.protein !== selectedProtein){
                    let target = d3.select("#mod_prot_id_heatmap_"+n);
                    target.style('stroke', 'red');
                    target.style('stroke-width', 2);
                }
            }
            if(selectedProtein !== d.protein){
                selectedProtein = d.protein;
            } else {
                selectedProtein = "";
            }
            renderSampleModPlot();
            renderSegmentPlot();
            selectedSample = [];
        }

        svg.selectAll()
            .data(proteinStructData, function(d){ return d.protein+":"+d.mod })
            .enter()
            .append("rect")
            .attr("x", function(d) {return x(d.protein) })
            .attr("y", function(d) { return y(d.mod) })
            .attr("rx",3)
            .attr("ry",3)
            .attr("width", x.bandwidth())
            .attr("height", y.bandwidth())
            .attr("id", function(d, i){return "mod_prot_id_heatmap_"+i })
            .style("fill", function(d) {return myColor(d.value)})
            .style("stroke", function(d){
                if(d.protein === selectedProtein){
                    return "red";
                }
                return "none";
            })
            .style("stroke-width", function(d){
                if(d.protein === selectedProtein){
                    return 2;
                }
                return 0;
            })
            .on("mouseover", mouseover)
            .on("mousemove", mousemove)
            .on("mouseleave", mouseleave)
            .on("click", mouseclick)

        // COLORBAR LEGEND

        // make sequential gradient 
        color_bar_width = 30;
        color_bar_height = height;
    
        var defs = svg.append("defs");
        var linearGradient = defs.append("linearGradient")
        .attr("id", "linear-gradient_mod_heat")
        .attr("x1", "0%")
        .attr("y1", "100%")
        .attr("x2", "0%")
        .attr("y2", "0%")
        .attr("spreadMethod", "pad");

        // get 100 points from my color and add to linearGradient
        for(let i = 0; i < 100; i++){
            value = minValue + (maxValue - minValue) * i / 100;
            color = myColor(value);
            linearGradient.append("stop")
            .attr("offset", i+"%")
            .attr("stop-color", color)
        }

        var rect = svg.append("rect")
        .attr("x", width +5)
        .attr("y", 0)
        .attr("width", color_bar_width)
        .attr("height", color_bar_height)
        .style("fill", "url(#linear-gradient_mod_heat)");

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
        
    });
}

$('document').ready(function(){
    renderProteinModPlot();
});
