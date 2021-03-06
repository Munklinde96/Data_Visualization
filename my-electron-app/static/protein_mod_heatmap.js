
function renderProteinModPlot(){
    // set the dimensions and margins of the graph
    var margin = {top: 30, right: 65, bottom: 30, left: 130},
        width = screen.width*0.325,
        height = 500 - margin.top - margin.bottom;
    // remove old svg
    d3.select("#mod_heatmap").select("svg").remove();
    
    // append the svg object to the body of the page
    var svg = d3.select("#mod_heatmap")
        .append("svg")
        .attr("width", width + margin.left + 2*margin.right)
        .attr("height", height + margin.top + 2*margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    //Read the data
    d3.json("http://127.0.0.1:5000/get-protein-mod-data?min_mod_count="+getMinModificationCount()+"&normalization_type="+getSelectedNormalization()+"&samples="+getSelectedSamples(), function(error, data) {
        // validate request
        if (error) throw error;
        // build modification and protein categories
        var maxValue = null;
        var minValue = 0;
        var proteinStructData = null;
        var sortingOrderProteins = getSortingOrderProteins();
        var normalizationType = getSelectedNormalization();

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
       
        var start = d3.hsl(210, 1, 0.90); // org color min 225??, 100%, 70%
        var end = d3.hsl(210, 0.3, 0.20);   // org color max 225??, 30%, 20%
        var myColor = d3.scaleSequential(d3.interpolateHsl(start, end))
            .domain([minValue, maxValue]);

        // Protein categories
        var proteins = d3.map(proteinStructData, function(d){return d.protein;}).keys()
        var initial_proteins = proteins;
        if ( Object.keys(sortingOrderProteins).length === 0) {
            sortProteinsBySumOfAllModifications(proteins, proteinStructData);
            sortingOrderProteins = proteins; 
            setSortingOrderProteins(proteins);
        } else {
            proteins = sortingOrderProteins
            for (var i = 0; i < proteins.length; i++) { // remove proteins that are not in initial_proteins 
                if (initial_proteins.indexOf(proteins[i]) == -1) {
                    proteins.splice(i, 1);
                    i--;
                }
            }
        }
        // Modification categories
        var modifications = d3.map(proteinStructData, function(d){return d.mod;}).keys();
        sortModificationsAndMoveUnmofifiedToTop(modifications);
        
        // Build X scales and axis:
        var x = d3.scaleBand()
        .range([ 0, width ])
        .domain(proteins)
        svg.append("g")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x))
        .selectAll("text")  
        .style("text-anchor", "end")
        .attr("transform", "rotate(-90)")
        .attr("dx", "-.8em")
        .attr("dy", "-0.6em");
        
        // Build Y scales and axis:
        var y = d3.scaleBand()
        .range([ height, 0 ])
        .domain(modifications)
        svg.append("g")
        .call(d3.axisLeft(y));

        // add the text label for the x axis
        svg.append("text")
        .attr("transform",
            "translate(" + (width/2) + " ," +
                            (height + margin.top + 30) + ")")
        .style("text-anchor", "middle")
        .text("Protein");

         // add label to y axis
         svg.append("text")
         .attr("transform", "rotate(-90)")
         .attr("y", 0 - margin.left)
         .attr("x",0 - (height / 2))
         .attr("dy", "1em")
         .style("text-anchor", "middle")
         .text("Modification Type");

        // add the subtitle
        var normalizationType = getSelectedNormalization();
        
        // Build Y scales and axis:
        var y = d3.scaleBand()
        .range([ height, 0 ])
        .domain(modifications)
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
        tooltip.html( d.value)
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY - 10) + "px");
        }
        var mouseleave = function(d) {
        tooltip.style("opacity", 0)
        }

        var mouseclick = function(d, i){
            for(let n = 0; n < proteinStructData.length; n++){
                // remove old
                if(getSelectedProtein() == proteinStructData[n].protein){
                    let target = d3.select("#mod_prot_id_heatmap_"+n);
                    target.style('stroke', 'none');
                }
                if(d.protein === proteinStructData[n].protein){
                    let target = d3.select("#mod_prot_id_heatmap_"+n);
                    target.style('stroke', 'red');
                    target.attr('rx', 2)
                    target.attr('ry', 2)
                    target.style('stroke-width', 2);
                }
            }
            setSelectedProtein(d.protein);            
            setSelectedSamples([]);
            renderSampleModPlot();
            setDocumentLabels();
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
                if(d.protein === getSelectedProtein()){
                    // egg white
                    return 'red';
                }
                return "none";
            })
            .style("stroke-width", function(d){
                if(d.protein === getSelectedProtein()){
                    return 2;
                }
                return 0;
            })
            .on("mouseover", mouseover)
            .on("mousemove", mousemove)
            .on("mouseleave", mouseleave)
            .on("click", mouseclick)

        // make Black line and small spacing above row for 'Unmodified'
        svg.append("line")
        .attr("x1", 0)
        .attr("y1", y("Unmodified"))
        .attr("x2", width)
        .attr("y2", y("Unmodified"))
        .style("stroke", "black")
        .style("stroke-width", 4);

        // make 'Unmodifed' tick label in y axis italic
        makeUnmodifiedLabelItalic(svg);

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
        .tickFormat( function(d){
            // if normalization is protein_intensity, use scientific notation
            if(getSelectedNormalization() === "protein_intensity"){
                return d3.format(".2e")(d);
            }
            else return d.toFixed(2);
            })
        // set tick for last value also
        legendAxis.tickValues(legendAxis.scale().ticks(5).concat(legendAxis.scale().domain()));

        // add the legend axis to the svg
        svg.append("g")
        .attr("transform", "translate(" + (width + color_bar_width + 5) + ",0)")
        .call(legendAxis);
        
        // add text label to colorbar
        svg.append("text")
        .attr("transform", "translate("+ (width + color_bar_width + 35 + color_bar_width) +","+(color_bar_height/2)+")rotate(-90)")
        .style("text-anchor", "middle")
        .text(function(){
            return "Modification Count";
        });
        
    });
}


function makeUnmodifiedLabelItalic(svg) {
    svg.selectAll(".tick text")
        .filter(function (d, i) {
            return d === "Unmodified";
        })
        .style("font-style", "italic");
}

function sortModificationsAndMoveUnmofifiedToTop(modifications) {
    modifications.sort(); // sort modifications alphabetically
    var unmodifiedIndex = modifications.indexOf('Unmodified'); // move 'Unmodified' to the top
    if (unmodifiedIndex != -1) {
        modifications.splice(unmodifiedIndex, 1);
        modifications.unshift('Unmodified');
    }
}

function sortProteinsBySpecificOrder(proteinStructData, order) {
    var sortedProteinStructData = [];
    for (var i = 0; i < order.length; i++) {
        for (var j = 0; j < proteinStructData.length; j++) {
            if (order[i] === proteinStructData[j].protein) {
                sortedProteinStructData.push(proteinStructData[j]);
            }
        }
    }
    return sortedProteinStructData;
}

function sortProteinsBySumOfAllModifications(proteins, proteinStructData) {
    proteins.sort(function (a, b) {
        var aSum = 0;
        var bSum = 0;
        for (var i = 0; i < proteinStructData.length; i++) {
            if (proteinStructData[i].protein == a) {
                aSum += proteinStructData[i].value;
            }
            if (proteinStructData[i].protein == b) {
                bSum += proteinStructData[i].value;
            }
        }
        return bSum - aSum;
    });
}


 renderProteinModPlot();
