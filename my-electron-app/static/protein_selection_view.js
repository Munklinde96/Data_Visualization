function removeProteinSelectionPlot(){
    d3.select("#peptide_selection_view").select("svg").remove();
    return;
}

function renderProteinSelectionPlot(){
    // remove old svg
    removeProteinSelectionPlot();
    
    //var p = document.getElementById("peptide_selection_view").getElementsByTagName("p")[0];
    //p.innerHTML = "Peptide Segments Plot - Protein: " + selectedProtein;
    //p.style.fontWeight = "bold";
    //p.style.fontStyle = "italic";

    
    d3.json("http://127.0.0.1:5000/get-segment-protein-data?protein="+getSelectedProtein()+"&samples="+getSelectedSamples()+"&start_pos="+proteinStartPos+"&end_pos="+proteinEndPos, function(error, data) {
        if (error) throw error;
        console.log("data is here");
        console.log(data);
        var rect_patches = data.peptide_patches;
        var mod_patches = data.mod_patches;
        var peptide_seq_prot = selectedSequence;
        console.log(peptide_seq_prot);
        var modification_color_map = data.modification_color_map;
        var min_peptide = data.min_peptide;
        var max_peptide = data.max_peptide;
        var min_intensity = min_peptide[0];
        var max_intensity = max_peptide[0];
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

    var segment_height = 6;
    var values_distance = 10;
    // set the dimensions and margins of the graph
    console.log("test 1");
    console.log("width is: "+data.peptide_patches[0][2]*10);
    console.log("test 2");
    var margin = {top: 30, right: 30, bottom: 30, left: 30},
        width = data.peptide_patches[0][2]*10,
        height = data.height*5;
        
    // append the svg object to the body of the page
    var svg = d3.select("#peptide_selection_view")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    // make scrolable svg
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    


    var rx = 3
    var ry = 3
    var stroke_width = 0.1
    var opacity = 0.7
    var opacity_mod = 0.8
    var opacity_mod_grey = 0.5
    // sample data for rectangles
    var rects = svg.selectAll("foo")
        .data(rect_patches)
        .enter()
        .append("rect")
        //.attr("x", d => d[0]*values_distance)
        .attr("x", d => 0)
        .attr("y", d=> d[1]*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        .attr("fill", d=> d[4])
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)
        
    
    function expo(x, f) {
        return Number.parseFloat(x).toExponential(f);
      }
    
    var subWidth = 0;
    if(rect_patches.length > 0){
        console.log(rect_patches[0][0])
        console.log(values_distance);
        subWidth = rect_patches[0][0]*values_distance;
        console.log(subWidth);
    }

    var mod_rects = svg.selectAll("boo")
        .data(mod_patches)
        .enter()
        .append("rect")
        .attr("x", d => d[0]*values_distance-subWidth)
        .attr("y", d=> (d[1])*segment_height)
        .attr("width", d=> d[2]*values_distance)
        .attr("height", d=> d[3]*segment_height)
        // .attr("fill", d=> d[4])
        .attr("fill", 'grey')
        .attr("rx", rx)
        .attr("ry", ry)
        .attr("opacity", opacity_mod)
        .attr("stroke", "black")
        .attr("stroke-width", stroke_width)

    var color_bar_width = 20;
    var color_bar_height = 180;

    // find difrence in magnintude and use as steps
    var max_exp = expo(max_intensity,1);
    var min_exp = expo(min_intensity,1);
    var max_exp_last_char = max_exp.substring(max_exp.length-1);
    var min_exp_last_char = min_exp.substring(min_exp.length-1);
    var steps = Math.abs(max_exp_last_char - min_exp_last_char)+1;

    var step = (Math.log(max_intensity) - Math.log(min_intensity)) / (steps - 1);
    log_steps = [];
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
    
    //rect.attr("transform", "translate(" + 0 + "," + (height/2 - margin.top + 100) + ")");
    });
}
