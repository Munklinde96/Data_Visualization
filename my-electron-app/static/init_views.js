const selectedProteinKey = "SELECTED_PROTEIN_KEY";
const selectedSamplesKey = "SELECTED_SAMPLES_KEY";
const selectedNormalizationKey = "SELECTED_NORMALIZATION_KEY";
const selectedMinModificationCountKey = "SELECTED_MIN_MODIFICATION_COUNT_KEY";
const tabSelectedKey = "TAB_SELECTED_KEY";

function getSelectedProtein(){
    var protein = window.localStorage.getItem(selectedProteinKey);
    console.log("protein: "+protein);
    return protein;
}

function setSelectedProtein(protein){
    window.localStorage.setItem(selectedProteinKey, protein);
    return;
}

function getSelectedSamples(){
    var selectedSamples = window.localStorage.getItem(selectedSamplesKey);
    if(selectedSamples && selectedSamples !== ""){
        var jsonSamples = JSON.parse(selectedSamples);
        console.log(jsonSamples);
        if(jsonSamples){
            jsonSamples.samples.sort();
            return jsonSamples.samples;
        }
    }
    return [];
}

function setSelectedSamples(samples){
    if(samples && samples.length !== 0){
        const Samples = {
            "samples": samples,
        }
        var jsonSamples = JSON.stringify(Samples);
        window.localStorage.setItem(selectedSamplesKey, jsonSamples);
        return;
    }
    window.localStorage.setItem(selectedSamplesKey, "");
}

function getSelectedNormalization(){
    var normalization = window.localStorage.getItem(selectedNormalizationKey);
    console.log("normalization: "+normalization);
    if(normalization){
        return normalization
    }
    return "";
}

function setSelectedNormalization(normalization){
    window.localStorage.setItem(selectedNormalizationKey, normalization);
    return;
}

function getMinModificationCount(){
    var modCount = Number(window.localStorage.getItem(selectedMinModificationCountKey));
    console.log("modCount: "+modCount);
    return modCount;
}

function setMinModificationCount(count){
    window.localStorage.setItem(selectedMinModificationCountKey, count.toString());
    return;
}

function getTabPressed(){
    var tabSelected = window.localStorage.getItem(tabSelectedKey);
    if (tabSelected === "true"){
        window.localStorage.setItem(tabSelectedKey, "false");
        return true;
    }
    return false;
}

function setTabPressed(){
    window.localStorage.setItem(tabSelectedKey, "true");
}

function setDocumentLabels(){
    console.log("setting document labels");
    // inputs
    try {
        document.getElementById("min_mod_count_input").value = getMinModificationCount();
        document.getElementById("normalization_select").value = getSelectedNormalization();
    } catch(e){}
    // state box
    try {
        document.getElementById("state_selected_protein").textContent = getSelectedProtein();
        document.getElementById("state_selected_normalization").textContent = getSelectedNormalization();
        document.getElementById("state_selected_min_count").textContent = getMinModificationCount();
        if(getSelectedSamples().length === 0){
            document.getElementById("state_selected_samples").textContent = "All samples";
        } else {
            document.getElementById("state_selected_samples").textContent = getSelectedSamples().toString();
        }
    } catch(e){}
}

(function() {
    
 })();