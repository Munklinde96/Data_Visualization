const selectedProteinKey = "SELECTED_PROTEIN_KEY";
const selectedSamplesKey = "SELECTED_SAMPLES_KEY";
const selectedNormalizationKey = "SELECTED_NORMALIZATION_KEY";
const selectedMinModificationCountKey = "SELECTED_MIN_MODIFICATION_COUNT_KEY";
const sortingOrderProteinsKey = "SORTING_ORDER_PROTEINS";
const tabSelectedKey = "TAB_SELECTED_KEY";

function getSortingOrderProteins(){
    var order = window.localStorage.getItem(sortingOrderProteinsKey);
    if(order && order !== ""){
        var jsonOrder = JSON.parse(order);
        if(jsonOrder){
            return jsonOrder.order;
        }
    }
    return [];
}

function setSortingOrderProteins(order){
    if(order.length != 0){
        const Order = {
            "order": order,
        }
        var jsonOrder = JSON.stringify(Order);
        window.localStorage.setItem(sortingOrderProteinsKey, jsonOrder);
        return;
    }
    window.localStorage.setItem(sortingOrderProteinsKey, "");
}

function getSelectedProtein(){
    var protein = window.localStorage.getItem(selectedProteinKey);
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