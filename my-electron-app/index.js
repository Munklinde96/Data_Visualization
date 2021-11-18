const { app, BrowserWindow, net } = require("electron");
const path = require("path");
const XMLHttpRequest = require('xhr2');

const backendAddr = "http://127.0.0.1:5000"

function initializeBackend(){
    const Http = new XMLHttpRequest();
    const url=backendAddr+'/';
    Http.open("GET", url);
    Http.send();
    return Http;
}

const loadMainWindow = () => {
    initializeBackend().onreadystatechange = (e) => {
            const mainWindow = new BrowserWindow({
                width : 1200,
                height: 800,
                webPreferences: {
                    nodeIntegration: true
                }
            });
            mainWindow.maximize();
            mainWindow.loadFile(path.join(__dirname, "index.html")).then((_) => {
                //console.log("site is ready");
            }).then((_) => {
                
            });   
    }
}

app.on("ready", loadMainWindow);

app.on("window-all-closed", () => {
    if (process.platform !== "darwin") {
      app.quit();
    }
});

app.on("activate", () => {
    if (BrowserWindow.getAllWindows().length === 0) {
        loadMainWindow();
    }
});  