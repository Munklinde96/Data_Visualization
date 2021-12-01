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
    process.env['ELECTRON_DISABLE_SECURITY_WARNINGS'] = 'true';
    initializeBackend().onreadystatechange = (e) => {
            const mainWindow = new BrowserWindow({
                width : 1200,
                height: 800,
                webPreferences: {
                    nodeIntegration: true,
                    webviewTag: true,                
                    contextIsolation: false
                }
            });
            mainWindow.maximize();
            mainWindow.loadFile(path.join(__dirname, "segment_plot.html"));
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