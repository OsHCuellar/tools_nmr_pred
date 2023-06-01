# tools-example

Web interface of the [nmr_pred](https://github.com/maltefranke/nmr_pred) code. The current version only allow the user to visualize and analyze the chemical diversity of the datased used to train the model.

# Setup
Get a local copy of the code
```
gh repo clone OsHCuellar/tools_nmr_pred
```
inside the main directory build the docker image
```
cd tools_nmr_pred
sudo docker-compose build
```
and run the container
```
sudo docker-compose up
```
You can acces to the web-tool by opening the direction that will be given in a browser.
