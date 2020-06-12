# 10 reasons to study a gene

This is the supplementary code for the paper "10 reasons to study a gene" by
Ionut Sebastian Mihai, Debojyoti Das, Gabija Maršalka and Johan Henriksson (corresponding).

*Laboratory of Molecular Infection Medicine Sweden (MIMS), 
Umeå Centre for Microbial Research, Department of Molecular Biology, Umeå University, Umeå, Sweden*

### Reproducing

The features are produced in feature_summarize.R. The smaller input files are already in the input/-folder. Bigger files need to be downloaded. There are files in
this directory describing this process. The generated features are otherwise provided with the paper supplementary.

The linear model is provided in model_linear.R - it also has some rendering for the nonlinear models.

The nonlinear model fitting is done by the Jupyter notebooks: model_neuralnetwork.ipynb and model_xgboost.ipynb

## Deployment of server

In the website directory, simply run "make" (assumes you have copied required files into its data folder; this likely needs our guidance).
All the dependencies are in requirements.txt

## Built With

* The majority of the analysis is in R using many fairly standard packages. SQLdf is a great package!
* Machine learning was done with Python using PyTorch and XGBoost
* There is some Greta code here for piece-wise linear model fitting not mentioned in the paper
* Webserver uses Dash
* Figures were made mainly in Inkscape
* Paper written using google docs and the paperpile reference manager

## License

This project is licensed under the MIT License.
