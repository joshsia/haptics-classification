# haptics-classification

## About

As humans, we use our sense of touch to help us classify and understand objects that we interact with. Touch is a useful alternative to vision when we can't see any object (due to poor lighting or occlusion) or an object property cannot be identified by sight (e.g. weight). This project investigates whether robots can also identify objects by touch using haptics data such as pressure, vibration, electrode impedances, etc in MATLAB. The objects used in this project are: acrylic, black foam, car sponge, flour sack, kitchen sponge and steel vase.

Hierarchical clustering was performed on the pressure, vibrations and temperature data using the standardized euclidean distance as the distance metric. The clusters found are shown below.

![hier-clustering](https://github.com/joshsia/haptics-classification/blob/main/results/clustering/pvt_cluster_3dplot_seuclidean.png)

An ensemble of trees (with bagging) was also trained to classify the objects based on electrode impedences. The confusion matrix of the model is shown below.

![conf-matrix](https://github.com/joshsia/haptics-classification/blob/main/results/clustering/elec_confusion_chart.png)

## Credits
The data used in this project was collected by the University of Pennsylvania's GRASP lab (PR2 robot which uses two *BioTac* tactile sensor fingertips). The dataset was cleaned by Ben Richardson, a PhD student at the Max Planck Institute for Intelligent Systems and was distributed by Dr Ad Spiers at Imperial College London as part of the Computer Vision and Pattern Recognition module.

This project was co-authored with Mr Yi Teng Voon.
