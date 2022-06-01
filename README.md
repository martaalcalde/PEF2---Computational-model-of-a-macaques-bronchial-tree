# PEF2---Computational-model-of-a-macaques-bronchial-tree
Codes for the final project of the subject Projects of Engineering Physics 2 of the degree Engineering Physics (UPC). The model has been built using MATLAB.
By Marta Alcalde and Núria Mercadé.

BT_Model1.m: Computational macaque bronchial tree reaching to the 20 pulmonary segments using Model 1. It generates a matrix BT_model1.mat, containing the information of each one of the branches, and Figure 4 of the Article.

BT_Model2.m: Computational macaque bronchial tree reaching to the 20 pulmonary segments using Model 2. It generates a matrix BT_model2-mat, containing the information of each one of the branches, and Figure 5 of the Article.

Bronchis.m: Fills the lungs iteratively using Vegué and Català method. Needs BT_model1.mat or BT_model2.mat. Generates a matrix called bronchis_model1.mat/bronchis_model2.mat depending on the model previously used.

FiguraBT.m: Generates figures 6 and 7 of the Article using bronchis_model1.mat or bronchis_model2.mat.

setments.mat: 3D matrix containing the information of CT images. It is filled with integers from 0 to 20. 0 indicates that the point is outside the lung. Integers from 1 to 20 indicate the segment where the point is found.
