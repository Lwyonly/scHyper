# scHyper: Reconstructing cell-cell communication through hypergraph neural networks 

## What is scHyper?

![workflow](../../桌面/研/workflow.jpg)

We developed scHyper, which aims to decode context-driven intercellular communication from scRNA-seq data by constructing hypergraph networks and applying hypergraph neural network models . (i) scHyper integrates the Cellinker and CellChatDB to obtain rich L-R interaction relationships and signaling pathway annotations. (ii) scHyper constructs intercellular communication networks in the form of hypergraphs, analyzing intercellular communication from a global perspective. (iii) scHyper uses hypergraph neural network models to learn intercellular communication patterns and capture nonlinear relationships and complex interaction patterns. (iv) scHyper implements nonparametric tests to evaluate the statistical significance of interactions between cell types and thus identifies intercellular communication relationships that are biologically significant.

## How to use scHyper?

First, we construct the intercellular communication tensor and hypergraph network. The relevant code can be found in the dataprocess.py. Through data processing, we generated the following three files，and the user needs to specify the save path, such as save_path='to/path/data/demo'.

------test_data.npz;

------train_data.npz;

------use_to_predict.npz.

Second, we train model by graph neural networks. To run the code:

```
cd model
python main_torch.py --data demo 
```

Users can experiment by modifying the following parameters:

The **--data** argument can take a folder that contains the three files from the first step. 

The **--save-path** argument specifies the path where the model will be saved.

The **--load-path** argument specifies the loading paths for the training set, test set and prediction arrays (from Step 1).

Finally, we used nonparametric tests to identify significant intercellular communications. The relevant code can be found in the dataprocess.py. 

## Tutorial

To facilitate understanding and use of scHyper, we provide a simple demo, including count and meta, which provides tutorials on data processing, results and visualization.



