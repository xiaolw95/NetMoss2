# NetMoss2
NetMoss2 is the new version of NetMoss, which is a tool developed for integrating large-scale data and identifying disease associated biomarkers based on network algorithm.         

In this new version, both single file and 

For more information, please see paper "Large-scale microbiome data integration enables robust biomarker identification" published in Nature Computational Science.    


## Contents  
- [Installation](#installation)     
- [Basic Usage](#basic-usage)     
- [Input](#input)     
- [Output](#output)     
- [Classification](#classification)    
- [Example](#example)       

## Installation    
Installation with `devtools`     
```
library(devtools)
install_github("xiaolw95/NetMoss2")
library(NetMoss2)
```

## Basic Usage     
The NetMoss2 function is used to calculate NetMoss score of significant bacteria between case and control groups. Users are demanded to provide four directories or files as follows:      
```
NetMoss(case_dir = case_dir,    
        control_dir = control_dir,    
        net_case_dir = net_case_dir,   
        net_control_dir = net_control_dir)   
```
`case_dir:`  the directory or a single file of case data.     
`control_dir:`  the directory or a single file of control data.      
`net_case_dir:`  the directory or a single file of case network.      
`net_control_dir:`  the directory or a single data of control network.      


## Input     
Abundance or network matrix should be included in the input.    

##### Abundance Table
`case_dir` or `control_dir` includes abundance matrix which refers to the relative abundance of case or contol samples, with the row as bacteria and the column as samples. Abundance file can be processed from raw sequence using [QIIME2](https://qiime2.org/), [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn) or other tools.       
| taxon_names   | sample1 | sample2 | sample3 |    
|  ---  |  ---  |  ---  |  ---  |       
|   taxon1    |    60   |    20   |   10    |       
|   taxon2    |    30   |    77   |   89    |    
|   taxon3    |    0    |    23   |   15    |      
|   ... ...   |         |         |         |          

##### Network Matrix
`net_case_dir` or `net_control_dir` includes network matrix which refers to the adjacency matrix of correltaion between the bacteria. Microbial correlation can be deduced from any tools for which [SparCC](https://github.com/bio-developer/sparcc) or [SPIEC-EASI](https://github.com/zdk123/SpiecEasi) are especially recommended.     

|          | taxon1 | taxon2 | taxon3 |      
|  ------  | -----  | -----  | -----  |      
|  taxon1  |    1   |  -0.3  |  0.5   |      
|  taxon2  |  -0.3  |    1   |  0.67  |      
|  taxon3  |   0.5  |  0.67  |    1   |      
|  ... ... |        |        |        |     

##### Network construction
For convenience, we also provide a `netBuild` function to build microbial networks from abundance tables. To use this function, users are asked to provide abundance directories (contain case and control abundance tables). Network matrix will be output to the same directories automatically. For single file usage, users are asked to provided the abundance matrix only.          
```
netBuild(case_dir = case_dir,
      control_dir = control_dir,
      method = "sparcc")
```
`case_dir:`  the directory or a single file of case data.      
`control_dir:`  the directory or a single file of control data.          
`method:` the method to build networks. "sparcc" and "pearson" strategy are provided to choose.


## Output
The output of the NetMoss is a table of NetMoss score for each taxon:     
| taxon_names | NetMoss_score |      p.val    |      p.adj    |    
|  ---------  | ------------  |  ------------ |  -----------  |       
|    taxon1   |      0.98     |  5.703785e-09 |  2.335836e-08 |       
|    taxon2   |      0.7      |  1.467413e-04 |  2.629116e-04 |    
|    taxon3   |      0.32     |  2.018237e-04 |  3.542211e-04 |     
|    ... ...  |               |               |               |    

`taxon_names:` the name of the bacteria.         
`NetMoss_score:`  the NetMoss of the bacteria gets.      
`p.val:` the P value for the NetMoss score.   
`p.adj:` the adjust P value for the NetMoss score.    


## Classification       
In this section, we provide a pipeline to classify case and control groups based on the NetMoss markers. Iterative training and 10-fold cross validation stpes are implemented to guarantee the markers contain network and abundance informations. For this reason, it will take a long time to process the real datasets which contain large samples. Please be patient.
```
netROC(case_dir = case_dir,
      control_dir = control_dir,
      marker = marker,
      metadata = metadata,
      plot.roc = T,
      train.num = 20)
```
`case_dir:` the directory of case datasets.     
`control_dir:` the directory of control datasets.    
`marker:` a table of combined markers identified by NetMoss.     
`metadata:`  a table of clinical informations for all studies.     
`plot.roc:`  a logical parameter. If TRUE then the combined ROC of the result of classification will be plotted.     
`train.num:`  a numerical parameter which refers to trainning times of the model. By default, it is set to 20.        

First of all, efficient markers should be selected manually from the NetMoss result by users. Generally, we recommend a less strict threshold for the sparse network.
Also, a metadata file contains disease or health information for each sample needs to be inculded. The format should be like this:     
|  sample_id |   type  | study |     
|  ------  | -----  | -----  |     
|  SRRXXXXX  | disease | study1 |      
|  SRRXXXXX  | disease | study2 |       
|  SRRXXXXX  | healthy | study1 |        
|  ... ... |        |        |  

After preparing the two files, classification can be realized using the function `netROC`:     
```
marker = data.frame(result[which(result$NetMoss_Score > 0.3),])       
rownames(marker) = marker$taxon_names        
metadata = read.table("metadata.txt",header = T,sep = '\t',row.names = 1)     
myROC = netROC(case_dir,control_dir,marker,metadata)     
```

The result of the classfication is a table includes true positive rate and false positive rate:     
| threhold |  TPR  |  FPR  |      
|  ------  | ----- | ----- |      
|     0    |   1   |   1   |       
|    0.01  |  0.97 | 0.99  |       
|    0.03  |  0.9  | 0.87  |        
|  ... ... |       |       |  

A combined ROC will be ploted if the parameter `plot.roc` is set to be true.     

<img src="https://github.com/xiaolw95/NetMoss/blob/main/NetMoss_ROC.png" width = "500px">   


## Example
### example for multiple files
We have provided a small dataset to test the function.     
1. Download from the testthat directory () directly. 
Or get the dataset using `git clone` commond in `Linux`:      
```
git clone https://github.com/xiaolw95/NetMoss2.git     
cd NetMoss2/tests/testthat
```

2. After getting the dataset, the NetMoss score can be easily calculated using the `NetMoss2` function:       
```
##setwd('path-to-testthat-directory')
case_dir = paste0(getwd(),"/case_dir")
control_dir = paste0(getwd(),"/control_dir")
net_case_dir = paste0(getwd(),"/net_case_dir")
net_control_dir = paste0(getwd(),"/net_control_dir")
result = NetMoss(case_dir = case_dir,    
        control_dir = control_dir,    
        net_case_dir = net_case_dir,   
        net_control_dir = net_control_dir) 
```   

### example for single file
If users only have single file for case and control groups, NetMoss2 can also be used to identify significant biomarkes.   
```
data(testData)
nodes_result = NetMoss(case_dir = mydata[[1]],
                       control_dir = mydata[[2]],
                       net_case_dir = mydata[[3]],
                       net_control_dir = mydata[[4]])
my_result = nodes_result[[1]]
```
