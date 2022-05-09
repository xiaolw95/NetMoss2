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
The NetMoss function is used to calculate NetMoss score of significant bacteria between case and control groups. Users are demanded to provide four directories or files as follows:      
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

##### Network Construction
For convenience, we also provide a `netBuild` function to build microbial networks from abundance tables. To use this function, users are asked to provide abundance directories (contain case and control abundance tables). Network matrix will be output to the same parent directories automatically. For single file usage, users are asked to provided the abundance matrix only.      

`NOTE:` The `netBuild` function will creat "net_case_dir" and "net_control_dir" directories and output the network results into them. If the same directories exist, files will be overwritten.     

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

#### Visualization
In this part, we provide a function to illustrate the results. The `netPlot` function will output two kind of images. The one is the visualization of the important NetMoss score, the other is a paired networks to demostrated the difference of structure between case and control networks.         
```
netPlot(result = nodes_result,
        tax = "g__Bacteroides")
```
`result:` the result from NetMoss function. This is a list contained NetMoss score and integrated networks.    
`tax:` the target taxon which the users are interested. The taxon name must be included in the input file.    

There are two types of visualization of NetMoss score: barplot and point plot.    

<img src="https://github.com/xiaolw95/NetMoss2/blob/main/pic/NetMoss_score.jpg" width = "800px">       

Also, there are two types of visualization to demostrated the the difference of structure between case and control networks.     

<img src="https://github.com/xiaolw95/NetMoss2/blob/main/pic/network1.jpg" width = "400px"> <img src="https://github.com/xiaolw95/NetMoss2/blob/main/pic/network2.jpg" width = "400px">      

## Classification       
In this section, we provide a pipeline to classify case and control groups based on the NetMoss markers. Iterative training and 10-fold cross validation stpes are implemented to guarantee the markers contain network and abundance informations. For this reason, it will take a long time to process the real datasets which contain large samples. Please be patient.
```
netROC(case_dir = case_dir,
      control_dir = control_dir,
      marker = marker,
      metadata = metadata,
      plot.roc = TRUE,
      train.num = 20)
```
`case_dir:` the directory or a single file of case data.     
`control_dir:` the directory or a single file of control data.    
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


#### Result
After preparing the two files, classification can be realized using the function `netROC`     

The result of the classfication is a table includes true positive rate and false positive rate:     
| threhold |  TPR  |  FPR  |      
|  ------  | ----- | ----- |      
|     0    |   1   |   1   |       
|    0.01  |  0.97 | 0.99  |       
|    0.03  |  0.9  | 0.87  |        
|  ... ... |       |       |  

A combined ROC will be ploted if the parameter `plot.roc` is set to be true.     

<img src="https://github.com/xiaolw95/NetMoss2/blob/main/pic/NetMoss_ROC.png" width = "800px">   


## Example
### example for multiple files
We have provided a small dataset to test the R package. In our testthat directory, both abundance files and network files, as well as the metadata, are included.       
1. Download from the testthat directory () directly. 
Or get the dataset using `git clone` commond in `Linux`:      
```
git clone https://github.com/xiaolw95/NetMoss2.git     
cd NetMoss2/tests/testthat
```

2. After getting the dataset, the NetMoss score can be easily calculated using the `NetMoss` function:       
```
#setwd('path-to-testthat-directory')   ####set the directory to testthat

#read directory
case_dir = paste0(getwd(),"/case_dir")
control_dir = paste0(getwd(),"/control_dir")
net_case_dir = paste0(getwd(),"/net_case_dir")
net_control_dir = paste0(getwd(),"/net_control_dir")

#construct networks  ####if files exist, skip
#netBuild(case_dir = case_dir,
#         control_dir = control_dir,
#         method = "sparcc")

#calculate NetMoss score
nodes_result = NetMoss(case_dir = case_dir,    
        control_dir = control_dir,    
        net_case_dir = net_case_dir,   
        net_control_dir = net_control_dir) 
result = nodes_result[[1]]     ####NetMoss score result

#plot networks
netPlot(nodes_result,tax = "g__Enterobacter")    ####image saved

#plot roc 
marker = data.frame(result[which(result$p.adj < 0.05),])
marker = data.frame(marker[which(marker$NetMoss_Score > 0.3),])   ####marker selection
rownames(marker) = marker$taxon_names
metadata = read.table("metadata.txt",header = T,sep = '\t',row.names = 1)
myROC = netROC(case_dir = case_dir,
               control_dir = control_dir,
               marker = marker,
               metadata = metadata,
               plot.roc = TRUE, 
               train.num = 20)    ####image saved
```   

### example for single file
If users only have single file for case and control groups, NetMoss2 can also be used to identify significant biomarkes.   
```
#setwd("your-directory")

#load dataset
data(testData)

#contruct networks    ####if files exist, skip
#netBuild(case_dir = mydata[[1]],
#         control_dir = mydata[[2]],
#         method = "sparcc")     

#calculate NetMoss score
nodes_result = NetMoss(case_dir = mydata[[1]],
                       control_dir = mydata[[2]],
                       net_case_dir = mydata[[3]],
                       net_control_dir = mydata[[4]])
result = nodes_result[[1]]   ####NetMoss score result

#plot networks
netPlot(nodes_result,tax = "g__Bacteroides")    ####image saved

#plot roc 
marker = data.frame(result[which(result$p.adj < 0.05),])
marker = data.frame(marker[which(marker$NetMoss_Score > 0.3),])   ####marker selection
rownames(marker) = marker$taxon_names
metadata = mydata[[5]]
myROC = netROC(case_dir =  mydata[[1]],
               control_dir =  mydata[[2]],
               marker = marker,
               metadata = metadata,
               plot.roc = TRUE, 
               train.num = 20)    ####image saved
```
