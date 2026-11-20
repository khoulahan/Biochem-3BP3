## Lab # 9 Machine Learning

## Table of Contents
1. [The Problem](#problem)
2. [Defining Training and Testing Sets](#define)
3. [Training a Classification Model](#training)
4. [Evaluating a Classification Model](#testing)
5. [Tuning Classification Models](#tuning)

<a name="intro"></a>
## Introduction

The goal of this lab is to build a machine learning classifier to predict breast cancer subtypes based on RNA-Seq data using the public Galaxy platform.

**Lectures** - [Lecture 10](https://github.com/agmcarthur/Biochem-3BP3/blob/master/Lectures/Houlahan_MachineLearning.pptx) Introduction to Machine Learning ([~42 minute video](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/mcarthua_mcmaster_ca/Ed6icPcJcIJKu21HWdkOtZYBlaH28uz29GnlFIFAZvEu9w))

**Flash Updates**
* *Random Forest* 
* *Logistic Regression* 
* *Natural Language Processing* 

**Demo Videos**
* [Logging onto the Cluster](https://mcmasteru365-my.sharepoint.com/:v:/g/personal/houlahke_mcmaster_ca/IQC2ftW_ZyPVS5qzMttT4FYvAX0pB_x6evXgf2O9YmDjYWI?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=WCQxQd) ~2 minutes


**Background Reading** (optional)
* Whalen et al. 2022. Navigating the pitfalls of machine learning in genomics. [Nat Rev Genet. 23:169-81](https://pubmed.ncbi.nlm.nih.gov/34837041/)
* Thennavan et al. 2021. Molecular analysis of TCGA breast cancer histologic types. [Cell Genomics. 1:100067](https://pubmed.ncbi.nlm.nih.gov/35465400/) 

**Computer Resources**
* The lab can be completed on the *cluster*. Refer back to Lab 1 and 2 to refresh using command line arguments in Linux and R.

**Grading**
* Questions are for your learning and are not graded
* Problems are worth 5 points each (-1 for each error)
* Submit your answers to the Problems, plus any supplmental multiple choice questions, on **A2L Quizzes** before the deadline
* An answer key to Questions and Problems will be provided on A2L after the deadline

<a name="problem"></a>
## The Problem

Breast cancer is a molecularly heterogeneous disease. Treatment decisions are dependent on the presence of three protein receptors: [HER2, ER and PR](https://www.cancer.org/cancer/types/breast-cancer/understanding-a-breast-cancer-diagnosis/breast-cancer-hormone-receptor-status.html). The presence or absence of these three receptors categorizes tumours into three subtypes. 

| Subtype | Receptor | Prevalence | Survival (5 year) | Treatment | 
|----|:----:|:----:|:----:|----|
| ER+ | ER+ HER2- | 80% |>90% | Hormone therapy | 
| HER2+ | HER2+ |15-25% |~85% | Anti-HER2 therapy | 
| Triple Negative (TNBC) | ER- PR- HER2- |20% |77% | Chemotherapy 

In lab 2, you attempted to assign subtypes to each breast tumour by stratifying tumours based on the mRNA abundance of just *ESR1* and *ERBB2*. Now, we are going to get a little more sophisticated. We are going to build a machine learning model that can differentiate ER+ versus TNBC tumours from transcriptomic profiling data (i.e., RNA-sequencing) to ensure patients receive the most effective treatment. 

### Data

As in lab 2, we are going to be use transcriptomic profiling of 48 genes in 848 tumours which have been clinically defined as ER+ or TNBC. The data can be found in `/workspace/lab/studentlab/lab2_R/breast_cancer_rna.txt`.

The data is a tab-delimited file where each row corresponds to a tumour and the columns correspond to either the sample ID (`CLID`), the subtype of the tumour (`subtype`) or the mRNA abundance of a gene (e.g., `EGFR`). Here the subtype is encoded as a binary variable where `1 = TNBC` and `0 = ER+`. As before, let's start by reading in the data into R. Remember to first log onto an interactive node before launching R. 


```Bash
# request an interactive node on the cluster with 2G that will run for 24 hours
salloc --mem=2G --time=24:00:00 --partition=classroom
srun --pty bash
```

Launch R and read in data.

```R
library(tidyverse)
library(caret)

# read in file
rna <- read_tsv("/workspace/lab/studentlab/lab2_R/breast_cancer_rna.txt")
```

**Problem #1. Is the dataset balanced (ER+ versus TNBC)? What do we need to consider when training a model on this dataset? (hint: look up the `table()` command in R)**

<a name="define"></a>
## Defining Training and Testing Sets

First, we want to set our training and testing cohorts. Here we are going to use 80% of our cohort for training while holding out 20% for testing. 

**Question #1. Here we are leaving out 20% of the data as a test set. What is another validation method could we use to evaluate our model?**

```R
# set your own random seed here by replacing the x's with a random set of numbers 
set.seed(xxx)

# make sure treats subtype as binary, set as factor
rna$subtype <- as.factor(rna$subtype)

# split data into train and test 
train_index <- createDataPartition(
    rna$subtype,
    p = 0.80
    )
train_data <- rna[train_index[[1]],]
test_data <- rna[-train_index[[1]],]
```

**Note:** Set a random seed. Some steps in machine learning are random and leverage random number generators. For example, randomly splitting cohort into training and testing sets. To ensure someone else can reproduce your analyses, it is important that this randomness is reproducible. You can set a random seed to ensure random number generators are initialized the same way each time. This will ensure even the randomness is reproducible.

**Question #2. What is the ratio of ER versus TNBC tumours in the training dataset? What is the ratio in the testing dataset? How would this differ if we used the non-stratified shuffle strategy?**

> Flash Update - Logistic regression

<a name="training"></a>
## Training a Classification Model

We are going to build a classifier using a logistic regression model. A logistic regression model is a type of generalized linear regression that expects the response variable to be binary (0/1), i.e. ER+ versus TNBC tumour.

Train a logistic regression model on the training dataset using the **Generalized Linear Models** tool. Keep in mind **not** to include the sample IDs (e.g. TCGA-XX-XXXX-XX) as a feature in the dataset as we want to train the model to predict subtype based on the mRNA abundance of the 48 genes only.

We are using the `caret` package in R (which was designed to mirror `scikit-learn ` in python). The caret package allows us to set up our train/test data as well as our cross validation strategy and then we can test multiple machine learning models with the same training infrastructure. Here, we are going to set up a simple 5-fold cross validation strategy.

```R
set.seed(xxx)
# set up cross validation strategy
cv_strategy <- trainControl(
    method = "cv",
    number = 5
    )
# fit logistic regression
log_fit1 <- train(
    subtype ~ ., 
    data = train_data[,-1], # removing the first column
    method = "LogitBoost",
    trControl = cv_strategy
    )
```

**Question #3. What would happen if we accidentally left the sample ids as a feature in the model training?**

<a name="testing"></a>
## Evaluating a Classification Model

Now that we have trained our model on our training set, we want to test it on the held-out 20% test. Remember to remove the sample ID and subtype columns when predicting on the test set.

```R
# predict on held out data 
log_preds <- predict(
    log_fit1, 
    newdata = test_data[,-c(1:2)],
    type = 'raw'
    )
# calculate confusion matrix
confusionMatrix(
    data = log_preds,
    reference = test_data$subtype,
    positive = '1',
    mode = "prec_recall"
)
```

**Question #4. What is the accuracy of your model?**

**Problem #2. Is accuracy the best metric to use to evaluate your model? What other metrics would you calculate? How do the additional metrics compare to accuracy? If they differ, why?**

**Question #5. How do your machine learning predictions compare to your subtype assignments in Lab 2?**

> Flash Update - Random Forest

<a name="tuning"></a>
## Tuning Classification Models

Of course, logistic regression is not the only classification model available. Try training and testing one other classification algorithm using the same steps outlined above but changing the `method` in the training step. Check out the `caret` [documentation](https://topepo.github.io/caret/train-models-by-tag.html) to find all the possible models caret supports.

**Question #6. Which classification algorithm did you choose and why? How did its performance compare to logistic regression?**

<a name="prediction"></a>
## Predicting on Microarray Data

Comprehensively evaluating model performance requires testing the model on a completely unseen cohort. This evaluation more accurately reflects the anticipated performance when the model is applied in a real-world setting. **You should never look at the performance of the model on this cohort before picking a final model.**

In our case, we are going to apply our model, trained on RNA-sequencing, to predict the breast cancer subtype of an independent cohort of 1,733 tumours profiled with microarrays. The mRNA abundance of the same 48 genes have been measured. The data can be found: `/workspace/lab/studentlab/lab9_machine_learning/breast_cancer_microarray_prediction.txt`.

The data is a tab-delimited file where each row corresponds to a tumour and the columns correspond to either the sample ID (`CLID`), the subtype of the tumour (`subtype`) or the mRNA abundance of a gene (e.g. `EGFR`). Again, here the subtype is encoded as a binary variable where `1 = TNBC` and `0 = ER+`.

Using the same approach as outlined in `Evaluating a Classification Model`, evaluate your final prediction model (either the logistic regression model or the model you trained in `Tuning Classification Models`) on the independent microarray data.

**Problem #3. How well did your model perform on the microarray data? What are some factors that may influence model performance?**