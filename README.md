########### Scripts for CAGI 5 PCM1 challenge assessment ##################
 

## REQUIREMENTS
```

R:
>install.packages('ROCR')
>install.packages('plotrix')
>install.packages('ggthemes')
>install.packages('ggplot2')
```
## USAGE

to run all analyses, you have just to run the following command: Rscript ./src/CAGI_assessment_main.R
      This script will check submission format, computes all statistics and makes plots
      This script expects a folder structure like this to run:
      ./src <- contains all scripts
      ./results <- will contain all performance tables and plots
      ./data <- has to contain 3 folders: experimental_value, submissions, template 
      
      These are Input needed:
            # experimental values file  in ./data/experimental_value
            # submission template in ./data/template
            # submission files in ./data/submissions

      Running the script these Output files will be generated in ./results
           # confusion matrix of each submission
           # heatmap of Correlation among submissions (MO or MO+WT p value)
           # heatmap containing all performance indices (BACC, MCC, F1, TPR, TNR)
           # plot with percentage of correct predictions for each variant
           # scatterplot with all submissions p value and SD (MO or MO+WT p value)
           # plots of the distribution of BACC/F1/MCC values of random predictions
           # table containing all performance indices (BACC, MCC, F1, TPR, TNR)
           # table containing prediction rank for all performance indices (BACC, MCC, F1, TPR, TNR)
           # table with performance index and rank
           # table with one confusion matrix for each submission
           
