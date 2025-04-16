These R scripts were written by Ron Geller to facilitate the analysis of antiviral or neutralizing sera concentration curves from the data obtained using the Incucyte software or a standard plate reader format. The scripts collapse data derived from neutralization, viability, or antiviral testing data into long format, followed by fitting of log.logistic models to find the effective or inhibitor concentration specified by the user (e.g. IC50, IC90). 

Ensure you have R (tested with 4.3.1), R-Studio (tested with 2024.12.1) and the needed R packages installed:

- drc (tested with 3.0-1)
- tidyverse (tested with 2.0.0)
- readxl (already part of tidyverse)
- gridExtra (tested with 2.3)
- scales (tested with 1.3.0)

To run the scripts:
1. Download the script to a new directory
2. Update the experimental_template excel file with your experimental setup. Do not change the position of the columns or rows. 
3. Update the conditions_file with the conditions of your experiment. The possible options for each parameter are described on the second page (information page). You must update the name of your data files and locations of the folders.
4. Start the R project. All paths are relative to this folder. 
5. Open the execute_ec file. Run the lines for each step to ensure your conditions are correct and change the needed parameters.

- Step 1 will collapse your data to long format, joining it with the information from the experimental_template, writing it to a folder called collapsed_data in the directory you indicate in the out.dir column of the conditions_file (default is ../results/). 

- Step 2  will join all of your data into a single csv file in the case you have more than one plate scanned. It will write the collapsed combined data into the directory you specify. 

- Step 3 will read the collapsed data (or the combined collapsed data if you have more than 1 plate), fit the desired model (usually LL.2() or LL.3() for a two parameter or three parameter log logistic function since we are using data that is standardized to untreated controls) and generate a graph if you want. You have to modify the data.file to match where your collapsed data is (or the combined collapsed data). You can also control the number of graphs generated per page (4x4 is the max to see well; set with the n.row and n.col parameters).
The output you get is a directory called ec_data with a file for each effective concentration level and selected model (e.g. res_EC50_LL.3.csv) and another directory called ec_graphs with the graph of the fit for each effective concentration level and model if you select the plot = T parameter in the get_ec function.

- Step 4 will combine all the EC tables into one. It is optional but useful if you have more than 1 effective concertation level or model. 


Feel free to modify and distribute the scripts as desired (e.g. different plate formats, etc.)