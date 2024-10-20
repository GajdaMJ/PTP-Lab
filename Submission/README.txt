Hello. This is the submission document for Adiabatic reactors 1.

CSTR_Code has all the code we used for the CSTR model. The file "CSTR_Model.py" has the simplified model. The file "CSTR_Model_temp_step_change.py" has the code for the step change. All the other files were used for creating plots for the report. The code should run but it might not work if the paths change. Make sure that the paths are correct. 




The PBR_Code folder contains all of the code relavant for the PBR reactor. In this folder, you can find: 
PBR_model_only_final_plot: This file runs the PBR code and plots only the final plots as a 2 by 4 figure. 
PBR_model_step_change_only_final_plot: This file runs the PBR code for the step change temperature, only plotting the final 2 by 4 figure. 
PBR_model_step_change: This file runs the PBR code for the step changes in temperature plotting everything. 
PBR_model: This file runs the basic PBR code plotting everything. 
PBR_number_of_tanks: This file is uses the sum of least squares to determine the ideal number of tanks in series to use. 

Most important code for PBR models are: PBR_model and PBR_model_step_change

All the data obtained is stored in the data folder.
