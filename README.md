# Fast_GNSS_ReceiverMATLAB
A fast and highly sensitive GNSS receiver with a significant capability of detecting very weak GNSS signals (reflected signals) [MATLAB version]
Here are the instructions for a successful run:
•	Create two folders; ‘mat’ and ‘data’
•	In the ‘data’ folder create two new folders; ‘data_in’ and ‘data_out’ (or move the 'data_in' folder to the 'data' folder if already created)
•	Move all downloaded files/folders to the ‘mat’ folder except these three:
o	caCode.mat
o	data_file_list.txt
o	NTLab .bin file or other datasets
•	Move the three exception files mentioned above to the ‘data_in’ folder
•	Check and modify parameter settings in config_sdr_params.m file in the ‘cfg’ folder
•	Check the data_file_list.txt and add the input data files
•	Run the init.m file
