'''
This script serves as a second form of testing code which verifies 
that the quantity of samples in the Xena beta release hub's datasets 
are accurate. The GDC will often update existing projects by adding 
new samples or redacting innacurate ones. In the future, this program 
can be used to determine which datasets need updates to their data. 

The program returns a 'detect_missing_datasets.tsv' file that displays
which datasets from which project need updating. The file contains 
three main categories: 'missing_datasets' for datasets with samples in 
the GDC, but do not exist in the Xena hub. 'datasets_with_wrong_sampleN' 
for datasets that have an incorrect number of samples. This may be due 
to GDC adding or redacting new samples, or that the dataset doesn't 
exist on the Xena hub. 'missing_datatypes' displays the 'data_type', 
'workflow_type','platform', and 'experimental_strategy' joined together 
with a '/' of datasets that are missing from the Xena hub or are 
completely  missing from '_XENA_GDC_DTYPE' which acts as a reference 
to create datasets. 

The program also creates a logging file called 'missing_submitter_ids.json'. 
This file contains every submitter ID that is missing from the Xena hub, as
well as ones removed from the GDC. They are organized by project ID and the 
dataset they belong to. 

Usage Instructions:

To run the script on all projects in the GDC, run the file on your command line without any arguments

Example:
	python3 DetectMissingDatasets.py

There is a test mode for if you would like to only run the program on 
one project. To use the test mode, just add the project ID of the 
project as a second argument in the command line.

To run the script on an individual project:
[This file name] [project ID]

Example:
	python3 DetectMissingDatasets.py TCGA-BRCA

'''


import xenaPython as xena
import pandas as pd
import requests
import json
import xena_dataset 
import logging
import time
from datetime import datetime
import sys

###########################Constants###################
#_XENA_GDC_DTYPE Used as a reference for which files belong to which dataset 
_XENA_GDC_DTYPE = xena_dataset.GDCOmicset._XENA_GDC_DTYPE
PROJECT_ENDPT = "https://api.gdc.cancer.gov/projects"
FILE_ENDPT = "https://api.gdc.cancer.gov/files"
UNUSED_DATA_TYPE = ["Slide Image", "Biospecimen Supplement", "Clinical Supplement", "Masked Intensities", "Pathology Report", "Isoform Expression Quantification", "Tissue Microarray Image"]
UNUSED_EXPERIMENTAL_STRATEGY = ["scRNA-Seq"]
TUMOR_DATA_TYPE = ["Copy Number Segment", "Masked Copy Number Segment", "Gene Level Copy Number", "Masked Somatic Mutation", "Allele-specific Copy Number Segment"]
HUB = "https://gdcbetarelease.xenahubs.net"
FILE_FIELDS = ["file_id","data_type","analysis.workflow_type","platform", "experimental_strategy", "cases.samples.submitter_id", "cases.samples.tissue_type", "cases.submitter_id"]
PROJECT_FIELDS = ["project_id", "released"]

#######################Logging#############################
start_time = time.time()
logging.basicConfig(filename = 'missing_submitter_ids.json', level=logging.INFO, format='%(message)s')
timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
logging.info(f'{timestamp}')
#######################Global###############
#Global variables initialized here
'''
Attributes:
	test (string): This variable contains the project ID when 'test_mode' is True. 
	test_mode (boolean): True when the script will only run on one project, specified by the user.
	project_id_list (list): Where all project IDs of projects that will be run on are contained.
	missing_submitter_ids(dict): submitter_ids that are missing from the Xena hub are saved here. 
	Organized by project ID and dataset. Used to form logging file.

	Below variables used to form TSV file. Each element represents the data for one project
	total_MDS(list): Datasets missing from the Xena hub. 
	total_CMDS (list): # of datasets missing from the Xena hub.
	total_DSMS (list): Datasets that are missing samples.
	total_CDSMS (list): # of datasets that are missing samples.
	total_MDT (list): Datatypes that are missing from the Xena hub.
	total_CMDT (list): # of datatypes that are missing from the Xena hub. 

'''
test = ''
test_mode = False
if len(sys.argv) > 1:
	test_mode = True
	test = sys.argv[1]


project_id_list = []
total_MDS = []
total_CMDS = []
total_DSMS = []
total_CDSMS = []
total_MDT = []
total_CMDT = []
missing_submitter_ids = {}

def project_request(fields):
	'''
	Makes an API call to the GDC for all project IDs in JSON format. Then converts it to a JSON-formmated dictionary.
	
	This API call is done using a post request and returns in JSON format.
	args:

	returns: 
		Dictionary: JSON-formatted dictionary that contains all project IDs from the GDC.
	'''

	filters = {}
	
	fields = ",".join(fields)
	params = {"filters": json.dumps(filters), "fields": fields,
	"format": "json",
	"size": "100"}
	response = requests.post(PROJECT_ENDPT,headers = {"Content-Type": "application/json"}, json = params)
	responseJson = json.dumps(response.json(),indent=2)

	return responseJson


def format_to_list(responseJson):
	'''
	This function extracts out the project IDs within the JSON-formatted dictionary returned by the GDC API.

	Since the project IDs are nested in multiple dictionaries and lists, they first need to be "unpeeled". 
	Only then can they be accessed and appended to a new list. 

	args:
		responseJson(dictionary): A dictionary that is created by the return of the GDC API call. 

	returns:
		list: List of all project IDs in the GDC.

	'''
	responseJson = unpeel(responseJson)
	project_id_list = []
	for item in responseJson:
		first_key = next(iter(item))
		first_value = item[first_key]
		project_id_list.append(first_value)
	return project_id_list

def unpeel(resJson):
	'''
	Removes unnecessary dictionaries from the JSON-formatted dictionary returned by the GDC API.

	All metadata is located in the dictionary corresponding to "data", and the list corresponding 
	to "hits" within that dictionary.

	args:
		responseJson(dictionary): A dictionary that is created by the return of a GDC API call. 
	returns:
		list: A list that contains dictionaries corresponding to each file or project's metadata. (NOT yet a list of metadata)
	'''
	resJson = json.loads(resJson)
	resJson = resJson.get("data")
	resJson = resJson.get("hits")
	return resJson


def file_request(project_id, fields):
	'''
	Calls the GDC API for all open access files in a project with the relevant file metadata.
	
	This API call is done using a post request and returns in JSON format.

	args:
		project_id(string): The project_id of the project that you are requesting file metadata for.
		fields(list): A list of types of metadata that will be requested from the GDC; currently this 
		function uses the constant 'FILE_FIELDS'.
	returns:
		list: A list that contains dictionaries corresponding to each file's metadata. (NOT yet a list of pure metadata)

	'''
	fields = ",".join(fields)
	filters = {
	"op":"and",
	"content":[
	{
		"op": "in",
		"content":{
		"field": "cases.project.project_id",
		"value": [project_id]
	}
		},
	{
		"op": "in",
		"content":{
		"field": "access",
		"value": ["open"]
		}
	}
	]
	
	}
	params = {"filters": json.dumps(filters), "fields": fields,
	"format": "json",
	"size": "100000"}
	response = requests.post(FILE_ENDPT,headers = {"Content-Type": "application/json"}, json = params)
	responseJson = json.dumps(response.json(),indent=2)
	responseJson = unpeel(responseJson)
	return responseJson



def create_data_set(project_id):
	'''
	This function creates a list of Data_set objects which each represent a potential data set 
	in the project. 

	Each dataset has a set of requirements for a sample to belong in that dataset. These 
	requirements are contained in '_XENA_GDC_DTYPE'. These requirements then become attributes
	of each Data_set object

	args:
		project_id(string): The project_id of the project that you are creating datasets for.
	returns:
		list of Data_set objects: A list of potential datasets that belong to a project.
	'''
	data_set_list = []
	for key, value in _XENA_GDC_DTYPE.items():
		data_type = ''
		workflow_type = ''
		platform = ''
		experimental_strategy = ''
		if "data_type" in value:
			data_type = value["data_type"]
		if "analysis.workflow_type" in value:
			workflow_type = value["analysis.workflow_type"]
		if 'platform' in value:
			platform = value["platform"]
		if 'experimental_strategy' in value:
			experimental_strategy = value["experimental_strategy"]
		data_set = Data_set(key, data_type,workflow_type, platform, experimental_strategy, project_id)
		data_set_list.append(data_set)
	return data_set_list

def file_metadata(file_list, data_set_list, project_id, missing_file_list):
	'''
	This function takes in the list of file metadata returned by the GDC API and creates a File object 
	for each file. Then each file is organized to a dataset, placed into a missing datatype list, or 
	thrown out depending on its metadata. 

	args:
		file_list(list): list where each element is a dictionary of one file's metadata. Returned from the 
		GDC API in JSON format.
		data_set_list(list of Data_set objects): A list of potential data sets in a project.
		project_id(string): The project_id of the project you are organizing files for. Only used to
		check if the project is "CPTAC-3", a special case. 
		missing_file_list(list of File objects): list of files which do not belong in any of the potential datasets


	'''
	for file in file_list:

		file_id = file["file_id"]

		submitter_id = ''
		data_type = ''
		workflow_type = ''
		platform = ''
		experimental_strategy = ''
		if "data_type" in file:
			data_type = file["data_type"]
		if "analysis" in file:
			workflow_type = file["analysis"]["workflow_type"]
		if 'platform' in file:
			platform = file["platform"]
		if 'experimental_strategy' in file:
			experimental_strategy = file["experimental_strategy"]
		if project_id == "CPTAC-3":
			if not (data_type in UNUSED_DATA_TYPE) and not ((data_type == "Copy Number Segment") and (workflow_type == "DNAcopy" )) and not(experimental_strategy in UNUSED_EXPERIMENTAL_STRATEGY):
				sample_info = file["cases"][0]["samples"]
				
				if data_type in TUMOR_DATA_TYPE:
					for i in range(len(sample_info)):
						if sample_info[i]["tissue_type"] == "Tumor":	
							submitter_id = file["cases"][0]["submitter_id"]
							break
				else:
					submitter_id = file["cases"][0]["submitter_id"]
		elif not (data_type in UNUSED_DATA_TYPE) and not ((data_type == "Copy Number Segment") and (workflow_type == "DNAcopy")) and not(experimental_strategy in UNUSED_EXPERIMENTAL_STRATEGY):
			sample_info = file["cases"][0]["samples"]
			if data_type in TUMOR_DATA_TYPE:
				for i in range(len(sample_info)):
					if sample_info[i]["tissue_type"] == "Tumor":	
						submitter_id = sample_info[i]["submitter_id"]
						break
			else:
				submitter_id = sample_info[0]["submitter_id"]
							
		if len(submitter_id) > 0:
			file = File(file_id, data_type, workflow_type, platform, experimental_strategy, submitter_id, data_set_list[0].project_id)
			missing_file_list = match_to_data_set(missing_file_list, data_set_list,file)
	return missing_file_list
			
				
		
		

def match_to_data_set(missing_file_list, data_set_list, file):
	'''
	This function matches an individual file to a dataset by checking if its metadata matches the requirements of the dataset.
	If a file does not match any, it goes into the missing_file_list.

	args:
		missing_file_list(list of File objects): list of files which do not belong in any of the potential datasets
		data_set_list(list of Data_set objects): A list of potential data sets in a project.
		file(File object): A File object that contains its relevant metadata.
	returns:
		missing_file_list(list of File objects): list of files which do not belong in any of the potential datasets
	'''
	missing = True
	for data_set in data_set_list:
		match_count = 0
		if (data_set.data_type == "") or (data_set.data_type.lower() == file.data_type.lower()):
			match_count+=1
		if (data_set.workflow_type =="") or (data_set.workflow_type.lower() == file.workflow_type.lower()):
			match_count+=1
		if (data_set.platform == "") or (data_set.platform.lower() == file.platform.lower()):
			match_count+=1
		if (data_set.experimental_strategy == "") or (data_set.experimental_strategy.lower() == file.experimental_strategy.lower()):
			match_count+=1
		
		
		if match_count == 4:
			
			data_set.files.append(file)
			missing = False
			return missing_file_list
			
	if missing:
		missing_file_list.append(file)		
		return missing_file_list



def remove_missing_duplicates(missing_file_list):
	'''
	This function removes duplicate files which have the same attributes. Note: When comparing two File objects, 
	it compares their attributes.

	args:
		missing_file_list(list of File objects): list of files which do not belong in any of the potential datasets 
		as well as those which are part of a dataset that doesn't exist in the Xena hub. 
	returns:
		list: set of file attributes. Each element is a unique file, with its metadata joined with a "/".
	'''
	unique_list = list(set(missing_file_list))
	for i in range(len(unique_list)):
		temp_list = [unique_list[i].data_type, unique_list[i].workflow_type, unique_list[i].platform, unique_list[i].experimental_strategy]
		temp_list = [s for s in temp_list if len(s) > 0]
		unique_list[i] = "/".join(temp_list)
	return unique_list



def prune_data_sets(data_set_list):
	'''
	This function removes any possible data_sets that have 0 elements, as they are not real datasets in the project.

	args:
		data_set_list(list of Data_set objects): A list of potential data sets in a project.
	returns:
		list of Data_set objects: A list of real data sets in a project.
	'''
	pruned_list = []
	for data_set in data_set_list:
		if len(data_set.files) > 0:
			pruned_list.append(data_set)
	return pruned_list


def xena_dataset(data_set):
	'''
	This function calls the xena API for a list of sample submitter IDs that belong to a specific dataset

	args:
		data_set(Data_set object): A data set that the function will use to call the Xena API with.
	returns:
		list: A list of submitter IDs of the samples of the dataset. 
	'''
	elements = [data_set.project_id, data_set.name,"tsv"]
	data_set_id = ".".join(elements)
	
	samples = xena.dataset_samples(HUB, data_set_id, None)

	return samples

def compare_datasets(samples, missing_file_list, data_set):
	'''
	This function checks if the dataset created from the GDC API call exists on the Xena hub. 
	If it doesn't exist, the files of the dataset are added to the missing_file_list.

	args:
		samples (list): A list of submitter IDs of the dataset returned by the xena API call.
		missing_file_list(list of File objects): list of files which do not belong in any of the potential datasets
		data_set(Data_set object): A data set that the function will use to get the submitter IDs of the files with.
	returns:
		boolean: This boolean returns True if the data set was missing.
		list of File objects: list of files which do not belong in any of the potential datasets 
		as well as those which are part of a dataset that doesn't exist in the Xena hub. 
	'''
	if len(samples) == 0:
		
		for file in data_set.files:
			missing_file_list.append(file)
		return True, missing_file_list
	return False, missing_file_list
	
def compare_samples(samples, data_set, project_id):
	'''
	This function compares the number of samples as well as the submitter IDs from the Xena hub and the GDCof a 
	specific data_set.

	args:
		samples (list): A list of submitter IDs of the dataset returned by the xena API call.
		data_set(Data_set object): A data set that the function will use to get the submitter IDs of the files with.
		project_id(string): The project ID for the project that you are comparing samples on.
	returns:
		boolean: This boolean returns True if the submitter IDs between the GDC and Xena are not the same.
		list: The second variable is a list which returns the list of submitter IDs missing from the Xena side.
		list: The third variable is a list which returns the list of submitter IDs missing from the GDC side.

	'''
	submitter_id_list = []
	if project_id == "BEATAML1.0-COHORT":
		for file in data_set.files:
			submitter_id_list.append(file.submitter_id[:-1])
	else:
		for file in data_set.files:
			submitter_id_list.append(file.submitter_id)

	
	submitter_id_list = set(submitter_id_list)
	samples = set(samples)
	missing_ids = submitter_id_list - samples
	missing_ids_GDC_dtype = samples - submitter_id_list
	missing_ids = list(missing_ids)
	missing_ids_GDC_dtype = list(missing_ids_GDC_dtype)
	logging.info(f"sample Count Xena: {len(samples)}")
	logging.info(f"sample Count GDC: {len(submitter_id_list)}")
	logging.info(f"sample Count that Xena is missing: {len(missing_ids)}")
	logging.info(f"sample Count that GDC is missing: {len(missing_ids_GDC_dtype)}")
	print(f"sample Count Xena: {len(samples)}")
	print(f"sample Count GDC: {len(submitter_id_list)}")
	print(f"sample Count that Xena is missing: {len(missing_ids)}")
	print(f"sample Count that GDC is missing: {len(missing_ids_GDC_dtype)}")

	if len(samples) == len(submitter_id_list):
		
		return False, missing_ids, missing_ids_GDC_dtype

	if len(samples) < len(submitter_id_list):
		
		return True, missing_ids, missing_ids_GDC_dtype
	print("ERROR: Xena dataset has more samples than GDC.")
	return True, missing_ids, missing_ids_GDC_dtype

def to_tsv(total_CMDS, total_MDS, total_CDSMS, total_DSMS, total_CMDT, total_MDT, project_id_list):
	'''
	This function creates a TSV file called "detect_missing_datasets.tsv" using the columns: "#_missing_datasets",
	 "missing_datasets", "#_datasets_with_wrong_sampleN", "datasets_with_wrong_sampleN", "#_missing_datatypes", 
	"missing_datatypes". The rows of this file are indexed by project ID. 

	args:
		Each element corresponds to a project for the below lists.
		total_CMDS(list): list of the number of missing datasets.
		total_MDS(list): list of the missing datasets.
		total_CDSMS(list): list of the number of datasets with incorrect sample count.
		total_DSMS(list): list of the datasets with incorrect sample count.
		total_CMDT(list): list of the number of datatypes missing from the Xena hub.
		total_MDT(list): list of the datatypes missing from the Xena hub.
		project_id_list(list): project IDs of projects that were tested when this script ran
	'''
	dataset_dict = {"#_missing_datasets": total_CMDS, 
		"missing_datasets": total_MDS, 
		"#_datasets_with_wrong_sampleN": total_CDSMS,
		"datasets_with_wrong_sampleN": total_DSMS,
		"#_missing_datatypes": total_CMDT,
		"missing_datatypes": total_MDT}
	dataset_df = pd.DataFrame(dataset_dict, index =project_id_list)
	dataset_df.to_csv("detect_missing_datasets.tsv", sep="\t", index=True)
	logging.info("\n")
	print("\n")
	print(dataset_df)
	print("\nDataset Dataframe saved to tsv file: detect_missing_datasets.tsv")

def to_logging(missing_submitter_ids):
	'''
	Adds to the logging file all of the sample submitter IDs that were missing from the Xena hub or GDC, organized 
	by project and dataset.

	args:
		missing_submitter_ids (dictionary): Contains all submitter IDs that were missing from the Xena hub or GDC.
	'''

	logging.info(json.dumps(missing_submitter_ids, indent=2))
	print("\nmissing_submitter_ids dictionary saved using logging module. File name: missing_submitter_ids.json")

def test_check(test_mode, test, project_id_list,missing_submitter_ids):
	'''
	This functions checks if the script is going to run one project or all projects.
	args: 
		test_mode(boolean): A boolean which is true when test mode is on
		test (string): The project ID of the project that the script is being run on.
		project_id_list(list): The list of project IDs the script will run on. 
		missing_submitter_ids(dict): This dictionary is where all submitter IDs missing from Xena will be organized.
	returns:
		list: List of all project IDs the script will run through
		dictionary: A dictionary with no values, but the keys correspond to projects.
	'''

	if test_mode:
		project_id_list = [test]
	else:
		project_id_list = project_request(PROJECT_FIELDS)

		project_id_list = format_to_list(project_id_list)
		
	missing_submitter_ids = {project_id: None for project_id in project_id_list}
	return project_id_list, missing_submitter_ids


class Data_set:
	# This class acts as a storage of metadata on a Dataset. Including its requirements, project Id, and the files that belong to it.
	

	def __init__(self,name, data_type, workflow_type, platform, experimental_strategy, project_id):
		self.data_type = data_type
		self.workflow_type = workflow_type
		self.platform = platform
		self.project_id = project_id
		self.name = name
		self.experimental_strategy = experimental_strategy
		self.files = []

	def add_file(file):
		files.append(file)
class File:
	# This class acts as a storage of metadata on a File, including its submitter_id, project_id, and metadata used to match a file to a dataset. 
	def __init__(self, file_id, data_type, workflow_type, platform,experimental_strategy ,submitter_id, project_id):
		self.data_type = data_type
		self.workflow_type = workflow_type
		self.platform = platform
		self.project_id = project_id
		self.file_id = file_id
		self.experimental_strategy = experimental_strategy
		self.submitter_id = submitter_id

	def __eq__(self, other):
		if not isinstance(other, File):
			return NotImplemented
		return self.data_type == other.data_type and self.workflow_type == other.workflow_type and self.platform == other.platform and self.experimental_strategy == other.experimental_strategy

	def __hash__(self):
		return hash((self.data_type, self.workflow_type, self.platform, self.experimental_strategy))



def test_project(project_id, 
	missing_submitter_ids, 
	total_CMDS,
	total_MDS,
	total_CDSMS,
	total_DSMS,
	total_CMDT,
	total_MDT):
	'''
	This function runs the file request to GDC and organizes them to correct datasets. Then it runs sample requests to Xena, and compares the submitter IDs of both.

	args:
		project_id (string): The project ID of the project that the function is running through.
		missing_submitter_ids (dictionary): submitter_ids that are missing from the Xena hub are saved here. Organized by project ID and dataset. 
		Used to form logging file

		Below variables used to form TSV file. Each element represents the data for one project
		total_MDS(list):Datasets missing from the Xena hub. 
		total_CMDS (list): # of datasets missing from the Xena hub.
		total_DSMS (list): Datasets that are missing samples.
		total_CDSMS (list): # of datasets that are missing samples.
		total_MDT (list): Datatypes that are missing from the Xena hub.
		total_CMDT (list): # of datatypes that are missing from the Xena hub. 

	'''
	logging.info(f"Project ID: {project_id} \n")
	print(f"Project ID: {project_id} \n")

	file_list = file_request(project_id, FILE_FIELDS)
	missing_file_list = []
	data_set_list = create_data_set(project_id)
	missing_file_list=  file_metadata(file_list, data_set_list, project_id, missing_file_list)
	data_set_list = prune_data_sets(data_set_list)
	logging.info(f"Number of samples which don't match any Xena data types: {len(missing_file_list)} \n")
	print(f"Number of samples which don't match any Xena data types: {len(missing_file_list)} \n")

	missing_data_sets = []
	data_sets_missing_sample = []
	missing_ids = {data_set.name: None for data_set in data_set_list}
	missing_ids_GDC = {data_set.name: None for data_set in data_set_list}


	for data_set in data_set_list:
		
		logging.info(f" \ndata type: {data_set.name}")
		print(f"\ndata type: {data_set.name}")

		samples = xena_dataset(data_set)

		isTrue ,missing_file_list = compare_datasets(samples, missing_file_list, data_set)
		
		if isTrue:
			missing_data_sets.append(data_set)
		isTrue, missing_ids_list, missing_ids_GDC_dtype = compare_samples(samples, data_set, project_id)
		if isTrue:
			data_sets_missing_sample.append(data_set)
		missing_ids[data_set.name] = missing_ids_list
		missing_ids_GDC[data_set.name] = missing_ids_GDC_dtype


	unique_list = []
	if len(missing_file_list) != 0:
		unique_list = remove_missing_duplicates(missing_file_list)
	for value in range(len(missing_file_list)):
		missing_file_list[value] = missing_file_list[value].submitter_id
	missing_file_list = list(set(missing_file_list))
	missing_ids["not_in_datasets"] = missing_file_list
	missing_ids["not_in_GDC"] = missing_ids_GDC
	missing_submitter_ids[project_id] = missing_ids

	logging.info(f"\nnumber of unique samples missing from datasets: {len(missing_file_list)}")
	logging.info(f"missing data types: {unique_list}")
	print(f"\nnumber of unique samples missing from datasets: {len(missing_file_list)}")
	print(f"missing data types: {unique_list}")
	
	names_MDS = []
	names_DSMS = []
	logging.info("\nmissing_datasets: ")
	print("\nmissing_datasets: ")
	for item in missing_data_sets:
		logging.info(item.name)
		print(item.name)
		names_MDS.append(item.name)
	logging.info("\ndatasets_missing_samples: ")
	print("\ndatasets_missing_samples: ")
	for item in data_sets_missing_sample:
		logging.info(item.name)
		print(item.name)
		names_DSMS.append(item.name)



	total_MDS.append(names_MDS)
	total_CMDS.append(len(names_MDS))
	total_DSMS.append(names_DSMS)
	total_CDSMS.append(len(names_DSMS))
	total_CMDT.append(len(unique_list))
	total_MDT.append(unique_list)


#Script runs in below lines.
project_id_list, missing_submitter_ids = test_check(test_mode, test, project_id_list,missing_submitter_ids)
for project_id in project_id_list:
	
	test_project(project_id, missing_submitter_ids, total_CMDS,total_MDS,total_CDSMS,total_DSMS,total_CMDT,total_MDT)

to_tsv(total_CMDS, total_MDS, total_CDSMS, total_DSMS, total_CMDT, total_MDT, project_id_list)
to_logging(missing_submitter_ids)
end_time = time.time()
runtime = end_time - start_time
print(f"\nRuntime: {runtime} seconds")
