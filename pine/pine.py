"""
PINE - a tool for visualizing protein-protein interactions.

Given a list of protein IDs, constructs in Cytoscape
  1) an interaction network among all proteins in the list
  2) pathway network of proteins from the list
optionally may include their corresponding PTM modifications, fold change, p-values or categories

"""

import sys
def eprint(*args, **kwargs):
  ''' Print to stderr instead of stdout '''
  print(*args, file=sys.stderr, **kwargs)

from pinepy2cytoscape.data.cyrest_client import CyRestClient
try:
  import requests
except ImportError:
  eprint("ImportError: Please make sure you are using Python3")
from urllib.parse import quote
import datetime as dt
from datetime import datetime
import urllib
import urllib.request as urllib2
import pandas as pd
import json
from io import StringIO
import csv
import getopt 
import os
from os import path
import time
import re
import logging
import gzip
import json
import math
import igraph
import traceback
import collections
import subprocess
import warnings
from collections import Counter
import socket
warnings.simplefilter(action='ignore', category=FutureWarning)

__author__ = "Niveda Sundararaman, James Go and Vidya Venkatraman"
__credits__ = ["Niveda Sundararaman", "James Go", "Vidya Venkatraman"]
__license__ = "Apache-2.0"
__version__ = "0.1.2"
__maintainer__ = "Niveda Sundararaman"
__email__ = "GroupHeartBioinformaticsSupport@cshs.org"
__status__ = "Production"

CYREST_PORTS = [8012, 8013]
CYREST_URL = None
CYREST_PORT = None

class CytoscapeError(Exception):
  ''' Exception for errors when making requests to Cytoscape '''
  def __init__(self, message):
    super().__init__(message)

class PineError(Exception):
  def __init__(self, message):
    super().__init__(message)

def setup_logger(name, log_file, level=logging.DEBUG, with_stdout=False):
  ''' Define logging criteria to generate logs to include warnings or reasons for dropping protein from analysis '''
  handler = logging.FileHandler(log_file,mode='w')
  logger = logging.getLogger(name)
  logger.setLevel(level)
  logger.addHandler(handler)
  if with_stdout:
    logger.addHandler(logging.StreamHandler(sys.stdout))
  return logger

def port_in_use(port):
  ''' Check if can connect to a localhost port. Credit: https://stackoverflow.com/a/52872579 '''
  with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    bound = s.connect_ex(('localhost', port)) == 0
    s.close()
    return bound

def request_retry(url, protocol, headers=None, data=None, json=None, timeout=300, timeout_interval=5, request_timeout=300):
  '''
  API request calls to implement app calls, network construction and style changes within Cytoscape
  Retries requests until a 200 status is returned or the timeout expires
  timeout is the maximum time of the function
  request_timeout is the maximum time for a single request.
  If request_timeout > timeout, then request_timeout will be the maximum time for the function
  '''
  if protocol == "GET":
    func = requests.get
  elif protocol == "PUT":
    func = requests.put
  elif protocol == "POST":
    func = requests.post
  elif protocol == "DELETE":
    func = requests.delete

  kwargs = {}
  kwargs["timeout"] = request_timeout
  if headers:
    kwargs["headers"] = headers
  if data:
    kwargs["data"] = data
  if json:
    kwargs["json"] = json

  timer = time.time()
  timer_end = timer + timeout
  while timer < timer_end:
    try:
      res = func(url, **kwargs)
      res.raise_for_status()
      return res
    except requests.exceptions.RequestException as e:
      time.sleep(timeout_interval)
      timer = time.time()
      last_exception = e
  raise CytoscapeError(str(last_exception))

def input_failure_comment(uniprot_query, merged_out_dict, input_id, comment):
  ''' Label an ID from the input with a failure comment '''
  uniprot_query[input_id] = {
    'Uniprot': "NA",
    'Primary': "NA",
    'Synonym': "NA",
    'Organism': "NA",
    'Reviewed': "NA",
    'Date_modified': "NA"
  }

  merged_out_dict[input_id] = {
    'Primary': '',
    'CommentGene': comment,
    'String': '',
    'Genemania': '',
    'ClueGO': ''
  }

def db_handling(db_file):
  ''' Function that reads FASTA input provided by user and associates proteinIDs to their sequence later used for PTM site determination '''
  database_dict = {}
  proteinID = ""
  seq = ""
  if db_file != "":
    with open(db_file) as database:
      for line in database:
        if line.startswith(">"):
          if "DECOY" in line:
            if proteinID:
              if '|' in proteinID:
                proteinID = proteinID.split('|')[1]
       
              if proteinID not in database_dict:
                database_dict[proteinID] = [seq]
              else:
                database_dict[proteinID].append(seq)
              proteinID = ""
              seq = ""
              continue
          else:
            if proteinID:
              if '|' in proteinID:
                proteinID = proteinID.split('|')[1]
              if proteinID not in database_dict:
                database_dict[proteinID] = [seq]
              else:
                database_dict[proteinID].append(seq)
          
            proteinID = line.strip(">\n\r").split(" ")[0]
            seq = ""
        
        else:           
          line_seq = line.strip("\n\r")
          if " " in line_seq:
            raise PineError("Fasta file sequences should not have spaces")
          seq += line_seq
            
      if proteinID:
        if '|' in proteinID:
          proteinID = proteinID.split('|')[1]
      
        if proteinID not in database_dict:
          database_dict[proteinID] = [seq]
        else:
          database_dict[proteinID].append(seq)
     
        proteinID = line.strip(">\n").split(" ")[0]
  return(database_dict)

def find(seq,database):
  ''' Function to obtain start and end positions of the peptide (without including the modifications) within the FASTA sequence '''
  seq = seq.replace("[","\[")
  seq = seq.replace("]","\]")
  seq = seq.replace("(","\(")
  seq = seq.replace(")","\)")
  seq = seq.replace("{","\{")
  seq = seq.replace("}","\}")
  start_index_list = []
  for m in re.finditer(seq, database):
    start_index_list.append(m.start())
  return start_index_list
           
def find_mod(seq): 
  ''' Function that uses FASTA input to find all the modification positions in the peptide '''
  result_dict = {} 
  i = 0
  index = 0
  for letter in seq:
    if letter == "[" or letter == "(" or letter == "{":
      start_point_site = i
      start_point = index
    if letter == "]" or letter == ")" or letter == "}":
      end_point = index
      amino_acid = seq[start_point-1:end_point+1]
      if amino_acid in result_dict:
        #if result_dict.has_key(amino_acid):
        result_dict[amino_acid].append(start_point_site-1)
      else:
        result_dict[amino_acid] = [start_point_site-1]
      i = i -(end_point - start_point) - 1
    i += 1
    index += 1
  return result_dict

def latest_ontologies(all_ontologies):
  """
  Takes a dict of onotologies from ClueGO and only returns the ontologies from the latest dataset
  """
  latest_dates = {}
  for ont_id in all_ontologies:
    ont = all_ontologies[ont_id]
    ont_type = None
    ont_date = None
    for key_val in ont.split(", "):
      if key_val.count("=") != 1:
        continue
      key, val = key_val.split("=")
      if key == "date":
        try:
          ont_date = datetime.strptime(val, "%d.%m.%Y")
        except ValueError:
          break
      elif key == "type":
        ont_type = val.lower()
    if ont_type is None or ont_date is None:
      continue
    if ont_type not in latest_dates or ont_date > latest_dates[ont_type]["date"]:
      latest_dates[ont_type] = {
        "date": ont_date,
        "ids": [ont_id],
      }
    elif latest_dates[ont_type]["date"] == ont_date:
      latest_dates[ont_type]["ids"].append(ont_id)

  final_onts = {}
  for ld in latest_dates.values():
    for idee in ld["ids"]:
      final_onts[idee] = all_ontologies[idee]
  return final_onts


def preprocessing(inp, type, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out, database_dict, include_list, db_file, enzyme, path_to_new_dir, logging_file, cy_fc_cutoff, cy_pval_cutoff, exclude_ambi, cy_settings_file):
  '''
  Reads input file and processes data based on type of analysis that will be performed
  Checks for all required columns and errors out the wrong analysis was chosen
  '''
  is_prot_col = False
  is_FC = False
  is_pval = False
  is_cat = False
  is_label_col = False
  is_pep_col = False
  each_protein_list = []
  prot_list = {}
  site_info_dict = {}
  max_FC_len = 0
  unique_labels = []
  each_category = []
  initial_length = 0
  repeat_prot_ids = []
  retain_prot_ids = []
  repeat_prot_ids_2 = {}
  list_of_duplicates = []
  to_return_unique_protids_length = 0
  ambigious_sites = {}
  ambigious_sites_ptms = {}
  pep_to_mod = {}
  duplicate_inc_ptm_proteins = []
  duplicate_inc_ptm_proteins_set = set()
  duplicate_ptm_proteins = []
  mapping_multiple_regions = {}
  pep_not_in_fasta = {}
  dropped_invalid_fc_pval = {}
  dropped_cutoff_fc_pval = {}
  dropped_invalid_site = {}
  dict_of_picked_pep = {}
  ambi_pep_count = 0
  initial_query_prot_count = 0
  initial_query_pep_count = 0
  unique_prot_pep = {}
  ctr = 0
  mult_mods_of_int = True  #False
  unique_unimods = []
  pep_to_prot_dict = {}
  dup_pep_list = []
  raw_category_set = set()
  raw_label_set = set()
  count_dup_pep = []
  collect_dup_pep = {}
  
  try:
    with open(inp,'r') as csv_file:
      '''
      Read input and collect columns based on type of analysis chosen
      Types:1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = Category; 5 = Singlefc-ptm; 6 = Multifc-ptm
      '''
      input_file = csv.reader(csv_file, delimiter=',')
      line_count = 0
      for row in input_file:
        skip_val = False    
        if line_count == 0:
          header = row
          # Get input headers
          for i in range(len(row)):                                                                               
            if "proteinid" in row[i].lower():
              protein = i
              is_prot_col = True
            elif "fc" in row[i].lower():
              FC=i
              is_FC = True  
            elif "adj.pval" in row[i].lower():
              pval = i
              is_pval = True
            elif "pval" in row[i].lower() and not is_pval:
              pval=i
              is_pval = True  
            elif "category" in row[i].lower():
              cat = i
              is_cat = True
            elif "label" in row[i].lower():
              label = i  
              is_label_col = True
            elif "peptide" in row[i].lower():
              peptide_col = i
              is_pep_col = True         
          
          # If required columns are missing from input type, print appropriate error messages
          if type == "1":
            if not (is_prot_col and is_FC): 
              eprint("Error: Columns 'ProteinID' and 'FC' required for singleFC run type")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_cat or is_label_col or is_pep_col):
              eprint("Error: SingleFC run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
          
          elif type == "2":
            if not (is_prot_col and is_FC and is_label_col):
              eprint("Error: Columns 'ProteinID', 'FC' and 'Label' required for multiFC run type")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_cat or is_pep_col):
              eprint("Error: MultiFC run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
          
          elif type == "3":
            if not (is_prot_col):
              eprint("Error: Columns 'ProteinID' required for list only run type")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_cat or is_pep_col or is_FC or is_label_col or is_pval):
              eprint("Error: List Only run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
          
          elif type == "4":
            if not (is_prot_col and is_cat):
              eprint("Error: Columns 'ProteinID' and 'Category' required for category run type")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_pep_col or is_FC or is_label_col or is_pval):
              eprint("Error: Category run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
              
          elif type == "5":
            if not (is_prot_col and is_pep_col and is_FC):
              eprint("Error: Columns 'ProteinID', 'Peptide' and 'FC' required for singlefc-ptm")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_cat or is_label_col):
              eprint("Error: Singlefc-ptm run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
              
          elif type == "6":
            if not (is_prot_col and is_pep_col and is_FC and is_label_col):
              eprint("Error: Columns 'ProteinID', 'Peptide', 'FC' and 'Label' required for multifc-ptm")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if (is_cat):
              eprint("Error: Multifc-ptm run chosen but other headers found. Please check that the run selected is correct")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
              
        else: 
          str_each_row = ",".join(row)      
          str_each_row = str_each_row.replace(",","")
          if not str_each_row:
            continue
          if len(header) != len(row):
            eprint("Error: Number of columns does not match number of header columns at line " + str(line_count + 1))
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
            sys.exit(1) 
          if is_cat:
            if row[cat].strip() == "":
              eprint("Error: Blank category in line " + str(line_count+1))
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1) 
            raw_category_set.add(row[cat])
          if is_label_col:
            if row[label].strip() == "":
              eprint("Error: Blank label in line " + str(line_count+1))
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1) 
            raw_label_set.add(row[label])
          
          # Check if column marked as proteinID in the input has valid Uniprot IDs
          check_match = re.match('^[A-Za-z0-9\-]+$', row[protein])
          if not check_match:
            eprint("Error: Invalid proteinID: " + row[protein] + " in line " + str(line_count+1))
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
            sys.exit(1) 
          
          # Check if column marked as peptide in the input has valid peptide terms. Must only include peptide sequences with modifications and its corresponding unimod number enclosed in brackets. Ex: SEDVLAQS[+80]PLPK
          if type == "5" or type == "6":            
            check_match = re.match('^[A-Za-z]{1,}([\[\(\{]{1}[^\[\(\{\)\]\}]{1,}[\]\}\)]{1}){0,}[A-Z]{0,}', row[peptide_col])
            if not check_match:
              eprint("Error: Invalid peptide: " + row[peptide_col] + " in line " + str(line_count+1))
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
              sys.exit(1)
            if row[protein] not in unique_prot_pep:
              unique_prot_pep.update({row[protein]:[row[peptide_col]]})
              initial_query_pep_count += 1
            else:
              if row[peptide_col] not in unique_prot_pep[row[protein]]:
                unique_prot_pep[row[protein]].append(row[peptide_col])
                initial_query_pep_count += 1
            if not is_pval:
              get_pval = 1.0 # If pvalue is not provided, default set to 1.0
            else:
              get_pval = row[pval]
            get_fc_val = row[FC]
                        
          if type == "1" or type == "2":
            initial_query_pep_count += 1
            skip_val = False
            if not is_pval:
              get_pval = 1.0
            else:
              # Discard all non-numeric fold changes and pvalues
              try:
                float(row[pval])
                if math.isinf(float(row[pval])):
                  skip_val = True 
                else:
                  get_pval = float(row[pval])
              except ValueError:
                skip_val = True
                
            if is_FC:
              try:
                float(row[FC])
                if math.isinf(float(row[FC])):
                  skip_val = True
                else:
                  get_fc_val = float(row[FC])
              except ValueError:
                skip_val = True
                
            if skip_val:
              if row[protein] not in each_protein_list:          
                dropped_invalid_fc_pval.update({row[protein]:None})
              continue
            elif not skip_val and row[protein] in dropped_invalid_fc_pval:
              del dropped_invalid_fc_pval[row[protein]]           
                    
          if type == "5" or type == "6":        
            all_mods_for_prot = []
            if include_list:
              modInSeq_dict = {}
              peptide = row[peptide_col]
              protein_list_id = row[protein]
              modInSeq_all_dict = find_mod(peptide)
              combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
              peptide_sub = re.sub(combined_pat, '', peptide)
              peptide_sub = peptide_sub.strip('"')		           
              if protein_list_id in database_dict:
                for each_seq_in_db_dict in database_dict[protein_list_id]:
                  # Check for position of peptide in FASTA sequence
                  seqInDatabase_list = find(peptide_sub,each_seq_in_db_dict)
                  if len(seqInDatabase_list) > 1:
                    # If peptide maps to multiple regions in FASTA, then either pick first option or drop peptide based on user's choice to retain or drop ambiguity
                    if protein_list_id not in mapping_multiple_regions:   
                      key_val = protein_list_id + "-" + peptide
                      each_pos_store_as_mult = []
                      for each_pos in seqInDatabase_list:
                        store_as_mult = []
                        for each_mod in modInSeq_all_dict:
                          for each_include_list in include_list:
                            combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
                            each_mod_key = re.sub(combined_pat, '', each_mod) 
                            each_include_list = re.sub(combined_pat, '', each_include_list)
                            if each_include_list == each_mod_key:
                              for each_pos_in_list in modInSeq_all_dict[each_mod]:
                                store_as_mult.append(each_include_list + str(each_pos_in_list+each_pos+1))
                        if store_as_mult:
                          each_pos_store_as_mult.append('/'.join(store_as_mult))
                      if key_val not in mapping_multiple_regions and each_pos_store_as_mult:
                        mapping_multiple_regions.update({key_val:each_pos_store_as_mult})
                        if protein_list_id not in ambigious_sites:
                          ambigious_sites.update({protein_list_id:[each_pos_store_as_mult[0]]})
                        else:
                          ambigious_sites[protein_list_id].append(each_pos_store_as_mult[0])
                  # Track peptides not in FASTA sequence for later logging
                  if len(seqInDatabase_list) == 0:
                    seqInDatabase = -1
                    if protein_list_id not in pep_not_in_fasta:   
                      pep_not_in_fasta.update({protein_list_id:[peptide_sub]})
                    else:
                      if peptide_sub not in pep_not_in_fasta[protein_list_id]:
                        pep_not_in_fasta[protein_list_id].append(peptide_sub)
                  elif len(seqInDatabase_list) == 1:                        
                    seqInDatabase = seqInDatabase_list[0]
                  else:
                    if exclude_ambi:
                      seqInDatabase = "Ambiguous"
                    else:
                      seqInDatabase = seqInDatabase_list[0]                 
                  break
                  
              else:
                # If protein not in FASTA, then print error message. This could indicate that the user has provided a different FASTA than what was used prior to PINE analysis.
                eprint("Error: ProteinID '" + str(protein_list_id) + "' not found in " + str(db_file) + ". Check fasta database file provided.")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                sys.exit(1)         
                
              if seqInDatabase!=-1 and seqInDatabase!="Ambiguous":           
                for k in include_list:
                  for key,value in modInSeq_all_dict.items():                   
                    combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
                    #key = re.sub('\+','',key)                               
                    #k = re.sub('\+','',k)
                    key_match = re.search(combined_pat,key)
                    k_match = re.search(combined_pat,k)
                    if not (key_match and k_match):
                      key1 = re.sub(combined_pat, '', key)
                      k1 = re.sub(combined_pat,'',k)
                    else:
                      key1 = key
                      k1 = k                   
                    if k1.lower() in key1.lower(): 
                      for each_val in value:
                        val = int(each_val)+int(seqInDatabase)+1
                        match_unimod = re.findall(combined_pat, key)
                        remove_brackets = r'|'.join(('\[', '\{', '\(', '\]', '\}', '\)'))
                        key_with_unimod = re.sub(combined_pat,'',k1) + "{" + re.sub(remove_brackets, '', match_unimod[0]) + "}" 
                        #if match_unimod[0] not in unique_unimods:
                          #unique_unimods.append(match_unimod[0])                      
                        if key_with_unimod in modInSeq_dict:                      
                          modInSeq_dict[key_with_unimod].append(val)                          
                        else:
                          modInSeq_dict[key_with_unimod] = [val]
                        all_mods_for_prot.append(key)         
                         
              # If a protein ID has multiple modification sites, separator "/" is used for representation 
              sites = ""
              sites_list = []
              for key,value in modInSeq_dict.items():
                for each_val in value:
                  sites += key + str(each_val) + "/"
              sites = sites.strip("/")
              if sites:            
                if type == "5":
                  if (protein_list_id, sites, peptide) in duplicate_inc_ptm_proteins_set:
                    duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, get_fc_val, get_pval))
                    continue
                  max_FC_len = 1
                  # Map proteinIDs and sites to their corresponding fold change and pvalues
                  if is_prot_col and is_FC:
                    if protein_list_id not in site_info_dict:
                      site_info_dict[protein_list_id] = {}
                      site_info_dict[protein_list_id][sites] = {}
                      site_info_dict[protein_list_id][sites].update({peptide:[[get_fc_val],[get_pval],all_mods_for_prot]})
                    elif sites not in site_info_dict[protein_list_id]:
                      site_info_dict[protein_list_id][sites] = {}
                      site_info_dict[protein_list_id][sites].update({peptide:[[get_fc_val],[get_pval],all_mods_for_prot]})
                    elif not peptide in site_info_dict[protein_list_id][sites]:
                      site_info_dict[protein_list_id][sites].update({peptide:[[get_fc_val],[get_pval],all_mods_for_prot]})
                    else:
                      peptide_list = site_info_dict[protein_list_id][sites][peptide]
                      if peptide_list[0][0] != get_fc_val or peptide_list[1][0] != get_pval:
                        duplicate_inc_ptm_proteins_set.add((protein_list_id, sites, peptide))
                        duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, get_fc_val, get_pval))
                        duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, peptide_list[0][0], peptide_list[1][0]))
                        del site_info_dict[protein_list_id][sites][peptide]
                        if len(site_info_dict[protein_list_id][sites]) == 0:
                          del site_info_dict[protein_list_id][sites]
                          if len(site_info_dict[protein_list_id]) == 0:
                            del site_info_dict[protein_list_id]
                            ctr +=1 
                      else:
                        duplicate_ptm_proteins.append((protein_list_id, sites, peptide, get_fc_val, get_pval))
                      continue
                  else:
                    eprint("Error: Required columns- ProteinID, Peptide, FC")
                    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                    sys.exit(1)
                
                if type == "6":  
                  if is_prot_col and is_FC and is_label_col:
                    if (protein_list_id, sites, peptide, row[label]) in duplicate_inc_ptm_proteins_set:
                      duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, row[label], get_fc_val, get_pval))
                      continue               
                    if protein_list_id not in site_info_dict:
                      site_info_dict[protein_list_id] = {}
                      site_info_dict[protein_list_id][sites] = {}
                      site_info_dict[protein_list_id][sites].update({peptide: [[get_fc_val],[get_pval],all_mods_for_prot,[row[label]]] })
                    elif sites not in site_info_dict[protein_list_id]:
                      site_info_dict[protein_list_id][sites] = {}
                      site_info_dict[protein_list_id][sites].update({peptide: [[get_fc_val],[get_pval],all_mods_for_prot,[row[label]]] })
                    elif peptide not in site_info_dict[protein_list_id][sites]:
                      site_info_dict[protein_list_id][sites].update({peptide: [[get_fc_val],[get_pval],all_mods_for_prot,[row[label]]] })
                    elif row[label] not in site_info_dict[protein_list_id][sites][peptide][3]:
                      site_info_dict[protein_list_id][sites][peptide][0].append(get_fc_val)
                      site_info_dict[protein_list_id][sites][peptide][1].append(get_pval)
                      site_info_dict[protein_list_id][sites][peptide][3].append(row[label])
                    else:
                      peptide_list = site_info_dict[protein_list_id][sites][peptide]
                      ix = peptide_list[3].index(row[label])
                      if peptide_list[0][ix] != get_fc_val or peptide_list[1][ix] != get_pval:
                        duplicate_inc_ptm_proteins_set.add((protein_list_id, sites, peptide, row[label]))
                        duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, row[label], get_fc_val, get_pval))
                        duplicate_inc_ptm_proteins.append((protein_list_id, sites, peptide, row[label], peptide_list[0][ix], peptide_list[1][ix]))
                        del peptide_list[0][ix]
                        del peptide_list[1][ix]
                        del peptide_list[3][ix]
                        if len(peptide_list[0]) == 0:
                          del site_info_dict[protein_list_id][sites][peptide]
                          if len(site_info_dict[protein_list_id][sites]) == 0:
                            del site_info_dict[protein_list_id][sites]
                            if len(site_info_dict[protein_list_id]) == 0:
                              del site_info_dict[protein_list_id]
                              ctr +=1 
                      else:
                        duplicate_ptm_proteins.append((protein_list_id, sites, peptide, row[label], get_fc_val, get_pval))
                      continue
                  else:
                    eprint("Error: Required columns- ProteinID, Peptide, FC, Label")
                    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                    sys.exit(1)
              else:
                if seqInDatabase!="Ambiguous": 
                  # If no modifications of interest are found in the peptide sequence, assign site = NA              
                  sites = "NA"
                  if type == "6":
                    if protein_list_id not in dropped_invalid_site:
                      dropped_invalid_site[protein_list_id] = {}
                      dropped_invalid_site[protein_list_id].update({"PeptideandLabel":[peptide + "-" + row[label]], "Site":[sites]})
                    else:
                      dropped_invalid_site[row[protein]]["PeptideandLabel"].append(peptide + "-" + row[label])
                      dropped_invalid_site[row[protein]]["Site"].append(sites)
                  elif type == "5":
                    if protein_list_id not in dropped_invalid_site:
                      dropped_invalid_site[protein_list_id] = {}
                      dropped_invalid_site[protein_list_id].update({"PeptideandLabel":[peptide], "Site":[sites]})
                    else:
                      dropped_invalid_site[protein_list_id]["PeptideandLabel"].append(peptide)
                      dropped_invalid_site[protein_list_id]["Site"].append(sites)        
              
          if row[protein] not in each_protein_list:
            each_protein_list.append(row[protein])
            if type == "1" or type == "2" or type == "3":          
              if type == "1":
                if is_prot_col and is_FC:           
                  prot_list.update({row[protein]:[[float(get_fc_val)],[get_pval]]})
                  max_FC_len = 1
                else:
                  eprint("Error: Required columns- ProteinID, FC")
                  remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                  sys.exit(1)
                  
              elif type == "2":
                if is_prot_col and is_FC and is_label_col:
                  prot_list[row[protein]] = {}
                  prot_list[row[protein]].update({row[label]:[float(get_fc_val),get_pval]})
                  if row[label] not in unique_labels:
                    unique_labels.append(row[label])                
                else:
                  eprint("Error: Required columns- ProteinID, FC, Label")
                  remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                  sys.exit(1)
                  
              elif type == "3":
                if is_prot_col:
                  continue
                else:
                  eprint("Error: Required columns- ProteinID")
                  remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                  sys.exit(1)
              
            elif type == "4":
              if (not is_prot_col) or (not is_cat):
                eprint("Error: Required columns- ProteinID, Category")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
                sys.exit(1)
                
              if row[cat] not in each_category:
                each_category.append(row[cat])
              prot_list.update({row[protein]:[row[cat]]})
                          
          else:  
            if type == "3":
              each_protein_list.append(row[protein])
              repeat_prot_ids.append(row[protein])
            
            elif type == "1":
              each_protein_list.append(row[protein])    
              if float(get_fc_val) == prot_list[row[protein]][0][0] and prot_list[row[protein]][1][0] == float(get_pval):
                retain_prot_ids.append(row[protein])
              else:                          
                repeat_prot_ids.append(row[protein])
              
            elif type == "4":
              each_protein_list.append(row[protein]) 
              if row[cat] not in each_category:
                each_category.append(row[cat])
              if row[cat] in prot_list[row[protein]]:
                retain_prot_ids.append(row[protein])
              else:
                prot_list[row[protein]].append(row[cat])
            
            elif type == "2":
              each_protein_list.append(row[protein])          
              if row[label] not in unique_labels:
                unique_labels.append(row[label])
              if row[label] in prot_list[row[protein]]:
                if prot_list[row[protein]][row[label]][0] == float(get_fc_val) and prot_list[row[protein]][row[label]][1] == float(get_pval):
                  retain_prot_ids.append(row[protein])
                else:
                  if row[protein] not in repeat_prot_ids_2:
                    repeat_prot_ids_2.update({row[protein]:[row[label]]})
                  else:
                    repeat_prot_ids_2[row[protein]].append(row[label])
              else:
                prot_list[row[protein]].update({row[label]:[float(get_fc_val),get_pval]})
        line_count+=1 
        
  except FileNotFoundError:
    eprint("Error: Input file is missing")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  except (UnicodeDecodeError, IndexError):
    eprint("Error: Input file must be in CSV (comma separated value) format")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  # Remove duplicate peptides across proteinIDs within the input
  if type == "5" or type == "6":
    mult_pep_to_prot_dict = {}
    combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
    
    if exclude_ambi:
      for each_key in mapping_multiple_regions:
        each_prot = each_key.split("-")[0]
        each_pep = each_key.split("-")[1]
        each_site =  '|'.join(mapping_multiple_regions[each_key])
        each_pep_seq_only = re.sub(combined_pat,'',each_pep)      
        if each_pep_seq_only in mult_pep_to_prot_dict:
          mult_pep_to_prot_dict[each_pep_seq_only]["Protein"].append(each_prot)
          mult_pep_to_prot_dict[each_pep_seq_only]["Site"].append(each_site)
          mult_pep_to_prot_dict[each_pep_seq_only]["Peptide"].append(each_pep)
        else:
          mult_pep_to_prot_dict[each_pep_seq_only] = {}
          mult_pep_to_prot_dict[each_pep_seq_only].update({"Protein":[each_prot], "Site":[each_site], "Peptide":[each_pep]})
      
    for prot_id in site_info_dict:
      bool = 0
      for site in site_info_dict[prot_id]:
        for peptide in site_info_dict[prot_id][site]:
          
          pep_seq_only = re.sub(combined_pat,'',peptide)
          if pep_seq_only in pep_to_prot_dict:
            pep_to_prot_dict[pep_seq_only]["Protein"].append(prot_id)
            pep_to_prot_dict[pep_seq_only]["Site"].append(site)
            pep_to_prot_dict[pep_seq_only]["Peptide"].append(peptide)           
          else:
            pep_to_prot_dict[pep_seq_only] = {}
            pep_to_prot_dict[pep_seq_only].update({"Protein":[prot_id], "Site":[site], "Peptide":[peptide]})
    
    for each_pep in pep_to_prot_dict:
      if len(list(set(pep_to_prot_dict[each_pep]["Protein"]))) > 1:
        count_dup_pep.append(each_pep)
        for each_prot, each_peptide, each_site in zip(pep_to_prot_dict[each_pep]["Protein"], pep_to_prot_dict[each_pep]["Peptide"], pep_to_prot_dict[each_pep]["Site"]):
          if each_prot not in collect_dup_pep:
            collect_dup_pep[each_prot] = {}
            collect_dup_pep[each_prot].update({"Peptide":[each_peptide], "Site":[each_site]})
          else:
            if each_peptide not in collect_dup_pep[each_prot]["Peptide"]:
              collect_dup_pep[each_prot]["Peptide"].append(each_peptide)
              collect_dup_pep[each_prot]["Site"].append(each_site)              
          del site_info_dict[each_prot][each_site][each_peptide]
          if len(site_info_dict[each_prot][each_site]) == 0:
            del site_info_dict[each_prot][each_site]
            if len(site_info_dict[each_prot]) == 0:
              del site_info_dict[each_prot]
              
      if each_pep in mult_pep_to_prot_dict:
        count_dup_pep.append(each_pep)
        for each_prot, each_peptide, each_site in zip(pep_to_prot_dict[each_pep]["Protein"], pep_to_prot_dict[each_pep]["Peptide"], pep_to_prot_dict[each_pep]["Site"]):
          if each_prot in site_info_dict:
            if each_site in site_info_dict[each_prot]:
              if each_peptide in site_info_dict[each_prot][each_site]:
                if each_prot not in collect_dup_pep:
                  collect_dup_pep[each_prot] = {}
                  collect_dup_pep[each_prot].update({"Peptide":[each_peptide], "Site":[each_site]})
                else:
                  if each_peptide not in collect_dup_pep[each_prot]["Peptide"]:
                    collect_dup_pep[each_prot]["Peptide"].append(each_peptide)
                    collect_dup_pep[each_prot]["Site"].append(each_site) 
                del site_info_dict[each_prot][each_site][each_peptide]
                if len(site_info_dict[each_prot][each_site]) == 0:
                  del site_info_dict[each_prot][each_site]
                  if len(site_info_dict[each_prot]) == 0:
                    del site_info_dict[each_prot]
                    
        for each_prot, each_peptide in zip(mult_pep_to_prot_dict[each_pep]["Protein"], mult_pep_to_prot_dict[each_pep]["Peptide"]):
          each_site = '|'.join(mapping_multiple_regions[each_prot+"-"+each_peptide])
          if each_prot not in collect_dup_pep:
            collect_dup_pep[each_prot] = {}
            collect_dup_pep[each_prot].update({"Peptide":[each_peptide], "Site":[each_site]})
          else:
            if each_peptide not in collect_dup_pep[each_prot]["Peptide"]:
              collect_dup_pep[each_prot]["Peptide"].append(each_peptide)
              collect_dup_pep[each_prot]["Site"].append(each_site) 
          del mapping_multiple_regions[each_prot+"-"+each_peptide]        

  if type == "6":
    # need to calculate unique labels outside loop since inconsistent labels are deleted
    unique_labels = set()
    for prot_id in site_info_dict:
      for site in site_info_dict[prot_id]:
        for peptide in site_info_dict[prot_id][site]:
          for label in site_info_dict[prot_id][site][peptide][3]:
            unique_labels.add(label)
    unique_labels = list(unique_labels)
  
  if len(unique_labels) > 6:
    eprint("Error: Number of unique labels should not exceed 6")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
    
  if len(each_category) > 6:
    eprint("Error: Number of categories should not exceed 6")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)

  if type == "4" and len(raw_category_set) < 2:
    eprint("Error: Input must contain at least 2 unique categories")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)

  if (type == "2" or type == "6") and len(raw_label_set) < 2:
    eprint("Error: Input must contain at least 2 unique labels")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)

  merged_out_dict_2 = {}
  if cy_debug:
    if type == "5" or type == "6":
      countpep = 0
      for getprot in site_info_dict:
        for getsite,getpep in site_info_dict[getprot].items():
          countpep += len(list(getpep.keys()))          
      initial_query_prot_count = len(list(site_info_dict.keys())) + len([x for x in dropped_invalid_fc_pval if x not in site_info_dict]) + len([x for x in dropped_invalid_site if x not in site_info_dict]) + ctr
      if exclude_ambi:
        initial_query_prot_count += len(mapping_multiple_regions)
      unique_peptides = set()
      for peps in unique_prot_pep.values():
        for pep in peps:
          unique_peptides.add(pep)
      logging.debug(f"Initial query: {len(unique_prot_pep)} unique protein IDs, {len(unique_peptides)} unique peptides")
      if dropped_invalid_site:
        site_len = 0
        pep_len = 0
        warning = []
        i = 0
        for each_key in dropped_invalid_site:
          warning_str = each_key
          if dropped_invalid_site[each_key] == None:
            continue
            
          warning.append(warning_str)
          
          for get_each_pep,get_each_site in zip(dropped_invalid_site[each_key]["PeptideandLabel"],dropped_invalid_site[each_key]["Site"]):
            bool_yes = False
            if "-" in get_each_pep:
              each_pep = (get_each_pep.split("-"))[0]
            else:
              each_pep = get_each_pep
            
            if each_key in merged_out_dict_2:
              if 'DroppedPeptide' not in merged_out_dict_2[each_key]:
                merged_out_dict_2[each_key].update({'DroppedPeptide':[each_pep]})
                merged_out_dict_2[each_key].update({'DroppedSite':[get_each_site]})
                bool_yes = True
              else:  
                if each_pep not in merged_out_dict_2[each_key]['DroppedPeptide']:               
                  merged_out_dict_2[each_key]['DroppedPeptide'].append(each_pep)
                  merged_out_dict_2[each_key]['DroppedSite'].append(get_each_site)
                  bool_yes = True
        
            else:
              merged_out_dict_2[each_key] = {}
              merged_out_dict_2[each_key].update({'DroppedPeptide':[each_pep]})
              merged_out_dict_2[each_key].update({'DroppedSite':[get_each_site]})
              bool_yes = True              
            
            if bool_yes:
              pep_len +=1            
              if 'Comment' not in merged_out_dict_2[each_key]:
                merged_out_dict_2[each_key].update({'Comment':"No site;"})
              else:
                merged_out_dict_2[each_key]['Comment'] += "No site;"            
            
        logging.debug("DISCARD WARNING - No sites available: " + str(pep_len) + " peptides")       
        
    else:
      initial_query_prots = each_protein_list + [x for x in dropped_invalid_fc_pval if x not in each_protein_list]
      initial_query_prot_count = len(initial_query_prots)
      uniq_initial_query_prot_count = len(set(initial_query_prots))
      logging.debug(f"Initial query: {uniq_initial_query_prot_count} unique protein IDs")
      
  initial_length = len(each_protein_list)
  to_return_unique_protids_length = len(set(each_protein_list))
  
  if type == "5" or type == "6":
    # need to recalculate each_protein_list since some inconsistent protein ids are deleted
    each_protein_list = list(site_info_dict.keys())

  if type == "1": 
    dropping_repeats = []
    initial_each_protein_list = each_protein_list    
    if repeat_prot_ids:
      unique_each_protein_list = [x for x in each_protein_list if x.lower() not in [name.lower() for name in repeat_prot_ids]]  
      dropping_repeats = [x for x in each_protein_list if x.lower() in [name.lower() for name in repeat_prot_ids]]        
      for x in list(prot_list):
        if x.lower() in [name.lower() for name in repeat_prot_ids]:
          del prot_list[x]
    
    else:
      unique_each_protein_list = each_protein_list
    if retain_prot_ids:
      each_protein_list = list(set(unique_each_protein_list))
    else:
      each_protein_list = unique_each_protein_list
    
    all_dropped = dropping_repeats + retain_prot_ids
    all_dropped = sorted(all_dropped)
    if cy_debug:
      if dropping_repeats:       
        logging.debug("DISCARD WARNING - Dropping proteins due to inconsistent fold changes or p-values between duplicates: " + str(len(list(set(dropping_repeats)))))
        list_of_duplicates = list(set(dropping_repeats))
        
  elif type == "3":
    if repeat_prot_ids:
      unique_each_protein_list = list(set(each_protein_list))            
      each_protein_list = unique_each_protein_list      

  elif type == "4":
    unique_each_protein_list = list(set(each_protein_list))          
    each_protein_list = unique_each_protein_list

  elif type == "2":
      unique_each_protein_list = list(set(each_protein_list))
      count_dropped = 0
      additional_dropped = []

      if repeat_prot_ids_2:  
        unique_each_protein_list = [x for x in unique_each_protein_list if x.lower() not in [name.lower() for name in repeat_prot_ids_2]]        
        for each_prot_2 in repeat_prot_ids_2.keys():
          for each_label_drop in repeat_prot_ids_2[each_prot_2]:
            count_dropped +=1              
            if each_prot_2 in prot_list:
              if each_label_drop in prot_list[each_prot_2]:
                del prot_list[each_prot_2][each_label_drop]
                additional_dropped.append(each_prot_2)
                count_dropped +=1
                if len(prot_list[each_prot_2]) == 0:
                  del prot_list[each_prot_2]                
      
      if cy_debug:
        repeat_prot_ids_key_append = []
        for list_each_dropped_protid in repeat_prot_ids_2:
          for value_of_key in repeat_prot_ids_2[list_each_dropped_protid]:
            repeat_prot_ids_key_append.append(list_each_dropped_protid)           
        list_of_duplicates = additional_dropped + repeat_prot_ids_key_append + retain_prot_ids
        list_of_duplicates = list(set(repeat_prot_ids_2))
        if list_of_duplicates:
          logging.debug("DISCARD WARNING - Dropping proteins due to inconsistent fold changes or p-values between duplicates: " + str(len(list(set(unique_each_protein_list))))) 
      each_protein_list = unique_each_protein_list
      max_FC_len = len(unique_labels)
      prot_list_rearrange = {}
      for each_protid in prot_list:
        for each_label in unique_labels:
          if each_label in prot_list[each_protid]:
            if each_protid not in prot_list_rearrange:
              prot_list_rearrange.update({each_protid:[[(prot_list[each_protid][each_label])[0]],[(prot_list[each_protid][each_label])[1]]]})
            else:
              prot_list_rearrange[each_protid][0].append((prot_list[each_protid][each_label])[0])
              prot_list_rearrange[each_protid][1].append((prot_list[each_protid][each_label])[1])
          else:
            if each_protid not in prot_list_rearrange:
              prot_list_rearrange.update({each_protid:[[0],[1.0]]})
            else:
              prot_list_rearrange[each_protid][0].append(0.0)
              prot_list_rearrange[each_protid][1].append(1.0)
      prot_list = {}
      prot_list = prot_list_rearrange
  
  elif type == "5" or type == "6":
    site_info_dict_rearrange = {}
    site_info_dict_rearrange1 = {}
    all_dropped_pep = {}
    all_dropped_warning = ""
    count_dropped = 0
    if type == "5":
      max_FC_len = 1
    else:
      max_FC_len = len(unique_labels)
    for each_protid in site_info_dict:
      for each_site in site_info_dict[each_protid]:
        non_unique = False
        pep_for_prot_site = ""
        keys = list(site_info_dict[each_protid][each_site].keys())
        if len(keys) > 1:
          all_mods = [v[2] for k,v in site_info_dict[each_protid][each_site].items()]
          unique_all_mods = [list(x) for x in set(tuple(x) for x in all_mods)]
          if len(unique_all_mods) > 1:
            non_unique = True
            for each_key_pep in site_info_dict[each_protid][each_site]:
              new_each_site = ""
              skip_val = False
              for each_modsite in site_info_dict[each_protid][each_site][each_key_pep][2]:
                new_each_site_num = re.sub("[^0-9]", "", each_modsite)
                match = re.match(r"([A-Za-z]+)([0-9]+)", each_site)
                if match:
                  if new_each_site:
                    new_each_site += "/"
                  new_each_site += match.group(1) + "{" + new_each_site_num + "}" + match.group(2)
              new_each_site_info = site_info_dict[each_protid][each_site][each_key_pep]
              i = 0              
              for each_fc,each_pval in zip(list(new_each_site_info[0]),list(new_each_site_info[1])): 
                try:
                  float(each_fc)
                  if math.isinf(float(each_fc)):
                    skip_val = True 
                except ValueError:
                  skip_val = True
                        
                try:
                  float(each_pval)
                  if math.isinf(float(each_pval)):
                    skip_val = True
                except ValueError:
                  skip_val = True
              
                if skip_val:
                  if type == "6":
                    del new_each_site_info[3][i]
                  del new_each_site_info[0][i]
                  del new_each_site_info[1][i]
                  i-=1
                  
                  if not each_site_info[0]:
                    if each_protid not in dropped_invalid_fc_pval:
                      dropped_invalid_fc_pval[each_protid] = {}
                      dropped_invalid_fc_pval[each_protid].update({"PeptideandLabel":[each_key_pep], "Site":[each_site]})
                    else:
                      dropped_invalid_fc_pval[each_protid]["PeptideandLabel"].append(each_key_pep)
                      dropped_invalid_fc_pval[each_protid]["Site"].append(each_site)
                i += 1

              if not (cy_fc_cutoff == 0.0 and cy_pval_cutoff == 1.0):
                cutoff_drop = inp_cutoff_ptms(cy_fc_cutoff,cy_pval_cutoff,new_each_site_info)              
                if cutoff_drop:
                  if each_protid not in dropped_cutoff_fc_pval:
                    dropped_cutoff_fc_pval[each_protid] = {}
                    dropped_cutoff_fc_pval[each_protid].update({"PeptideandLabel":[each_key_pep], "Site":[each_site]})
                  else:
                    dropped_cutoff_fc_pval[each_protid]["PeptideandLabel"].append(each_key_pep)
                    dropped_cutoff_fc_pval[each_protid]["Site"].append(each_site)
                  if type == "5":  
                    del new_each_site_info[0][0]
                    del new_each_site_info[1][0]
                  else:
                    new_each_site_info[0].clear()
                    new_each_site_info[1].clear()
                    
              if new_each_site_info[0]:              
                if each_protid not in site_info_dict_rearrange1:
                  site_info_dict_rearrange1[each_protid] = {}
                  site_info_dict_rearrange1[each_protid][new_each_site] = {}
                  site_info_dict_rearrange1[each_protid][new_each_site].update({each_key_pep:new_each_site_info})
                else:
                  if new_each_site not in site_info_dict_rearrange1[each_protid]:
                    site_info_dict_rearrange1[each_protid][new_each_site] = {}
                    site_info_dict_rearrange1[each_protid][new_each_site].update({each_key_pep:new_each_site_info})
                  else:
                    if each_key_pep not in site_info_dict_rearrange1[each_protid][new_each_site]:
                      site_info_dict_rearrange1[each_protid][new_each_site].update({each_key_pep:new_each_site_info}) 
                
          else:
            if exclude_ambi:
              each_site_info = []
              if each_protid not in merged_out_dict:
                merged_out_dict[each_protid] = {}
              for each_pep in site_info_dict[each_protid][each_site]:
                ambi_pep_count += 1
                count_dropped+=1
                if each_protid in all_dropped_pep:
                  all_dropped_pep[each_protid].append(each_pep)
                  if each_site not in ambigious_sites_ptms[each_protid]:
                    ambigious_sites_ptms[each_protid].append(each_site)
                    ambigious_sites[each_protid].append(each_site)                
                else:
                  all_dropped_pep.update({each_protid:[each_pep]})
                  ambigious_sites_ptms.update({each_protid:[each_site]})
                  ambigious_sites.update({each_protid:[each_site]})
                  
                if 'DroppedPeptide' not in merged_out_dict[each_protid]:
                  merged_out_dict[each_protid].update({'DroppedPeptide':[each_pep]})
                  merged_out_dict[each_protid].update({'DroppedSite':[each_site]})
                else:              
                  merged_out_dict[each_protid]['DroppedPeptide'].append(each_pep)
                  merged_out_dict[each_protid]['DroppedSite'].append(each_site)
                  
                if 'Comment' not in merged_out_dict[each_protid]:
                  merged_out_dict[each_protid].update({'Comment': "Ambiguous site;"})
                else:
                  merged_out_dict[each_protid]['Comment'] += "Ambiguous site;"
                  
            else:              
              each_site_info, dropped_pep, picked_pep = ptm_scoring(site_info_dict[each_protid][each_site], enzyme, include_list)
              pep_for_prot_site = picked_pep
              ambi_pep_count += len(dropped_pep) + 1
              if each_protid in dict_of_picked_pep:
                dict_of_picked_pep[each_protid].append(picked_pep)
              else:
                dict_of_picked_pep.update({each_protid:[picked_pep]})
            
              if each_protid not in merged_out_dict:
                merged_out_dict[each_protid] = {}

              if 'PickedPeptide' not in merged_out_dict[each_protid]:
                merged_out_dict[each_protid].update({'PickedPeptide':[picked_pep+"**"], 'PickedSite':[each_site+"**"]})
              else:
                merged_out_dict[each_protid]['PickedPeptide'].append(picked_pep+"**")
                merged_out_dict[each_protid]['PickedSite'].append(each_site+"**")
              
              for each_dropped_pep in dropped_pep:
                count_dropped+=1
                if each_protid in all_dropped_pep:
                  all_dropped_pep[each_protid].append(each_dropped_pep)
                  if each_site not in ambigious_sites_ptms[each_protid]:
                    ambigious_sites_ptms[each_protid].append(each_site)
                    ambigious_sites[each_protid].append(each_site)                
                else:
                  all_dropped_pep.update({each_protid:[each_dropped_pep]})
                  ambigious_sites_ptms.update({each_protid:[each_site]})
                  ambigious_sites.update({each_protid:[each_site]})
                
                if 'DroppedPeptide' not in merged_out_dict[each_protid]:
                  merged_out_dict[each_protid].update({'DroppedPeptide':[each_dropped_pep]})
                  merged_out_dict[each_protid].update({'DroppedSite':[each_site]})
                else:              
                  merged_out_dict[each_protid]['DroppedPeptide'].append(each_dropped_pep)
                  merged_out_dict[each_protid]['DroppedSite'].append(each_site)
                     
                if 'Comment' not in merged_out_dict[each_protid]:
                  merged_out_dict[each_protid].update({'Comment': "Ambiguous site;"})
                else:
                  merged_out_dict[each_protid]['Comment'] += "Ambiguous site;"
                       
        else:
          each_site_info = site_info_dict[each_protid][each_site][keys[0]]
          pep_for_prot_site = keys[0]
          if each_protid in dict_of_picked_pep:
            dict_of_picked_pep[each_protid].append(keys[0])
          else:
            dict_of_picked_pep.update({each_protid:[keys[0]]})
          
          if each_protid not in merged_out_dict:
            merged_out_dict[each_protid] = {}
              
          if 'PickedPeptide' not in merged_out_dict[each_protid]:
            merged_out_dict[each_protid].update({'PickedPeptide':[keys[0]], 'PickedSite':[each_site]})
          else:
            merged_out_dict[each_protid]['PickedPeptide'].append(keys[0])
            merged_out_dict[each_protid]['PickedSite'].append(each_site)
            
        if not non_unique and each_site_info:
          skip_val = False
          i = 0
          for each_fc,each_pval in zip(list(each_site_info[0]),list(each_site_info[1])):
            try:
              float(each_fc)
              if math.isinf(float(each_fc)):
                skip_val = True 
            except ValueError:
              skip_val = True
                        
            try:
              float(each_pval)
              if math.isinf(float(each_pval)):
                skip_val = True
            except ValueError:
              skip_val = True
           
            if skip_val:
              if type == "6":
                del each_site_info[3][i]
              del each_site_info[0][i]
              del each_site_info[1][i]
              i -=1
              if not each_site_info[0]: 
                if each_protid not in dropped_invalid_fc_pval:
                  dropped_invalid_fc_pval[each_protid] = {}
                  dropped_invalid_fc_pval[each_protid].update({"PeptideandLabel":[pep_for_prot_site], "Site":[each_site]})
                else:
                  dropped_invalid_fc_pval[each_protid]["PeptideandLabel"].append(pep_for_prot_site)
                  dropped_invalid_fc_pval[each_protid]["Site"].append(each_site)             
            i+=1

          if each_site_info[0] and not (cy_fc_cutoff == 0.0 and cy_pval_cutoff == 1.0):
            cutoff_drop = inp_cutoff_ptms(cy_fc_cutoff,cy_pval_cutoff,each_site_info)              
            if cutoff_drop:
              if each_protid not in dropped_cutoff_fc_pval:
                dropped_cutoff_fc_pval[each_protid] = {}
                dropped_cutoff_fc_pval[each_protid].update({"PeptideandLabel":[pep_for_prot_site], "Site":[each_site]})
              else:
                dropped_cutoff_fc_pval[each_protid]["PeptideandLabel"].append(pep_for_prot_site)
                dropped_cutoff_fc_pval[each_protid]["Site"].append(each_site)
              if type == "5":      
                del each_site_info[0][0]
                del each_site_info[1][0]
              else:
                each_site_info[0].clear()
                each_site_info[1].clear()              
              
          if each_site_info[0]:
            if each_protid not in site_info_dict_rearrange:
              site_info_dict_rearrange[each_protid] = {}
              site_info_dict_rearrange[each_protid].update({each_site:each_site_info})
            else:
              site_info_dict_rearrange[each_protid].update({each_site:each_site_info})
    
    if cy_debug:
      if dropped_invalid_site:
        i = 0
        for each_key in dropped_invalid_site:
          if dropped_invalid_site[each_key] == None:
            continue
          
          for get_each_pep,get_each_site in zip(dropped_invalid_site[each_key]["PeptideandLabel"],dropped_invalid_site[each_key]["Site"]):
            bool_yes = False
            if "-" in get_each_pep:
              each_pep = (get_each_pep.split("-"))[0]
            else:
              each_pep = get_each_pep
            
            if each_key in merged_out_dict:
              if 'DroppedPeptide' not in merged_out_dict[each_key]:
                merged_out_dict[each_key].update({'DroppedPeptide':[each_pep]})
                merged_out_dict[each_key].update({'DroppedSite':[get_each_site]})
                bool_yes = True
              else:  
                if each_pep not in merged_out_dict[each_key]['DroppedPeptide']:               
                  merged_out_dict[each_key]['DroppedPeptide'].append(each_pep)
                  merged_out_dict[each_key]['DroppedSite'].append(get_each_site)
                  bool_yes = True
        
            else:
              merged_out_dict[each_key] = {}
              merged_out_dict[each_key].update({'DroppedPeptide':[each_pep]})
              merged_out_dict[each_key].update({'DroppedSite':[get_each_site]})
              bool_yes = True              
            
            if bool_yes:
              if 'Comment' not in merged_out_dict[each_key]:
                merged_out_dict[each_key].update({'Comment':"No site;"})
              else:
                merged_out_dict[each_key]['Comment'] += "No site;"
                
    if site_info_dict_rearrange1:
      for each_prot_id in site_info_dict_rearrange1:
        for each_site in site_info_dict_rearrange1[each_prot_id]:
          if len(site_info_dict_rearrange1[each_prot_id][each_site].keys()) > 1:
            if exclude_ambi:
              again_each_site_info = []
              if each_prot_id not in merged_out_dict:
                merged_out_dict[each_prot_id] = {}
              for each_pep in site_info_dict_rearrange1[each_prot_id][each_site]:
                ambi_pep_count += 1
                count_dropped+=1
                if each_prot_id in all_dropped_pep:
                  all_dropped_pep[each_prot_id].append(each_pep)
                  if each_site not in ambigious_sites_ptms[each_prot_id]:
                    ambigious_sites_ptms[each_prot_id].append(each_site)
                    ambigious_sites[each_prot_id].append(each_site)                
                else:
                  all_dropped_pep.update({each_prot_id:[each_pep]})
                  ambigious_sites_ptms.update({each_prot_id:[each_site]})
                  ambigious_sites.update({each_prot_id:[each_site]})
                  
                if 'DroppedPeptide' not in merged_out_dict[each_prot_id]:
                  merged_out_dict[each_prot_id].update({'DroppedPeptide':[each_pep]})
                  merged_out_dict[each_prot_id].update({'DroppedSite':[each_site]})
                else:              
                  merged_out_dict[each_prot_id]['DroppedPeptide'].append(each_pep)
                  merged_out_dict[each_prot_id]['DroppedSite'].append(each_site)
                  
                if 'Comment' not in merged_out_dict[each_prot_id]:
                  merged_out_dict[each_prot_id].update({'Comment': "Ambiguous site;"})
                else:
                  merged_out_dict[each_prot_id]['Comment'] += "Ambiguous site;"
                
            else:
              again_each_site_info, dropped_pep, picked_pep = ptm_scoring(site_info_dict_rearrange1[each_prot_id][each_site], enzyme, include_list)
              pep_for_prot_site = picked_pep
              ambi_pep_count += len(dropped_pep) + 1
              if each_prot_id in dict_of_picked_pep:
                dict_of_picked_pep[each_prot_id].append(picked_pep)
              else:
                dict_of_picked_pep.update({each_prot_id:[picked_pep]})
            
              if each_prot_id not in merged_out_dict:
                merged_out_dict[each_prot_id] = {}

              if 'PickedPeptide' not in merged_out_dict[each_prot_id]:
                merged_out_dict[each_prot_id].update({'PickedPeptide':[picked_pep+"**"], 'PickedSite':[each_site+"**"]})
              else:
                merged_out_dict[each_prot_id]['PickedPeptide'].append(picked_pep+"**")
                merged_out_dict[each_prot_id]['PickedSite'].append(each_site+"**")
              
              for each_dropped_pep in dropped_pep:
                count_dropped+=1
                if each_prot_id in all_dropped_pep:
                  all_dropped_pep[each_prot_id].append(each_dropped_pep)
                  if each_site not in ambigious_sites_ptms[each_prot_id]:
                    ambigious_sites_ptms[each_prot_id].append(each_site)
                    ambigious_sites[each_prot_id].append(each_site)                
                else:
                  all_dropped_pep.update({each_prot_id:[each_dropped_pep]})
                  ambigious_sites_ptms.update({each_prot_id:[each_site]})
                  ambigious_sites.update({each_prot_id:[each_site]})
                
                if 'DroppedPeptide' not in merged_out_dict[each_prot_id]:
                  merged_out_dict[each_prot_id].update({'DroppedPeptide':[each_dropped_pep]})
                  merged_out_dict[each_prot_id].update({'DroppedSite':[each_site]})
                else:              
                  merged_out_dict[each_prot_id]['DroppedPeptide'].append(each_dropped_pep)
                  merged_out_dict[each_prot_id]['DroppedSite'].append(each_site)
                     
                if 'Comment' not in merged_out_dict[each_prot_id]:
                  merged_out_dict[each_prot_id].update({'Comment': "Ambiguous site;"})
                else:
                  merged_out_dict[each_prot_id]['Comment'] += "Ambiguous site;"
          else:
            
            each_key_pep = list(site_info_dict_rearrange1[each_prot_id][each_site].keys())[0]
            if each_prot_id in dict_of_picked_pep:
              dict_of_picked_pep[each_prot_id].append(each_key_pep)
            else:
              dict_of_picked_pep.update({each_prot_id:[each_key_pep]})
             
            if each_prot_id not in merged_out_dict:
              merged_out_dict[each_prot_id] = {}
                      
            if 'PickedPeptide' not in merged_out_dict[each_prot_id]:
              merged_out_dict[each_prot_id].update({'PickedPeptide':[each_key_pep], 'PickedSite':[each_site]})
            else:
              merged_out_dict[each_prot_id]['PickedPeptide'].append(each_key_pep)
              merged_out_dict[each_prot_id]['PickedSite'].append(each_site)          
          
            again_each_site_info = site_info_dict_rearrange1[each_prot_id][each_site][each_key_pep]
 
          if again_each_site_info:
            if each_prot_id not in site_info_dict_rearrange:
              site_info_dict_rearrange[each_prot_id] = {}
              site_info_dict_rearrange[each_prot_id].update({each_site:again_each_site_info})
            else:
              site_info_dict_rearrange[each_prot_id].update({each_site:again_each_site_info})
                        
    site_info_dict = site_info_dict_rearrange    
    each_protein_list = list(site_info_dict.keys())
    if cy_debug and all_dropped_pep:
      warn = []
      ambi_site_len = 0
      for each_prot in ambigious_sites_ptms:
        warn.append(each_prot + "(" + ",".join(ambigious_sites_ptms[each_prot]) + ")")
        ambi_site_len += len(ambigious_sites_ptms[each_prot])
      if exclude_ambi:
        logging.warning("AMBIGUITY WARNING - Dropping all ambigious sites: " + str(ambi_pep_count) + " peptides")
      else:
        logging.warning("AMBIGUITY WARNING - Dropping all but one representation of ambigious sites: " + str(ambi_pep_count-ambi_site_len) + " peptides")

    if cy_debug and len(duplicate_inc_ptm_proteins) > 0:
      duplicate_inc_full_dropped = set()
      for d in duplicate_inc_ptm_proteins:
        prot_pep = (d[0], d[2])
        if d[0] not in site_info_dict or d[1] not in site_info_dict[d[0]]:
          duplicate_inc_full_dropped.add(prot_pep)
      if len(duplicate_inc_full_dropped) > 0:
        msg_inc = "; ".join([f"(Protein: {d[0]}, Peptide: {d[1]})" for d in duplicate_inc_full_dropped])
        len_inc = len(duplicate_inc_full_dropped)
        logging.warning(f"DISCARD WARNING - Dropping peptides due to inconsistent fold changes or p-values between duplicates: {len_inc}")
        logging.warning(f"Dropped peptides: {msg_inc}")
 
  if type == "6":
    site_info_dict_rearrange = {}
    for each_protid_site in site_info_dict:
      for each_site in site_info_dict[each_protid_site]:
        for each_label in unique_labels:
          if each_label in site_info_dict[each_protid_site][each_site][3]:
            index_of = site_info_dict[each_protid_site][each_site][3].index(each_label)
            if each_protid_site not in site_info_dict_rearrange:
              site_info_dict_rearrange[each_protid_site] = {}              
              site_info_dict_rearrange[each_protid_site].update({each_site:[[(site_info_dict[each_protid_site][each_site])[0][index_of]],[(site_info_dict[each_protid_site][each_site])[1][index_of]]]})
            else:
              if each_site not in site_info_dict_rearrange[each_protid_site]:
                site_info_dict_rearrange[each_protid_site].update({each_site:[[(site_info_dict[each_protid_site][each_site])[0][index_of]],[(site_info_dict[each_protid_site][each_site])[1][index_of]]]})
              else: 
                site_info_dict_rearrange[each_protid_site][each_site][0].append((site_info_dict[each_protid_site][each_site])[0][index_of])
                site_info_dict_rearrange[each_protid_site][each_site][1].append((site_info_dict[each_protid_site][each_site])[1][index_of])
          else:
            if each_protid_site not in site_info_dict_rearrange:
              site_info_dict_rearrange[each_protid_site] = {}
              site_info_dict_rearrange[each_protid_site].update({each_site:[[0],[1.0]]})
            else:
              if each_site not in site_info_dict_rearrange[each_protid_site]:
                site_info_dict_rearrange[each_protid_site].update({each_site:[[0],[1.0]]})
              else:
                site_info_dict_rearrange[each_protid_site][each_site][0].append(0.0)
                site_info_dict_rearrange[each_protid_site][each_site][1].append(1.0)
    site_info_dict = {}
    site_info_dict = site_info_dict_rearrange 
    
  if cy_debug:
    if mapping_multiple_regions:
      warning = []
      for each_mult in mapping_multiple_regions:
        if "-" in each_mult:
          prot_only = each_mult.split("-")[0]
          pep_only = each_mult.split("-")[1]
          site_only = '|'.join(mapping_multiple_regions[each_mult])
        else:
          prot_only = each_mult
          pep_only = ""
          site_only = '|'.join(mapping_multiple_regions[each_mult])

        warning.append(each_mult + "(" + ','.join(str(x) for x in mapping_multiple_regions[each_mult]) + ")")
 
        if prot_only in merged_out_dict:
          if 'PickedPeptide' in merged_out_dict[prot_only]:          
            if pep_only and pep_only in merged_out_dict[prot_only]['PickedPeptide']:         
              indexOf = merged_out_dict[prot_only]['PickedPeptide'].index(pep_only)
              merged_out_dict[prot_only]['PickedPeptide'][indexOf] = merged_out_dict[prot_only]['PickedPeptide'][indexOf] + "**"
              merged_out_dict[prot_only]['PickedSite'][indexOf] = merged_out_dict[prot_only]['PickedSite'][indexOf] + "**"
                
          if 'DroppedPeptide' in  merged_out_dict[prot_only]:
            if pep_only and pep_only in merged_out_dict[prot_only]['DroppedPeptide']:         
              indexOf = merged_out_dict[prot_only]['DroppedPeptide'].index(pep_only)
              merged_out_dict[prot_only]['DroppedPeptide'][indexOf] = merged_out_dict[prot_only]['DroppedPeptide'][indexOf] + "**"
              
          if exclude_ambi:
            if 'DroppedPeptide' not in merged_out_dict[prot_only]:
              if pep_only:
                merged_out_dict[prot_only].update({'DroppedPeptide':[pep_only+"**"],'DroppedSite':[site_only],'Comment':"Map to multiple FASTA regions;"})
            else:
              if pep_only and pep_only not in merged_out_dict[prot_only]['DroppedPeptide']:
                merged_out_dict[prot_only]["DroppedPeptide"].append(pep_only+"**")
                merged_out_dict[prot_only]["DroppedSite"].append(site_only)
                merged_out_dict[prot_only]['Comment'] += "Map to multiple FASTA regions;"
                                        
        else:
          if exclude_ambi:
            if pep_only:
              merged_out_dict[prot_only] = {}
              merged_out_dict[prot_only].update({'DroppedPeptide':[pep_only+"**"],'DroppedSite':[site_only],'Comment':"Map to multiple FASTA regions;"})
                      
      if warning:
        if exclude_ambi:
          logging.warning("AMBIGUITY WARNING - Peptides mapped to multiple regions in FASTA: " + str(len(mapping_multiple_regions)) + " peptides")
          #logging.warning("Dropping queries: " + ','.join(warning))
        else:        
          logging.warning("AMBIGUITY WARNING - Peptides mapped to multiple regions in FASTA, first picked: " + ','.join(warning))
          
    if count_dup_pep and collect_dup_pep:
      logging.debug("DISCARD WARNING - Duplicate peptides across protein ids: " + str(len(list(set(count_dup_pep)))) + " peptides")
      for each_prot in collect_dup_pep:
        coll_pep = collect_dup_pep[each_prot]["Peptide"]
        coll_site = collect_dup_pep[each_prot]["Site"]
        if each_prot in merged_out_dict:
          if 'DroppedPeptide' not in merged_out_dict[each_prot]:
            if coll_pep:
              merged_out_dict[each_prot].update({'DroppedPeptide':coll_pep,'DroppedSite':coll_site,'Comment':"Duplicate peptide across protein IDs;"})
          else:
            if coll_pep and coll_pep not in merged_out_dict[each_prot]['DroppedPeptide']:
              merged_out_dict[each_prot]["DroppedPeptide"].extend(coll_pep)
              merged_out_dict[each_prot]["DroppedSite"].extend(coll_site)
              merged_out_dict[each_prot]['Comment'] += "Duplicate peptide across protein IDs;"
                  
    if dropped_invalid_fc_pval:
      site_len = 0
      pep_len = 0
      warning = []
      i = 0

      for each_key in dropped_invalid_fc_pval:
        warning_str = each_key
        if dropped_invalid_fc_pval[each_key] == None:
          if each_key not in merged_out_dict:
            merged_out_dict[each_key] = {} 
          if 'CommentGene' not in merged_out_dict[each_key]:
            merged_out_dict[each_key].update({'CommentGene':"Invalid FC/PVal;"})
          else:
            merged_out_dict[each_key]['CommentGene'] += "Invalid FC/PVal;"
          continue
        if "PeptideandLabel" in dropped_invalid_fc_pval[each_key]:
          if "Site" in dropped_invalid_fc_pval[each_key]:
            warning_str += "(" + ",".join(set(dropped_invalid_fc_pval[each_key]["Site"])) + ")"
            site_len += len(set(dropped_invalid_fc_pval[each_key]["Site"]))
            pep_len += len(dropped_invalid_fc_pval[each_key]["Site"])
           
        warning.append(warning_str)
        
        es = 0
        for each_pepandlabel in dropped_invalid_fc_pval[each_key]["PeptideandLabel"]:
          if "Site" in dropped_invalid_fc_pval[each_key]:
            e_s = dropped_invalid_fc_pval[each_key]["Site"][es]
          if "-" in dropped_invalid_fc_pval[each_key]["PeptideandLabel"]:
            each_pep = (each_pepandlabel.split("-"))[0]
          else:
            each_pep = each_pepandlabel
          
          if each_key in merged_out_dict:
            if 'PickedPeptide' in merged_out_dict[each_key]:
              if each_pep in merged_out_dict[each_key]['PickedPeptide']:
                indexOf = merged_out_dict[each_key]['PickedPeptide'].index(each_pep)
              elif each_pep+"**" in merged_out_dict[each_key]['PickedPeptide']:
                indexOf = merged_out_dict[each_key]['PickedPeptide'].index(each_pep+"**")
                each_pep = each_pep+"**"
              else:
                indexOf = -1
              if indexOf >=0:
                del merged_out_dict[each_key]['PickedPeptide'][indexOf]
                del merged_out_dict[each_key]['PickedSite'][indexOf] 
          else:
            merged_out_dict[each_key] = {}
            
          if 'DroppedPeptide' not in merged_out_dict[each_key]:
            merged_out_dict[each_key].update({'DroppedPeptide':[each_pep]})
            if "Site" in dropped_invalid_fc_pval[each_key]:
              merged_out_dict[each_key].update({'DroppedSite':[e_s]})
          else:              
            merged_out_dict[each_key]['DroppedPeptide'].append(each_pep)
            if "Site" in dropped_invalid_fc_pval[each_key]:
              merged_out_dict[each_key]['DroppedSite'].append(e_s)
        
          if 'Comment' not in merged_out_dict[each_key]:
            merged_out_dict[each_key].update({'Comment':"Invalid FC/PVal;"})
          else:
           merged_out_dict[each_key]['Comment'] += "Invalid FC/PVal;"
          es+=1
          
      if site_len > 0:
        logging.debug("DISCARD WARNING - Invalid Fold Change and P-Value terms: " + str(pep_len) + " peptides")
      else:
        excluded_from_count = list(prot_list.keys())
        logging.debug("DISCARD WARNING - Invalid Fold Change and P-Value terms: " + str(len(dropped_invalid_fc_pval)))
    
    if dropped_cutoff_fc_pval:
      site_len = 0
      pep_len = 0
      warning = []
      i = 0
      for each_key in dropped_cutoff_fc_pval:
        warning_str = each_key
        if dropped_cutoff_fc_pval[each_key] == None:
          continue
        if "PeptideandLabel" in dropped_cutoff_fc_pval[each_key]:
          if "Site" in dropped_cutoff_fc_pval[each_key]:
            warning_str += "(" + ",".join(set(dropped_cutoff_fc_pval[each_key]["Site"])) + ")"
            site_len += len(set(dropped_cutoff_fc_pval[each_key]["Site"]))
            pep_len += len(dropped_cutoff_fc_pval[each_key]["Site"])            
        warning.append(warning_str)
        
        es = 0
        for each_pepandlabel in dropped_cutoff_fc_pval[each_key]["PeptideandLabel"]:
          e_s = dropped_cutoff_fc_pval[each_key]["Site"][es]
          if "-" in each_pepandlabel:
            each_pep = (each_pepandlabel.split("-"))[0]
          else:
            each_pep = each_pepandlabel
          
          if each_key not in merged_out_dict:
            merged_out_dict[each_key] = {}
          if 'PickedPeptide' in merged_out_dict[each_key]:
            if each_pep in merged_out_dict[each_key]['PickedPeptide']:
              indexOf = merged_out_dict[each_key]['PickedPeptide'].index(each_pep)
            elif each_pep+"**" in merged_out_dict[each_key]['PickedPeptide']:
              indexOf = merged_out_dict[each_key]['PickedPeptide'].index(each_pep+"**")
              each_pep = each_pep+"**"
            else:
              indexOf = -1
            if indexOf >=0:
              del merged_out_dict[each_key]['PickedPeptide'][indexOf]
              del merged_out_dict[each_key]['PickedSite'][indexOf]          
                        
          if 'DroppedPeptide' not in merged_out_dict[each_key]:
            merged_out_dict[each_key].update({'DroppedPeptide':[each_pep]})
            merged_out_dict[each_key].update({'DroppedSite':[e_s]})
          else:              
            merged_out_dict[each_key]['DroppedPeptide'].append(each_pep)
            merged_out_dict[each_key]['DroppedSite'].append(e_s)
        
          if 'Comment' not in merged_out_dict[each_key]:
            merged_out_dict[each_key].update({'Comment':"FC/PVal cutoff not met;"})
          else:
           merged_out_dict[each_key]['Comment'] += "FC/PVal cutoff not met;"
          es+=1 
      if site_len > 0:
        logging.debug("DISCARD WARNING - Fold Change and P-Value cutoff not met: " + str(pep_len) + " peptides")  
                      
    if pep_not_in_fasta:
      warning = []
      count_warn = 0
      for each_mult in pep_not_in_fasta:
        if each_mult in each_protein_list:
          warning.append(each_mult + "(" + ','.join(pep_not_in_fasta[each_mult]) + ")")
          count_warn +=1 
      if warning:
        logging.debug("DISCARD WARNING - Peptides not found in FASTA: " + str(count_warn))
  
  if cy_debug:
    if type == "5" or type == "6":
      countpep = []
      for getprot,getsite in site_info_dict.items():
        countpep.extend(list(getsite.keys()))
      logging.debug("Remaining query: " + str(len(list(site_info_dict.keys()))) + " unique proteins IDs, " + str(len(list(set(countpep)))) + " unique peptides")
      each_protein_list = list(site_info_dict.keys())
    
  dup_prot_ids_to_return = []
  if type == "3":
    dup_prot_ids_to_return = repeat_prot_ids # nofc
  elif type == "4":
    dup_prot_ids_to_return = retain_prot_ids # category
  elif type == "1" or type == "2":
    dup_prot_ids_to_return = list_of_duplicates # multifc, singlefc
    each_protein_list = list(prot_list.keys())
  
  #if (len(unique_unimods)) > 1:
    #mult_mods_of_int = True 
  
  return(each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict, to_return_unique_protids_length, site_info_dict, ambigious_sites, unique_labels,dup_prot_ids_to_return, mult_mods_of_int)

def ptm_scoring(site_dict, enzyme, include_list):
  ''' PTM scoring algorithm for ambiguous sites: For proteins having different peptides for the same modification site, one representation of site is picked by choosing the peptide having no miscleavages or modifications other than the mod of interest '''
  enzyme_info = {'trypsin':{'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP', 'RP']}, 
                 'trypsin_p':{'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : []}, 
                 'lys_n':{'terminus' : 'N' , 'cleave' : ['K'], 'exceptions' : []}, 
                 'asp_n':{'terminus' : 'N' , 'cleave' : ['B','D'], 'exceptions' : []}, 
                 'arg_c':{'terminus' : 'C' , 'cleave' : ['R'], 'exceptions' : ['RP']}, 
                 'chymotrypsin':{'terminus' : 'C' , 'cleave' : ['F','Y','W','L'], 'exceptions' : ['FP','YP','WP','LP']}, 
                 'lys_c' : {'terminus' : 'C' , 'cleave' : ['K'], 'exceptions' : ['KP']}}
                 
  all_peptides = list(site_dict.keys())
  top_score = [0] * len(all_peptides)
  index = 0
  for each_peptide in all_peptides: 
    # Other Mods
    total_mods = 0
    om = 0
    all_mods_dict = find_mod(each_peptide)
    combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
    each_peptide_without_mods = re.sub(combined_pat,'',each_peptide)
    for k in include_list:
      for key,value in all_mods_dict.items():         
        key = re.sub('\+','',key)                               
        k = re.sub('\+','',k)
        if "c[" in key.lower() or "c{" in key.lower() or "c(" in key.lower():
          total_mods += 1
          continue
        key_match = re.search(combined_pat,key)
        k_match = re.search(combined_pat,k)
        if not (key_match and k_match):
          key = re.sub(combined_pat, '', key)
          k = re.sub(combined_pat,'',k)           
        if k.lower() in key.lower():
          total_mods += 1   
    om = len(all_mods_dict) - total_mods

    # Miscleavage  
    miscleave = 0
    exception = 0
    for each_cleave in enzyme_info[enzyme.lower()]['cleave']:
      miscleave += each_peptide_without_mods.count(each_cleave)
    if miscleave > 0:
      miscleave -=1
      
    # Exclude exceptions
    for each_except in enzyme_info[enzyme.lower()]['exceptions']:
      exception += each_peptide_without_mods.count(each_except)
    
    top_score[index] += miscleave - exception + om
    index += 1  
  # Pick lowest score
  min_score = min(top_score)
  indices = [i for i, x in enumerate(top_score) if x == min_score]
  largest_avg_FC = 0
  index = 0
  if len(indices) > 1:
    for each_index in indices:
      each_pick_pep = all_peptides[each_index]
      n = 0
      avg_FC = 0     
      for each_site_pos in site_dict[each_pick_pep][0]:
        try:
          float(each_site_pos)
          avg_FC += abs(each_site_pos)
          n+= 1
        except:
          continue
      if n!= 0:    
        avg_FC = avg_FC/n
        if abs(avg_FC) >= largest_avg_FC:     
          largest_avg_FC = abs(each_site_pos)
          index = each_index
  else:
    index = indices[0]
  pick_pep = all_peptides[index]
  dropped_pep = [x for x in all_peptides if not x == pick_pep]
  
  return(site_dict[pick_pep],dropped_pep,pick_pep)  
    
def inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict):
  ''' Filter input based on user defined fold change or p-value cutoffs '''
  queries_dropped = []
  for each_prot in list(prot_list):
    delete_each_prot = False
    for each_fc_val,each_pval in zip(prot_list[each_prot][0],prot_list[each_prot][1]):
      if not (abs(float(each_fc_val)) >= abs(float(cy_fc_cutoff)) and float(each_pval) <= float(cy_pval_cutoff)):
        delete_each_prot = True
      else:
        delete_each_prot = False
        break
        
    if delete_each_prot:
      del prot_list[each_prot]
      unique_each_protein_list.remove(each_prot)
      queries_dropped.append(each_prot)
      if each_prot not in merged_out_dict:
        merged_out_dict[each_prot] = {}
      if 'Primary' not in merged_out_dict[each_prot]:
        merged_out_dict[each_prot].update({'Primary':'' , 'String':'', 'Genemania':'', 'ClueGO':''})
      if 'CommentGene' in merged_out_dict[each_prot]:
        merged_out_dict[each_prot]['CommentGene'] += 'FC/Pval cutoff not met;'
      else:
        merged_out_dict[each_prot].update({'CommentGene':'FC/Pval cutoff not met;'})
      
  if cy_debug:
    logging.debug("DISCARD WARNING - Fold Change and P-Value cutoff not met: " + str(len(queries_dropped)))
    
  return(unique_each_protein_list, prot_list, merged_out_dict)

def inp_cutoff_ptms(cy_fc_cutoff,cy_pval_cutoff,each_site_info):
  ''' Filter input based on user defined fold change or p-value cutoffs for PTM based analysis '''
  delete_site = False
  for each_fc_val,each_pval in zip(each_site_info[0],each_site_info[1]):
    try:         
      if not (abs(float(each_fc_val)) >= abs(float(cy_fc_cutoff)) and float(each_pval) <= float(cy_pval_cutoff)):
        delete_site = True
      else:
        delete_site = False
        break
    except:
      break
  return(delete_site)
  
def uniprot_api_call(each_protein_list, prot_list, type, cy_debug, logging, merged_out_dict, species, cy_session, cy_out, cy_cluego_out, cy_cluego_in, path_to_new_dir, logging_file, site_info_dict, exclude_ambi, cy_settings_file):
  ''' 
  API call to Uniprot to correlate protein ids to their corresponding gene names. Ambiguity cases dealt with:
    1. Query protein ID mapping to multiple Uniprot IDs
    2. Query protein ID mapping to multiple primary genes
    3. Isoforms in the query which will map to same protein ID
  '''
  uniprot_query = {}
  each_primgene_list = []
  each_protid_list = []
  no_uniprot_val = []
  no_primgene_val = []
  prot_with_mult_primgene = []
  query_term = ' '.join(each_protein_list)
  url = 'https://www.uniprot.org/uploadlists/'
  ambigious_gene = []
  all_isoforms = []
  count_iso = 0
  duplicate_canonical = []
  retain_dup_can = []
  params = {
  'from':'ACC+ID',
  'to':'ACC',
  'format':'tab',
  'query':query_term,
  'columns': 'id,genes(PREFERRED),genes(ALTERNATIVE),organism,reviewed,last-modified'
  }
  try:
    response = requests.post(url, data=params)
    response.raise_for_status()
    decode = response.text
    list1=decode.split('\n')
    list1 = list1[1:]
    req_ids = []
    for l in list1:
      if l:
        req_ids.extend(l.split("\t")[6].split(","))
    id_counts = Counter(req_ids)
    duplicate_mapping_ids = set([x for x in id_counts if id_counts[x] > 1])
    seen_duplicate_mapping_ids = {}
  except:
    eprint("Error: Uniprot not responding. Please try again later")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  isoform_warning = ""
  for each_list1 in list1:
    remaining_isoforms = []
    is_mult_prim_gene_bool = False
    is_isoform_gene_bool = False
    if each_list1 != '':
      uniprot_list = each_list1.split('\t')
      uniprot_protid = uniprot_list[0]
      prot = uniprot_list[6]
      split_prot_list = prot.split(',')
      comment_merged = ""     
      if len(split_prot_list) > 1:
        is_isoform_gene_bool = True
        if exclude_ambi:
          if len(set([x.split("-")[0] for x in split_prot_list])) > 1: # more than one canonical ID after accounting for isoforms
            duplicate_canonical.append((uniprot_protid, split_prot_list))
            for spl in split_prot_list:
              input_failure_comment(uniprot_query, merged_out_dict, spl, "Multiple input IDs map to a single Uniprot ID")            
          else:
            remaining_isoforms = split_prot_list
            for spl in split_prot_list:
              input_failure_comment(uniprot_query, merged_out_dict, spl, "Isoform")
            all_isoforms.extend(remaining_isoforms)            
          each_prot = ""
        else:
          if len(set([x.split("-")[0] for x in split_prot_list])) > 1: # more than one canonical ID after accounting for isoforms
            get_list_dup_can = []
            is_retained_dup_can = ""
            for each_dup_can in split_prot_list:
              if each_dup_can == uniprot_protid:
                is_retained_dup_can = each_dup_can
                retain_dup_can.append(each_dup_can)
              else:
                get_list_dup_can.append(each_dup_can)
                input_failure_comment(uniprot_query, merged_out_dict, each_dup_can, "Multiple input IDs map to a single Uniprot ID")
            each_prot = is_retained_dup_can
            duplicate_canonical.append((uniprot_protid, get_list_dup_can))
          else:
            remaining_isoforms = []
            get_each_picked = ""
            for each_iso in split_prot_list:
              if "-" not in each_iso:
                get_each_picked = each_iso
              else:
                remaining_isoforms.append(each_iso)

            if not get_each_picked:
              remaining_isoforms = split_prot_list[1:]
              get_each_picked = split_prot_list[0]
            each_prot = get_each_picked
            count_iso += len(remaining_isoforms)
            for spl in remaining_isoforms:
              input_failure_comment(uniprot_query, merged_out_dict, spl, "Isoform")
            isoform_warning += "(" + ','.join([each_prot] + remaining_isoforms) + "),"
            all_isoforms.extend(remaining_isoforms)
      else:
        each_prot = split_prot_list[0]      

      if each_prot in duplicate_mapping_ids:
        if each_prot in seen_duplicate_mapping_ids:
          seen_duplicate_mapping_ids[each_prot].append({"protein": uniprot_protid, "gene": uniprot_list[1]})
          each_prot = ""
        else:
          seen_duplicate_mapping_ids[each_prot] = [{"protein": uniprot_protid, "gene": uniprot_list[1]}]
          input_failure_comment(uniprot_query, merged_out_dict, each_prot, "Input ID maps to multiple Uniprot IDs")
          each_prot = ""

      if each_prot == "":
        continue
      
      primary_gene = uniprot_list[1]
      if each_prot and each_prot not in merged_out_dict:
        merged_out_dict[each_prot] = {}
      
      # Pick first gene in case of multiple primary genes for single uniprot ids    
      if ";" in primary_gene:
        prot_with_mult_primgene.append(each_prot + "(" + primary_gene + ") ")
        is_mult_prim_gene_bool = True
        if not exclude_ambi:
          primary_gene = primary_gene.split(";")[0]
         
        if primary_gene and primary_gene not in ambigious_gene:
          ambigious_gene.append(primary_gene)
          
      if remaining_isoforms:
        if primary_gene not in ambigious_gene:
          ambigious_gene.append(primary_gene) 
          
      synonym_gene = uniprot_list[2]
      synonym_gene = synonym_gene.split(" ")
      organism_name = uniprot_list[3]
      uniprot_query[each_prot] = {}
      uniprot_query[each_prot].update({'Uniprot':uniprot_protid,'Primary':primary_gene,'Synonym':synonym_gene,'Organism':organism_name,'Reviewed':uniprot_list[4],'Date_modified':uniprot_list[5]})
      # Do not add empty primary gene to list
      if primary_gene and ";" not in primary_gene:
        each_primgene_list.append(primary_gene)
      elif not primary_gene:
        no_primgene_val.append(each_prot)
        comment_merged = "Primary gene unavailable;"
      elif ";" in primary_gene and exclude_ambi:
        comment_merged = "Multiple primary genes;"
      if uniprot_protid:
        each_protid_list.append(uniprot_protid)
      
      if not exclude_ambi and is_isoform_gene_bool and primary_gene:
        store_prim_gene = primary_gene + "**"       
      elif not exclude_ambi and is_mult_prim_gene_bool:
        store_prim_gene = primary_gene + "**"
      else:
        store_prim_gene = primary_gene
                 
      merged_out_dict[each_prot].update({'Primary':store_prim_gene, 'String':'', 'Genemania':'', 'ClueGO':''})
      if 'CommentGene' in merged_out_dict[each_prot]:
        merged_out_dict[each_prot]['CommentGene'] += comment_merged
      else:
        merged_out_dict[each_prot].update({'CommentGene':comment_merged})
      
  should_exit = False
  for each_prot_in_input in each_protein_list:
    if each_prot_in_input not in uniprot_query:
      no_uniprot_val.append(each_prot_in_input)
      uniprot_query[each_prot_in_input] = {}
      uniprot_query[each_prot_in_input].update({'Uniprot':"NA",'Primary':"NA",'Synonym':"NA",'Organism':"NA",'Reviewed':"NA",'Date_modified':"NA"})
      merged_out_dict[each_prot_in_input] = {}
      comment_merged = "Uniprot query not mapped;" 
      merged_out_dict[each_prot_in_input].update({'Primary':'', 'CommentGene':comment_merged, 'String':'', 'Genemania':'', 'ClueGO':''})
        
  date_modified = [re.sub("-","",dict['Date_modified']) for dict in uniprot_query.values() if dict['Date_modified'] != "NA"]
  
  if cy_cluego_in:
    try:
      settings_file = os.path.join(os.path.dirname(cy_cluego_in), "timestamp.json")
      if not os.path.exists(settings_file):
        eprint("Error: Settings file is missing for this run. Please restart your run")
        should_exit = True
      else:
        with open(settings_file) as f:
          settings = json.load(f)
        match = re.search(r"^(\d{4})-(\d{2})-(\d{2})T", settings["timestamp"])
        run_timestamp = float(match.group(1) + match.group(2) + match.group(3))
        for each_date_modified in date_modified:
          if float(each_date_modified) > float(run_timestamp):
            should_exit = True
            break
        if should_exit:
          eprint("Error: Uniprot has updated since the original run. Please restart your run")
    except json.JSONDecodeError:
      eprint("Error: Could not parse settings file. Please restart your run")
      should_exit = True
    except:
      eprint("Error: Could not parse settings file. Please restart your run")
      should_exit = True
     
  if should_exit:
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1) 
    
  organisms = [dict['Organism'] for dict in uniprot_query.values() if dict['Organism'] != "NA" and dict['Organism'] != ""]
  unique_organisms = list(set(organisms))
  if len(unique_organisms) > 1:
    count_each_organism = 0
    organism_of_count = ""
    list_of_prots = []
    for count_org in unique_organisms:
      val = [dict['Uniprot'] for dict in uniprot_query.values() if dict['Organism'] == count_org]
      count_val = len([dict['Uniprot'] for dict in uniprot_query.values() if dict['Organism'] == count_org])
      if count_each_organism == 0:
        count_each_organism = count_val
        organism_of_count = count_org
        list_of_prots = val
      else:
        if count_val < count_each_organism:
          count_each_organism = count_val
          organism_of_count = count_org 
          list_of_prots = val
    eprint("Error: Protein list is of more than 1 organism: " + ','.join(unique_organisms) + ". Please check Uniprot IDs " + ",".join(list_of_prots) + " of species " + organism_of_count)
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
    
  if unique_organisms:
    if not species.lower() in unique_organisms[0].lower():
      eprint("Error: Species mismatch. Species parameter provided is " + species + " and species of protein list is " + unique_organisms[0])
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)
   
  if cy_debug:
    if all_isoforms:
      if exclude_ambi:
        logging.warning("AMBIGUITY WARNING - Isoforms in query: " + str(len(all_isoforms)))
        logging.warning("Dropping queries: " + ','.join(all_isoforms))
      else:
        logging.warning("AMBIGUITY WARNING - Isoforms in query: " + str(count_iso))
        logging.warning("Dropping all but first query: " + isoform_warning.rstrip(","))

    duplicate_canonical_set = set()
    for dc in duplicate_canonical:
      for dc1 in dc[1]:
        duplicate_canonical_set.add(dc1)
    if duplicate_canonical:
      duplicate_canonical_str = ", ".join([ x[0] + " (" + ", ".join(x[1]) + ")" for x in duplicate_canonical])
      logging.warning(f"AMBIGUITY WARNING - Multiple query IDs mapping to single ID: {len(duplicate_canonical_set)}")
      logging.warning(f"Dropping queries: {duplicate_canonical_str}")
      if retain_dup_can:
        logging.warning("Retaining active queries: " + ",".join(retain_dup_can))

    if len(seen_duplicate_mapping_ids) > 0:
      sdmi_list = []
      for sdmi in seen_duplicate_mapping_ids:
        if sdmi in duplicate_canonical_set:
          continue # skip ids that were dropped from duplicate_canonical
        sdmi_val = seen_duplicate_mapping_ids[sdmi]
        sdmi_list.append(sdmi + "(" + ", ".join([x["protein"] for x in sdmi_val]) + ")")
      duplicate_mapping_str = ", ".join(sdmi_list)
      if len(sdmi_list) > 0:
        logging.warning(f"AMBIGUITY WARNING - Single query IDs mapping to multiple IDs: {len(sdmi_list)}")
        logging.warning(f"Dropping queries: {duplicate_mapping_str}")
        
    if prot_with_mult_primgene:
      if exclude_ambi:
        logging.warning("AMBIGUITY WARNING - Uniprot Multiple primary genes in query: " + str(len(prot_with_mult_primgene)))
        #logging.warning("Dropping queries: " + ','.join(prot_with_mult_primgene))
      else:
        logging.warning("AMBIGUITY WARNING - Uniprot Multiple primary genes, first picked: " + ','.join(prot_with_mult_primgene))
  
    if no_uniprot_val:
      logging.debug("DISCARD WARNING - Uniprot query not mapped: " + str(len(no_uniprot_val)))
      #logging.warning("Dropping queries: " + ','.join(no_uniprot_val))

    if no_primgene_val:
      logging.debug("DISCARD WARNING - Uniprot Primary gene unavailable: " + str(len(no_primgene_val)))
      #logging.warning("Dropping queries: " + ','.join(no_primgene_val))
  
  if type == "1" or type == "2":
    for each_in_list in prot_list:
      if each_in_list in uniprot_query:
        uniprot_query[each_in_list].update({"FC":prot_list[each_in_list][0]})
        uniprot_query[each_in_list].update({"PVal":prot_list[each_in_list][1]})
  
  elif type == "4":
    for each_in_list in prot_list:
      if each_in_list in uniprot_query:
        uniprot_query[each_in_list].update({"Category":prot_list[each_in_list]})      

  return(uniprot_query,each_primgene_list,merged_out_dict,ambigious_gene)

def get_dbsnp_classification(uniprot_query, primgene_list, prot_list, merged_out_dict):
  ''' Not used. Obtain SNP data for the list of protein ids through an API call to Uniprot '''
  variants = {}
  url = 'https://www.uniprot.org/uploadlists/'
  query_term = "gene_exact:" + ' or gene_exact:'.join(primgene_list)
  list_of_genes = " ".join(primgene_list)
  len=0
  params = {
  'from':'GENENAME',
  'to':'ACC',
  'format':'tab',
  'query':list_of_genes,
  'columns': 'genes(PREFERRED),organism,feature(NATURAL VARIANT)'
  }
  data = urllib.parse.urlencode(params).encode("utf-8")
  request = urllib2.Request(url, data)
  response = urllib2.urlopen(request)
  page = response.read()
  decode = page.decode("utf-8")
  list1=decode.split('\n')
  list1 = list1[1:]
  for each_list1 in list1:
    if each_list1 != '':
      uniprot_list = each_list1.split('\t')
      primary_gene = uniprot_list[0].split(";")
      organism = uniprot_list[1]
      if "human" in organism.lower():
        natural_variant = uniprot_list[2]       
        for uniprot_id in uniprot_query:  
          if uniprot_query[uniprot_id]['Primary'].lower() in [x.lower() for x in primary_gene]:
            # Get all variants 
            get_all_variants = natural_variant.split(" VARIANT")
            for each_variant in get_all_variants:
              m = re.search(r"\d+\s(\d+)\s([A-Z]{1}) -> ([A-Z]{1})",each_variant)
              if m:
                dbsnp = ""
                disease = ""
                ftid = ""
                position = m.group(1)
                AA_change_from = m.group(2)
                AA_change_to = m.group(3)
                if "unknown pathological significance" in each_variant:
                  classification = "Unclassified"
                else:
                  is_disease = re.search(r"\(([A-Z0-9,\sa-z]+)", each_variant)         
                  if (is_disease and not "dbSNP" in is_disease.group(1)):
                    is_unclassified = re.search(r"[a-z]+", ((is_disease.group(1)).replace('and','')).replace('in',''))
                    if is_unclassified:
                      classification = "Unclassified"
                    else:
                      classification = "Disease"
                      disease = ((is_disease.group(1)).replace(',',';')).replace('in ','')
                  else:
                    classification = "Polymorphism"
                dbsnp_match = re.search(r"dbSNP:(rs[0-9]+)",each_variant)
                if dbsnp_match:
                  dbsnp = dbsnp_match.group(1)
                ftid_match = re.search(r"/FTId=(VAR_[0-9]+)",each_variant)
                if ftid_match:
                  ftid = ftid_match.group(1)
    
                if uniprot_id in variants:
                  variants[uniprot_id]['Position'].append(position)
                  variants[uniprot_id]['FTID'].append(ftid)
                  variants[uniprot_id]['AA_change_from'].append(AA_change_from)
                  variants[uniprot_id]['AA_change_to'].append(AA_change_to)
                  variants[uniprot_id]['Classification'].append(classification)
                  variants[uniprot_id]['Disease'].append(disease)
                  variants[uniprot_id]['dbSNP'].append(dbsnp)
                else:
                  variants[uniprot_id] = {}
                  variants[uniprot_id].update({'Position':[position], 'FTID':[ftid], 'AA_change_from':[AA_change_from], 'AA_change_to':[AA_change_to], 'Classification':[classification], 'Disease':[disease], 'dbSNP':[dbsnp]})
   
  all_prot_site_snps = {}
  for each_prot in variants:
    if each_prot in prot_list:
      for each_pos in prot_list[each_prot]:
        combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}','[A-Za-z]'))
        combined_pat2 = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}','[0-9]'))
        each_pos_sub = re.sub(combined_pat, '', each_pos)
        each_aa_sub = re.sub(combined_pat2, '', each_pos)
        if str(each_pos_sub) in variants[each_prot]['Position']:
          indexOf = variants[each_prot]['Position'].index(str(each_pos_sub))
          if (each_aa_sub.lower() == variants[each_prot]['AA_change_from'][indexOf].lower()):
            if each_prot in all_prot_site_snps:
              all_prot_site_snps[each_prot].append(each_pos)
            else:
              all_prot_site_snps.update({each_prot:[each_pos]})
            if each_prot in merged_out_dict:
              if variants[each_prot]['Disease'][indexOf]:
                is_disease = variants[each_prot]['Disease'][indexOf]
              else:
                is_disease = "NA"              
              if "AA Change" in merged_out_dict[each_prot]:
                merged_out_dict[each_prot]["AA Change"].append(variants[each_prot]['Position'][indexOf] + "(" + variants[each_prot]['AA_change_from'][indexOf] + "->" + variants[each_prot]['AA_change_to'][indexOf] + ")")
                merged_out_dict[each_prot]["Classification"].append(variants[each_prot]['Classification'][indexOf])
                merged_out_dict[each_prot]["Disease"].append(is_disease)
              else:
                merged_out_dict[each_prot].update({"AA Change":[variants[each_prot]['Position'][indexOf] + "(" + variants[each_prot]['AA_change_from'][indexOf] + "->" + variants[each_prot]['AA_change_to'][indexOf] + ")"], "Classification":[variants[each_prot]['Classification'][indexOf]], "Disease":[is_disease]})
                        
  return(all_prot_site_snps, variants)
  
def get_query_from_list(uniprot_query, list):
  ''' Retrieve protein id for a gene name. Used as key-value mapping in different parts of the code ''' 
  dropped_list = {}
  for key in uniprot_query:
    for each_dupe_gene in list:
      if uniprot_query[key]['Primary'].lower() == each_dupe_gene.lower():
        dropped_list.update({key:each_dupe_gene})
        break
  return(dropped_list)

def check_dup_preferred_gene(mapping, interaction, unique_vertex, search, uniprot_query, cy_debug, logging, merged_out_dict):
  ''' Function to detect and discard all duplicates in preferred gene name for interaction databases String and GeneMANIA '''
  dupe_preferred_list = []
  retain_pf_nodes = []
  dropped_primgene_nodes = []
  dropped_pf_nodes = []
  unique_each_preferred_list = []
  dropped_interaction_list = []
  get_dupe_query_primgenes = []
  get_dupe_query_proids = []
  mapping_values = list(mapping.values())
  dupe_in_interaction_only = []
  discarded_interactions = []
  unique_nodes_of_dropped_interactions = []
  warning = []
  
  if len(mapping_values) != len(set(mapping_values)):
    dupe_preferred_list = remove_list_duplicates(mapping_values)
       
    for key, val in mapping.items():
      for each_dupe in dupe_preferred_list:
        if val.lower() == each_dupe.lower() and key.lower() == val.lower():
          retain_pf_nodes.append(each_dupe)
        elif val.lower() == each_dupe.lower() and key.lower() != val.lower():
          dropped_primgene_nodes.append(key)
          if each_dupe not in dropped_pf_nodes:
            dropped_pf_nodes.append(each_dupe)
    
    actual_dropped_nodes = [x.lower() for x in dropped_pf_nodes if x.lower() not in [name.lower() for name in retain_pf_nodes]]
    lower_dupe_preferred_list = actual_dropped_nodes
    #Drop interactions with dropped dupe preferred list nodes
    for each_interaction in interaction:
      node1 = each_interaction.split(" ")[0]
      node2 = each_interaction.split(" ")[1]
      if node1.lower() not in lower_dupe_preferred_list and node2.lower() not in lower_dupe_preferred_list:
        dropped_interaction_list.append(each_interaction)
        if node1.lower() not in [name.lower() for name in unique_each_preferred_list]:
          unique_each_preferred_list.append(node1)
        if node2.lower() not in [name.lower() for name in unique_each_preferred_list]:
          unique_each_preferred_list.append(node2)
        
      else:
        if node1.lower() in lower_dupe_preferred_list and node1.lower() not in [name.lower() for name in dupe_in_interaction_only]:
          dupe_in_interaction_only.append(node1)
        if node2.lower() in lower_dupe_preferred_list and node2.lower() not in [name.lower() for name in dupe_in_interaction_only]:
          dupe_in_interaction_only.append(node2)
        discarded_interactions.append(each_interaction)          
        if node1.lower() not in [name.lower() for name in unique_nodes_of_dropped_interactions]:
          unique_nodes_of_dropped_interactions.append(node1)
        if node2.lower() not in [name.lower() for name in unique_nodes_of_dropped_interactions]:
          unique_nodes_of_dropped_interactions.append(node2)
        
    if dropped_primgene_nodes: 
      get_dupe_query_proids = get_query_from_list(uniprot_query,dropped_primgene_nodes)
    
    get_dupe_query_primgenes_2 = {}
    for each_drop in unique_nodes_of_dropped_interactions:
      if each_drop not in unique_each_preferred_list and each_drop not in dupe_in_interaction_only:
        for k,v in mapping.items():
          if v.lower() == each_drop.lower():
            get_dupe_query_primgenes_2.update({k:v})
    get_dupe_query_proids_2 = get_query_from_list(uniprot_query,list(get_dupe_query_primgenes_2.keys()))
        
    if len(get_dupe_query_proids) != 0:     
      for each_get_dupe_primgene in dropped_primgene_nodes:
        each_get_dupe_ids = [key for key, val in get_dupe_query_proids.items() if val.lower() == each_get_dupe_primgene.lower()][0]
        warning.append(each_get_dupe_ids + "(" + get_dupe_query_proids[each_get_dupe_ids] + ") ") 
    
    if len(get_dupe_query_proids_2) != 0:    
      for each_get_dupe_ids in get_dupe_query_proids_2:
        warning.append(each_get_dupe_ids + "(" + get_dupe_query_proids_2[each_get_dupe_ids] + ") ")  
      
  else:
    return(interaction, unique_vertex, dupe_preferred_list, merged_out_dict, warning)
    
  return(dropped_interaction_list, unique_each_preferred_list, dupe_preferred_list, merged_out_dict, warning)
  
def overwrite_values(overwrite_preferred, interaction, unique_vertex, mapping, merged_out_dict):
  ''' If uniprot gene name and preferred gene names differ only by '-' symbol, categorize gene as primary '''  
  for key,val in mapping.items():
    if val.lower() in overwrite_preferred:
      mapping[key] = overwrite_preferred[val.lower()] 
      
  for n, i in enumerate(unique_vertex):
    if i.lower() in overwrite_preferred:
      unique_vertex[n] = overwrite_preferred[i.lower()]

  for n, each_interaction in enumerate(interaction):
    node1 = each_interaction.split(" ")[0]
    node2 = each_interaction.split(" ")[1]
    if node1.lower() in overwrite_preferred or node2.lower() in overwrite_preferred:
      new_node1 = overwrite_preferred[node1.lower()] if node1.lower() in overwrite_preferred else node1
      new_node2 = overwrite_preferred[node2.lower()] if node2.lower() in overwrite_preferred else node2
      interaction[n] = new_node1 + " " + new_node2
  return(interaction, unique_vertex, mapping, merged_out_dict)

def drop_ei_if_query(interaction,unique_vertex,genes_before_initial_drop,each_inp_list):
  ''' Drop all External Interactors if it is in the original input list '''
  dropped_nodes = []
  dropped_interactions = []
  remaining_nodes = []
  remaining_interactions = []
  for each_node in unique_vertex:
    if each_node.lower() in [name1.lower() for name1 in genes_before_initial_drop] and each_node.lower() not in [name2.lower() for name2 in each_inp_list]:
      dropped_nodes.append(each_node)
      for each_interaction in interaction:
        node1 = each_interaction.split(" ")[0]
        node2 = each_interaction.split(" ")[1]       
        if node1.lower() == each_node.lower() or node2.lower() == each_node.lower():
          if each_interaction not in dropped_interactions:
            dropped_interactions.append(each_interaction)
    else:
      remaining_nodes.append(each_node)

  remaining_interactions = [x for x in interaction if x.lower() not in [name.lower() for name in dropped_interactions]]    
 
  return(remaining_interactions,remaining_nodes)
  
def create_string_cytoscape(uniprot_query,each_inp_list, species, limit, score, cy_debug, logging, merged_out_dict, genes_before_initial_drop, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file):
  ''' 
  Obtain gene-gene interactions from String database:
    1. Get mapping table and categorize genes as primary, secondary, external synonym or external interactor.
    2. Duplicates in preferred gene names are also dropped 
    3. Unique list of interactions is returned.
  '''    
  string_db_out = {}
  string_mapping = {}
  string_interaction = {}
  string_unique_vertex = {}
  string_list_input = "\n".join(each_inp_list)
  
  data = {
  'identifiers': string_list_input,
  'species': species,
  'echo_query': 1,
  "limit" : 1
  }
  try:
    response = requests.post('https://string-db.org/api/tsv-no-header/get_string_ids', data=data)
    response.raise_for_status()
    out = pd.read_table(StringIO(response.text), header=None)
    query = list(out[0].values.flatten())
    preferred = list(out[5].values.flatten())
    overwrite_preferred = {}
    for i in range(len(query)):
      if (query[i].replace("-","").lower() == preferred[i].replace("-","").lower()) and (query[i].lower() != preferred[i].lower()):
        overwrite_preferred.update({preferred[i].lower():query[i]})
      string_mapping.update({query[i]:preferred[i]}) 

    data = {
    'identifiers': string_list_input,
    'species': species,
    'additional_network_nodes': str(limit),
    'required_score': str(score)
    }
    response = requests.post('http://string-db.org/api/tsv-no-header/interactions', data=data)
    response.raise_for_status()
    string_db_out.update(pd.read_table(StringIO(response.text), header=None))
    string1 = list(string_db_out[2].values.flatten())
    string2 = list(string_db_out[3].values.flatten())
  
    string_interaction,string_unique_vertex = remove_duplicates_within(string1,string2)
 
    no_mapping = [i for i in each_inp_list if i.lower() not in [x.lower() for x in string_mapping]]
    no_mapping_prots =  get_query_from_list(uniprot_query, no_mapping)

    string_interaction,string_unique_vertex,string_dupe_preferred,merged_out_dict,dup_pref_warning = check_dup_preferred_gene(string_mapping,string_interaction,string_unique_vertex,"String", uniprot_query, cy_debug, logging, merged_out_dict)
    
    no_interactions = [i for i in list(string_mapping.keys()) if string_mapping[i].lower() not in [x.lower() for x in string_unique_vertex] and string_mapping[i].lower() not in [y.lower() for y in string_dupe_preferred]]
    no_interactions_prots = get_query_from_list(uniprot_query, no_interactions)
    
    if overwrite_preferred:
      string_interaction,string_unique_vertex,string_mapping,merged_out_dict = overwrite_values(overwrite_preferred, string_interaction, string_unique_vertex, string_mapping, merged_out_dict)
  
    if len(genes_before_initial_drop) != len(each_inp_list):
      string_interaction,string_unique_vertex = drop_ei_if_query(string_interaction,string_unique_vertex,genes_before_initial_drop,each_inp_list)
    
    if cy_debug: 
      if len(no_mapping) != 0:
        logging.debug("DISCARD WARNING - String queries not mapped: " + str(len(no_mapping) + len(dup_pref_warning)))
        warning = []
        for each in no_mapping_prots:
          warning.append(each + "(" + no_mapping_prots[each] + ") ")
        if dup_pref_warning:
          logging.debug("Dropping queries: " + ','.join(warning) + "," + ','.join(dup_pref_warning))
        else:
          logging.debug("Dropping queries: " + ','.join(warning))
    
      if len(no_interactions) != 0:
        logging.debug("DISCARD WARNING - String queries with no interactions: " + str(len(no_interactions)))
        warning = []
        for each in no_interactions_prots:
          warning.append(each + "(" + no_interactions_prots[each] + ") ")
        logging.debug("Dropping queries: " + ','.join(warning))
  except requests.exceptions.HTTPError:
    eprint("Error: String is not responding. Please try again later")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  except:
    return(string_interaction, string_unique_vertex, string_mapping, merged_out_dict)
  return(string_interaction, string_unique_vertex, string_mapping, merged_out_dict)
  
def remove_duplicates_within(node1,node2):
  ''' Remove duplicate interactions ''' 
  unique_interaction_list = []
  unique_vertex_list = []
  for i in range(len(node1)):     
    inter = node1[i] + " " + node2[i]
    backwards = node2[i] + " " + node1[i]
    if (inter not in unique_interaction_list) and (backwards not in unique_interaction_list):
      unique_interaction_list.append(inter)
    
    if node1[i] not in unique_vertex_list:
      unique_vertex_list.append(node1[i])
    if node2[i] not in unique_vertex_list:
      unique_vertex_list.append(node2[i])
  
  return(unique_interaction_list,unique_vertex_list)

def create_genemania_interactions(uniprot_query,each_inp_list,species,limit,att_limit,cy_debug,logging,merged_out_dict,cy_session, cy_out, cy_cluego_out, genes_before_initial_drop, path_to_new_dir, logging_file, cy_settings_file):
  '''
  Obtain gene-gene interactions from Genemania database via Cytoscape App request call:
    1. Get mapping table and categorize genes as primary, secondary, external synonym or external interactor.
    2. Duplicates in preferred gene names are also dropped.
    3. Unique list of interactions is returned.
  '''
  genemania_interaction = []
  genemania_unique_vertex = []
  genemania_mapping = {}
  genemania_node = {}
  overwrite_preferred = {}
  
  join_genes = '|'.join(each_inp_list)
  #body = dict(attrLimit=str(att_limit), geneLimit=str(limit), genes=join_genes, organism=species)
  body = dict(attrLimit=str(att_limit), geneLimit=str(limit), genes=join_genes, organism=species, offline=True)
  try:
    get_genemania = request_retry(f'{CYREST_URL}/v1/commands/genemania/search', 'POST', json=body, timeout=3, request_timeout=10000)
  except CytoscapeError as e:
    if "500 Server Error: Internal Server Error" in str(e):
      return (genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict)
    else:
      raise e
   
  try:    
    uploaded_list = get_genemania.json()
    current_network_suid = str(uploaded_list['data']['network'])
    request = f'{CYREST_URL}/v1/networks/' + str(uploaded_list['data']['network']) + '/tables/defaultedge'
    resp = request_retry(request, 'GET', json=body)
    edge_info =resp.json()
    request = f'{CYREST_URL}/v1/networks/' + str(uploaded_list['data']['network']) +'/tables/defaultnode'
    resp = request_retry(request, 'GET', json=body)
    node_info = resp.json()
    request_retry(f"{CYREST_URL}/v1/networks/" + str(current_network_suid), 'DELETE')
  except:
    return(genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict)
    
  for each_node in node_info['rows']:
    genemania_node.update({each_node['shared name']:each_node['gene name']})   
    if 'query term' in each_node:
      if (each_node['query term'].replace("-","").lower() == each_node['gene name'].replace("-","").lower()) and (each_node['query term'].lower() != each_node['gene name'].lower()):
        overwrite_preferred.update({each_node['gene name']:each_node['query term']})
      genemania_mapping.update({each_node['query term']:each_node['gene name']})
         
  for each_edge in edge_info['rows']:
    edges = each_edge['shared name'].split('|')
    interaction = genemania_node[edges[0]] + " " + genemania_node[edges[1]]      

    if interaction not in genemania_interaction:
      split_interaction = interaction.split(" ")
      backwards = split_interaction[1] + " " + split_interaction[0]
      if backwards not in genemania_interaction:       
        genemania_interaction.append(interaction)
      
      if split_interaction[1] not in genemania_unique_vertex :
        genemania_unique_vertex.append(split_interaction[1])
      if split_interaction[0] not in genemania_unique_vertex :
        genemania_unique_vertex.append(split_interaction[0])
  
  no_mapping = [i for i in each_inp_list if i.lower() not in [x.lower() for x in genemania_mapping]]
  no_mapping_prots =  get_query_from_list(uniprot_query, no_mapping)
  
  no_interactions = [i for i in list(genemania_mapping.keys()) if genemania_mapping[i].lower() not in [x.lower() for x in genemania_unique_vertex] ]
  no_interactions_prots = get_query_from_list(uniprot_query, no_interactions)
  
  #genemania_interaction,genemania_unique_vertex,genemania_dupe_preferred,merged_out_dict = check_dup_preferred_gene(genemania_mapping,genemania_interaction,genemania_unique_vertex,"Genemania", uniprot_query, cy_debug, logging, merged_out_dict)
  if overwrite_preferred:
    genemania_interaction,genemania_unique_vertex,genemania_mapping,merged_out_dict = overwrite_values(overwrite_preferred, genemania_interaction,genemania_unique_vertex,genemania_mapping,merged_out_dict)
  
  if len(genes_before_initial_drop) != len(each_inp_list):
    genemania_interaction,genemania_unique_vertex = drop_ei_if_query(genemania_interaction,genemania_unique_vertex,genes_before_initial_drop,each_inp_list)  
    
  if cy_debug:
    if len(no_mapping) != 0:  
      logging.debug("DISCARD WARNING - Genemania queries not mapped: " + str(len(no_mapping)))
      warning = []
      for each in no_mapping_prots:
        warning.append(each + "(" + no_mapping_prots[each] + ") ")
      logging.debug("Dropping queries: " + ','.join(warning))
    
    if len(no_interactions) != 0:
      logging.debug("DISCARD WARNING - Genemania queries with no interactions: " + str(len(no_interactions)))
      warning = []
      for each in no_interactions_prots:
        warning.append(each + "(" + no_interactions_prots[each] + ") ")
      logging.debug("Dropping queries: " + ','.join(warning))
  
  return(genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict)

def categorize_gene(unique_vertex, mapping, uniprot_query):
    ''' 
    Categorize genes from interaction databases as primary gene, secondary gene, external synonym, external interactor
      1. Primary gene : Uniprot gene name is called by the same name by interaction database
      2. Secondary gene: Uniprot gene name is called by a different name by interaction database and is considered a synonym gene by Uniprot
      3. External synonym: Uniprot gene name is called by a different name by interaction database and is not considered a synonym gene by Uniprot
      4. External Interactor: Gene is part of interactions but is not a part of the input list
    '''
    categories = {}
    for node in unique_vertex:     
      line = ""
      node_category = ""
      query = ""
      if node in list(mapping.values()):
        if node not in categories:
          # If gene name = preferred name in string mapping table -> Primary gene OR Secondary gene OR external synonym
          for prot_id, prot_info in uniprot_query.items():
            # If gene name in string mapping table AND gene name is a primary gene according to uniprot   
            if prot_info['Primary'].lower() == node.lower():
              node_category = "Primary gene"
              categories.update({node:node_category})
              break
                         
            elif node.lower() in (name.lower() for name in prot_info['Synonym']):
              # If gene name in string mapping table AND gene name is a synonym according to uniprot AND gene name not primary gene -> Secondary gene
              node_category = "Secondary gene"
              categories.update({node:node_category})
        
          if node not in categories:        
            # If gene name in string mapping table AND gene name not a primary gene or synonym according to uniprot -> External synonym
            node_category = "External Synonym"
            categories.update({node:node_category})          
        
        # Query = node; for primary gene  & query = ""; for external interactor
        count_of_occurences = list(mapping.values()).count(node)
        query = list(mapping.keys())[list(mapping.values()).index(node)]
    
      else:
        # If gene name not in string mapping table -> external interactor
        node_category = "External Interactor"
        categories.update({node:node_category})

    return(categories)

def get_search_dicts(interaction, categories, logging, cy_debug, uniprot_query, mapping, merged_out_dict, search):
  ''' For list of interactions from interaction databases, retain only interactions categorized as primary gene-primary gene & primary gene-external interactors '''
  search_dict = {}
  categorized_interactions = {}
  dropped_nodes = []
  dropped_interaction_count = 0 
  primary_nodes = []
  count_retained_interactions = 0
  
  for each_interaction in interaction:
    nodes_in_list =  each_interaction.split(" ")
    node1 = nodes_in_list[0]
    node2 = nodes_in_list[1]
    ## Store all interactions based on category of interaction (P-P, P-S, S-P, P-ES, ES-P, P-EI, EI-P, S-ES, ES-S, S-EI, EI-S, ES-EI, EI-ES, S-S, ES-ES, EI-EI
    interaction_cat = categories[node1] + "-" + categories[node2]
    reverse_interaction_cat = categories[node2] + "-" + categories[node1]
    if interaction_cat in categorized_interactions:
      categorized_interactions[interaction_cat].append(each_interaction)
    elif reverse_interaction_cat in categorized_interactions:
      categorized_interactions[reverse_interaction_cat].append(each_interaction)
    else:
      categorized_interactions.update({interaction_cat:[each_interaction]})
    
    # Get only interactions P-P, P-EI, EI-P, EI-EI individually for string and genemania
    if (categories[node1] == "Primary gene" or categories[node1] == "External Interactor") and (categories[node2] == "Primary gene" or categories[node2] == "External Interactor"):
      count_retained_interactions += 1
      if categories[node1] == "Primary gene" and node1 not in primary_nodes:
        primary_nodes.append(node1)
      if categories[node2] == "Primary gene" and node2 not in primary_nodes:
        primary_nodes.append(node2)
        
      if node1 in search_dict:
        if node2 not in search_dict[node1]:
          search_dict[node1].append(node2)
      else:
        search_dict.update({node1:[node2]})
        
      if node2 in search_dict:
        if node1 not in search_dict[node2]:
          search_dict[node2].append(node1)
      else:
        search_dict.update({node2:[node1]})
    else:
      if categories[node1] != "External Interactor" and node1 not in dropped_nodes:
        dropped_nodes.append(node1)
      if categories[node2] != "External Interactor" and node2 not in dropped_nodes:
        dropped_nodes.append(node2)
      dropped_interaction_count += 1
      
  filtered_dropped_nodes = [x for x in dropped_nodes if x not in primary_nodes]
  get_dropped_query_primgenes = {}
  for each_drop in filtered_dropped_nodes:
    for k,v in mapping.items():
      if v.lower() == each_drop.lower():
        get_dropped_query_primgenes.update({k:v})
  filtered_dropped_nodes_prot = get_query_from_list(uniprot_query, get_dropped_query_primgenes.keys())
  
  if cy_debug:
    for each_categorical_interaction in categorized_interactions:
      logging.debug(str(each_categorical_interaction) + " interactions: " + str(len(categorized_interactions[each_categorical_interaction])))
    if len(filtered_dropped_nodes) !=0:
      logging.debug("Only retaining interactions of category: Primary gene-Primary gene, Primary gene-External Interactor, External Interactor-External Interactor")
      logging.debug("DISCARD WARNING -" + str(search) + " dropped interaction category query nodes: " + str(len(filtered_dropped_nodes)))
      warning = []
      for each in filtered_dropped_nodes_prot:
        warning.append(each + "(" + filtered_dropped_nodes_prot[each] + "->" + get_dropped_query_primgenes[filtered_dropped_nodes_prot[each]] + ") ")
      logging.debug("Dropping queries: " + ','.join(warning))
    logging.debug(str(search) + " Total query nodes: " + str(len(primary_nodes)))
    logging.debug(str(search) + " Total interactions: " + str(count_retained_interactions))
    
  return(search_dict, merged_out_dict)  

def get_merged_interactions(filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, type):
  ''' Merge interactions from String and Genemania '''
  unique_merged_interactions = set(unique_merged_interactions)
  unique_nodes = set(unique_nodes)
  for each_node1 in filtered_dict:
    for each_node2 in filtered_dict[each_node1]:
      interaction = each_node1 + " " + each_node2
      backwards = each_node2 + " " + each_node1
      if (interaction.lower() not in unique_merged_interactions) and (backwards.lower() not in unique_merged_interactions):
        unique_merged_interactions.add(interaction.lower())
      
      if each_node2.lower() not in unique_nodes:
        unique_nodes.add(each_node2.lower())
        
    if each_node1.lower() not in unique_nodes:
      unique_nodes.add(each_node1.lower())
      
  return(list(unique_merged_interactions), list(unique_nodes))
  
def get_everything_together(each,uniprot_query, uniprot_list, max_FC_len, each_category, type, site_dict, ambigious_sites, ambigious_genes):
  ''' For genes in merged interaction list, append other info: FC, pval, category etc '''  
  not_in_list = 1
  if not uniprot_list:
    uniprot_list.update({'name':[each]})
    length_each = len(each)
    uniprot_list.update({'length':[length_each]})
    not_in_list = 0
    if each.lower() in [a.lower() for a in ambigious_genes]:
      ambi_gene = each + "**"
      uniprot_list.update({'ambigious_genes': [ambi_gene]})
    else:
      uniprot_list.update({'ambigious_genes': [each]})
  else:
    uniprot_list['name'].append(each)
    length_each = len(each)
    uniprot_list['length'].append(length_each)
    if each.lower() in [a.lower() for a in ambigious_genes]:
      ambi_gene = each + "**"
      uniprot_list['ambigious_genes'].append(ambi_gene)
    else:
      uniprot_list['ambigious_genes'].append(each)
   
  prot_as_list = get_query_from_list(uniprot_query, [each])
  if prot_as_list:
    prot_val = list(prot_as_list.keys())[0]
  
  if type == "1" or type == "2":
    if prot_as_list:
      is_significant = 0
      is_FC_true = 0.0
      if not_in_list == 0:
        count_FC_uniprot = 1
        for each_FC,each_pval in zip(uniprot_query[prot_val]['FC'],uniprot_query[prot_val]['PVal']):
          if float(each_pval) < 0.05:
            is_significant = 1
            each_pval = float(each_pval)
          term_FC = 'FC' + str(count_FC_uniprot)
          term_pval = 'pval' + str(count_FC_uniprot)
          uniprot_list.update({term_FC:[float(each_FC)],term_pval:[each_pval]})
          if (float(each_FC) > 0 or float(each_FC) < 0):
            is_FC_true = 1.0
          count_FC_uniprot+=1
        uniprot_list.update({'query':[each],'significant':[is_significant],'FC_exists':[is_FC_true]})
        not_in_list +=1
      else:
        count_FC_uniprot = 1
        for each_FC,each_pval in zip(uniprot_query[prot_val]['FC'],uniprot_query[prot_val]['PVal']):    
          if float(each_pval) < 0.05:
            is_significant = 1
            each_pval = float(each_pval)
          term_FC = 'FC' + str(count_FC_uniprot)
          term_pval = 'pval' + str(count_FC_uniprot)
          uniprot_list[term_FC].append(float(each_FC))
          uniprot_list[term_pval].append(each_pval)
          if (float(each_FC) > 0 or float(each_FC) < 0):
            is_FC_true = 1.0
          count_FC_uniprot+=1
        uniprot_list['query'].append(each)
        uniprot_list['significant'].append(is_significant)
        uniprot_list['FC_exists'].append(is_FC_true)
    
    else:
      if not_in_list == 0:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list.update({term_FC:[0.0], term_pval:[1.0]})
        uniprot_list.update({'query':['NA'],'significant':[0],'FC_exists':[0]})
        not_in_list +=1
      else:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list[term_FC].append(0.0)
          uniprot_list[term_pval].append(1.0)
        uniprot_list['query'].append('NA')
        uniprot_list['significant'].append(0)
        uniprot_list['FC_exists'].append(0)
  
  elif type == "5" or type == "6":
    if prot_as_list:
      countiter = 0     
      for each_site in site_dict[prot_val]:
        is_significant = 0
        is_FC_true = 0.0
        
        ambigious = ""
        if prot_val in ambigious_sites:
          if each_site in ambigious_sites[prot_val]:
            ambigious = each_site + "**"
        if not ambigious:
          ambigious = each_site
          
        if countiter!=0:
          uniprot_list['name'].append(each)
          if each.lower() in [a.lower() for a in ambigious_genes]:
            ambi_gene = each + "**"
            uniprot_list['ambigious_genes'].append(ambi_gene)
          else:
            uniprot_list['ambigious_genes'].append(each)
          uniprot_list['length'].append(len(each))
          
        if not_in_list == 0:
          count_FC_uniprot = 1
          for each_FC,each_pval in zip(site_dict[prot_val][each_site][0],site_dict[prot_val][each_site][1]):
            if float(each_pval) < 0.05:
              is_significant = 1
              each_pval = float(each_pval)
            term_FC = 'FC' + str(count_FC_uniprot)
            term_pval = 'pval' + str(count_FC_uniprot)
            uniprot_list.update({term_FC:[float(each_FC)],term_pval:[each_pval]})
            if (float(each_FC) > 0 or float(each_FC) < 0):
              is_FC_true = 1.0
            count_FC_uniprot+=1
          uniprot_list.update({'query':[each],'significant':[is_significant],'FC_exists':[is_FC_true],'site':[each_site+"-"+each], 'ambigious_site':[ambigious]})
          not_in_list +=1
        else:
          count_FC_uniprot = 1
          for each_FC,each_pval in zip(site_dict[prot_val][each_site][0],site_dict[prot_val][each_site][1]):    
            if float(each_pval) < 0.05:
              is_significant = 1
              each_pval = float(each_pval)
            term_FC = 'FC' + str(count_FC_uniprot)
            term_pval = 'pval' + str(count_FC_uniprot)
            uniprot_list[term_FC].append(float(each_FC))
            uniprot_list[term_pval].append(each_pval)
            if (float(each_FC) > 0 or float(each_FC) < 0):
              is_FC_true = 1.0
            count_FC_uniprot+=1
          uniprot_list['query'].append(each)
          uniprot_list['significant'].append(is_significant)
          uniprot_list['FC_exists'].append(is_FC_true)
          uniprot_list['site'].append(each_site+"-"+each)
          uniprot_list['ambigious_site'].append(ambigious)
        countiter +=1
        
    else:
      if not_in_list == 0:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list.update({term_FC:[0.0], term_pval:[1.0]})
        uniprot_list.update({'query':['NA'],'significant':[0],'FC_exists':[0],'site':['NA'],'ambigious_site':['NA']})
        not_in_list +=1
      else:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list[term_FC].append(0.0)
          uniprot_list[term_pval].append(1.0)
        uniprot_list['query'].append('NA')
        uniprot_list['significant'].append(0)
        uniprot_list['FC_exists'].append(0)  
        uniprot_list['site'].append('NA')
        uniprot_list['ambigious_site'].append('NA')
  
  else:
    if max_FC_len == 0:
      if prot_as_list:
        if not_in_list == 0:
          uniprot_list.update({'query':[each]})
          not_in_list+= 1
        else:
          uniprot_list['query'].append(each)
      else:
        if not_in_list == 0:
          uniprot_list.update({'query':["NA"]})
          not_in_list+= 1
        else:
          uniprot_list['query'].append("NA")
  
  if type == "4":
    is_category_true = 0.0
    if prot_as_list:
      for i in each_category:
        if i in uniprot_query[prot_val]['Category']:
          if i not in uniprot_list:
            uniprot_list.update({i:[1.0]})
          else:
            uniprot_list[i].append(1.0)
          is_category_true = 1.0  
        else:
          if i not in uniprot_list:
            uniprot_list.update({i:[0.0]})
          else:
            uniprot_list[i].append(0.0)
    else:
      for i in each_category:
        if i not in uniprot_list:
          uniprot_list.update({i:[0.0]})
        else:
          uniprot_list[i].append(0.0)
          
    if "category_true" not in uniprot_list:
      uniprot_list.update({"category_true":[is_category_true]})
    else:
      uniprot_list["category_true"].append(is_category_true)
    
  return(uniprot_list)
      
def remove_list_duplicates(list):
  ''' Return only elements from the list where the lower case version is seen multiple times '''
  dupe_gene_list = []
  set_list = []
  for x in list:
    if x.lower() in [set_name.lower() for set_name in set_list]:
      if x.lower() not in [dupe_name.lower() for dupe_name in dupe_gene_list]:
        dupe_gene_list.append(x)
    else:
      set_list.append(x)
  return(dupe_gene_list)

def cluego_filtering(unique_nodes, cluego_mapping_file, uniprot_query, cy_debug, logging, merged_out_dict, unique_each_primgene_list, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file):
  '''
  Perform ClueGO input mapping list of genes via a Cytoscape App request call:
    1. Only retain ClueGo primary genes
    2. Drop all duplicate preferred genes
  '''    
  # Assumption made: ClueGO file is sorted in increasing order of uniqueid 
  acceptable_genes_list = {}
  each_preferred_gene = []   
  with gzip.open(cluego_mapping_file,'rt') as csv_file:
    input_file = csv.reader(csv_file, delimiter='\t')   
    line_count = 0
    is_sym_col = False
    is_uniq_col = False
    for row in input_file:
      if line_count == 0:
        for i in range(len(row)):
          if "symbolid" in row[i].lower():
            symbolid = i
            is_sym_col = True
          if "uniqueid#entrezgeneid" in row[i].lower():
            uniqueid = i
            is_uniq_col = True            
      else:
        if not is_sym_col or not is_uniq_col:
          eprint("Error: Required columns in ClueGO mapping file: SymbolID and UniqueID#EntrezGeneID")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
          sys.exit(1)
        primgene = row[symbolid]
        secgene = []        
        if "|" in row[symbolid]:
          secgene = primgene.split("|")[1:]
          primgene = primgene.split("|")[0]
        if primgene.lower() in (name.lower() for name in unique_nodes):
          if primgene.lower() not in acceptable_genes_list:
            acceptable_genes_list[primgene.lower()] = {}
            acceptable_genes_list[primgene.lower()].update({'GeneCategory':"Primary",'GeneUniqueID':row[uniqueid],"GenePreferredName":primgene})
          elif primgene.lower() in acceptable_genes_list and acceptable_genes_list[primgene.lower()]["GeneCategory"] != "Primary":
            acceptable_genes_list[primgene.lower()]["GeneCategory"] = "Primary"
            acceptable_genes_list[primgene.lower()]["GeneUniqueID"] = row[uniqueid]
            acceptable_genes_list[primgene.lower()]["GenePreferredName"] = primgene
        for each_secondary_gene in secgene:
          if each_secondary_gene.lower() in (name.lower() for name in unique_nodes):
            if each_secondary_gene.lower() not in acceptable_genes_list:
              acceptable_genes_list[each_secondary_gene.lower()] = {}
              acceptable_genes_list[each_secondary_gene.lower()].update({"GeneCategory":"Other","GeneUniqueID":row[uniqueid],"GenePreferredName":primgene})                          
      line_count +=1
  
  not_found_in_cluego = [i for i in unique_nodes if i.lower() not in acceptable_genes_list ] 
  not_found_query = [i for i in not_found_in_cluego if i.lower() in [name.lower() for name in unique_each_primgene_list]]
  not_found_query_prot = get_query_from_list(uniprot_query, not_found_query)
  not_found_ei = [i for i in not_found_in_cluego if i.lower() not in [name.lower() for name in unique_each_primgene_list]]

  if len(not_found_query) > 0 or len(not_found_ei) > 0:    
    warning1 = []
    warning2 = not_found_ei
    for each in not_found_query_prot:
      warning1.append(each + "(" + not_found_query_prot[each] + ") ")
      if 'CommentGene' in merged_out_dict[each]:
        merged_out_dict[each]['CommentGene'] += "Cluego- not found;"
      else:
        merged_out_dict[each].update({'CommentGene':"Cluego- not found;"})
    if cy_debug and (warning1 or warning2):
      logging.debug("DISCARD WARNING - ClueGO query + External Interactor unavailable: " + str(len(not_found_query)) + " + " + str(len(not_found_ei)))
      #if warning1:
        #logging.warning("Dropped queries: " + ','.join(warning1))
      #if warning2:
        #logging.warning("Dropped External Interactor: " + ','.join(warning2))
 
  # Drop queries with duplicate preferred gene list
  for each_acceptable in acceptable_genes_list:
    each_preferred_gene.append(acceptable_genes_list[each_acceptable]['GenePreferredName'])
  
  get_dupe_query_primgenes = {}
  filtered_unique_list = [i for i in acceptable_genes_list if i.lower() not in [name.lower() for name in get_dupe_query_primgenes]]
  
  drop_non_primary = {}
  for each_node in filtered_unique_list:
    if acceptable_genes_list[each_node]['GeneCategory'] != 'Primary':
      drop_non_primary.update({each_node:acceptable_genes_list[each_node]['GenePreferredName']})
  
  not_found_query = [i for i in drop_non_primary if i.lower() in [name.lower() for name in unique_each_primgene_list]]
  get_drop_query_proids = get_query_from_list(uniprot_query, not_found_query)
  not_found_ei = [i for i in drop_non_primary if i.lower() not in [name.lower() for name in unique_each_primgene_list]]
  
  warning = []
  for each_in_list in get_drop_query_proids:
    if 'CommentGene' in merged_out_dict[each_in_list]:
      merged_out_dict[each_in_list]['CommentGene'] += "ClueGO non-primary query;"
    else:
      merged_out_dict[each_in_list].update({'CommentGene':"ClueGO non-primary query;"})
    warning.append(each_in_list + "(" + get_drop_query_proids[each_in_list] + "->" + drop_non_primary[get_drop_query_proids[each_in_list]] + ") ")
    
  if cy_debug:
    if len(not_found_query) != 0 or len(not_found_ei) != 0:
      logging.debug("DISCARD WARNING - ClueGO non-primary query + External Interactor genes: " + str(len(not_found_query)) + " + " + str(len(not_found_ei)))
      #if warning:
        #logging.warning("Dropping queries: " + ','.join(warning))
      #if not_found_ei:
        #logging.warning("Dropping External Interactor: " + ','.join(not_found_ei))  
  
  filtered_preferred_unique_list = [acceptable_genes_list[i]['GenePreferredName'] for i in filtered_unique_list if i.lower() not in [x.lower() for x in drop_non_primary]]
  
  return(filtered_preferred_unique_list,merged_out_dict)

def calc_protein_change_sf(df, uniprot_list, type):
  cluego_dict = df.to_dict()
  add_col_to_end = {}
  for each_term in cluego_dict:
    if "associated genes found" in each_term.lower():
      for each_col in cluego_dict[each_term]:
        up_or_down = {}
        total_genes = 0
        calc_up = 0
        calc_down = 0
        list_of_genes = cluego_dict[each_term][each_col].strip('][').split(', ') 
        for each_gene in list_of_genes:
          if type == "5":
            if each_gene.lower() in uniprot_list['name']:
              indices = [i for i, x in enumerate(uniprot_list['name']) if x == each_gene.lower()]
              for each_index in indices:
                total_genes += 1                  
                FC_val_each_gene = (uniprot_list['FC1'])[each_index]
                if FC_val_each_gene > 0:
                  calc_up += 1
                elif FC_val_each_gene < 0:
                  calc_down -= 1
            else:
              total_genes+=1         
              
          elif type == "1":
            total_genes += 1          
            if each_gene.lower() in uniprot_list["name"]:
              indexOf = uniprot_list["name"].index(each_gene.lower())
              FC_val_each_gene = (uniprot_list['FC1'])[indexOf]
              if FC_val_each_gene > 0:
                calc_up += 1
              elif FC_val_each_gene < 0:
                calc_down -= 1                
        
        if (calc_up/total_genes*100) > 60: 
          percent_val = str(round((calc_up/total_genes*100),2))
        elif (abs(calc_down)/total_genes*100) > 60:
          percent_val = str(round((calc_down/total_genes*100),2))
        else: 
          percent_val = "0"
        up_or_down.update({"Percent":percent_val})
        
        add_new_col = json.dumps([{'label': 'Status', 'percent': v} for k,v in up_or_down.items()])
        #add_new_col = json.dumps([{'label': k, 'percent': v} for k,v in new_dict.items()])
        add_col_to_end[each_col] = add_new_col
  cluego_dict.update({'Status':add_col_to_end})
  df = pd.DataFrame.from_dict(cluego_dict)
  return(df)

def calc_protein_change_mf(df, uniprot_list, type, max_FC_len, unique_labels):
  cluego_dict = df.to_dict()
  add_col_to_end = {}
  for each_term in cluego_dict:
    if "associated genes found" in each_term.lower():
      for each_col in cluego_dict[each_term]:
        up_or_down = {}
        list_of_genes = cluego_dict[each_term][each_col].strip('][').split(', ') 
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          val_term_FC = unique_labels[i-1]
          total_genes = 0
          calc_up = 0
          calc_down = 0
          for each_gene in list_of_genes:           
            if type == "6":   
              if each_gene.lower() in uniprot_list['name']:
                indices = [i for i, x in enumerate(uniprot_list['name']) if x == each_gene.lower()]
                for each_index in indices:
                  total_genes += 1                  
                  FC_val_each_gene = (uniprot_list['FC1'])[each_index]
                  if FC_val_each_gene > 0:
                    calc_up += 1
                  elif FC_val_each_gene < 0:
                    calc_down -= 1
              else:
                total_genes+=1
                
            elif type == "2":
              total_genes += 1            
              if each_gene.lower() in uniprot_list["name"]:
                indexOf = uniprot_list["name"].index(each_gene.lower())
                FC_val_each_gene = (uniprot_list[term_FC])[indexOf]
                if FC_val_each_gene > 0:
                  calc_up += 1
                elif FC_val_each_gene < 0:
                  calc_down -= 1   
                  
          if (calc_up/total_genes*100) > 60: 
            percent_val = str(round((calc_up/total_genes*100),2))
          elif (abs(calc_down)/total_genes*100) > 60:
            percent_val = str(round((calc_down/total_genes*100),2))
          else: 
            percent_val = "0"
          up_or_down.update({val_term_FC:percent_val})
          
        add_new_col = json.dumps([{'label': k, 'percent': v} for k,v in up_or_down.items()])
        add_col_to_end[each_col] = add_new_col
        
  cluego_dict.update({'Status':add_col_to_end})
  df = pd.DataFrame.from_dict(cluego_dict)
  return(df)        
                  
SEP = "/"
HEADERS = {'Content-Type': 'application/json'}

CLUEGO_BASE_PATH = "/v1/apps/cluego/cluego-manager"

def writeLines(lines, out_file, type, uniprot_list, max_FC_len, unique_labels):
  ''' Write the lines to a file, removing duplicates of specific columns '''
  df = pd.read_csv(StringIO(lines), sep='\t')
  df = df.drop_duplicates(subset=['GOTerm','Ontology Source', 'Nr. Genes', 'Associated Genes Found'])
  if type == "1" or type == "5":
    df = calc_protein_change_sf(df, uniprot_list, type)
  if type == "2" or type == "6":
    df = calc_protein_change_mf(df, uniprot_list, type, max_FC_len, unique_labels)
    
  df.to_csv(out_file, header=True, index=False, sep='\t', mode='w')
 
def writeBin(raw,out_file):
  ''' Write binary data to a file '''
  file = open(out_file,'wb')
  file.write(raw)
  file.close()
  
def writeLog(lines,out_file):
    file = open(out_file,'w')
    for line in lines:
        file.write(line)
    file.close()
    
def cluego_run(organism_name,output_cluego,merged_vertex,group,select_terms, leading_term_selection, reference_file,cluego_pval, cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file, cy_cluego_log_out, type, uniprot_list, max_FC_len, unique_labels):
  '''
  Obtain ClueGO annotations based on user settings for the list of genes via a Cytoscape App request call. Output is a list of annotation terms along with associated genes and corresponding term p-value
  ''' 
  if group.lower() == "global":
    if len(merged_vertex) <= 500:
      min_number_of_genes_per_term = 3
    else:
      min_number_of_genes_per_term = 20
    min_percentage_of_genes_mapped = 0
    min_go = 1
    max_go = 6
    kappa = 0.5
  elif group.lower() == "medium":
    if len(merged_vertex) <= 500:
      min_number_of_genes_per_term = 3
      min_percentage_of_genes_mapped = 2
    else:
      min_number_of_genes_per_term = 20
      min_percentage_of_genes_mapped = 4
    min_go = 7
    max_go = 11
    kappa = 0.4
  else:
    min_number_of_genes_per_term = 1
    min_percentage_of_genes_mapped = 50
    min_go = 7
    max_go = 15
    kappa = 0.4
    
  #### Select the ClueGO Organism to analyze ####
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"organisms"+SEP+"set-organism"+SEP+str(organism_name), "PUT", headers=HEADERS)
  
  ## Use custom reference file
  if reference_file:
    if not path.exists(reference_file):
      eprint("Error: Path to ClueGO Reference file " + reference_file + " does not exist")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)
    response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"stats/Enrichment%2FDepletion%20(Two-sided%20hypergeometric%20test)/Bonferroni%20step%20down/false/false/true/"+reference_file, "PUT")
  
  # Set the number of Clusters
  number = 1
  max_input_panel_number = number
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"cluster"+SEP+"max-input-panel"+SEP+str(max_input_panel_number), "PUT")
  
  cluster = 1
  gene_list = json.dumps(merged_vertex)
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"cluster"+SEP+"upload-ids-list"+SEP+quote(str(cluster)), "PUT", data=gene_list, headers=HEADERS)
  
  # 2.5 Set analysis properties for a Cluster
  input_panel_index = 1
  node_shape = "Ellipse" # ("Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle")
  cluster_color = "#ff0000" # The color in hex, e.g. #F3A455
  
  no_restrictions = False # "True" for no restricions in number and percentage per term
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"cluster"+SEP+"set-analysis-properties"+SEP+str(input_panel_index)+SEP+node_shape+SEP+quote(cluster_color)+SEP+str(min_number_of_genes_per_term)+SEP+str(min_percentage_of_genes_mapped)+SEP+str(no_restrictions), "PUT", headers=HEADERS)
  
  #Select visual style
  visual_style = "ShowClusterDifference"
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"cluster"+SEP+"select-visual-style"+SEP+visual_style, "PUT", headers=HEADERS)
  
  ## 3.1 Get all available Ontologies
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"get-ontology-info", "GET", headers=HEADERS)
  ontology_info = latest_ontologies(response.json())
  
  i = 0
  list_ontology = []
  for each_ontology in ontology_info:
    get_each = ontology_info[each_ontology]
    m = re.search(r"name=([A-Za-z\-\.0-9]+)\,",get_each)
    if m:
      ontologies = m.group(1)
      if select_terms.lower() == "biological process" or select_terms.lower() == "all":
        if "BiologicalProcess" in ontologies:
          list_ontology.append(each_ontology+";"+"Ellipse")
      if select_terms.lower() == "cellular component" or select_terms.lower() == "all":
        if "CellularComponent" in ontologies:
          list_ontology.append(each_ontology+";"+"Ellipse")
      if select_terms.lower() == "molecular function" or select_terms.lower() == "all":
        if "MolecularFunction" in ontologies:
          list_ontology.append(each_ontology+";"+"Ellipse")
      if select_terms.lower() == "pathways" or select_terms.lower() == "all":
        if "Human-diseases" in ontologies or "KEGG" in ontologies or "Pathways" in ontologies or "WikiPathways" in ontologies or "CORUM" in ontologies:
          list_ontology.append(each_ontology+";"+"Ellipse")
    i+=1

  ####Select Ontologies
  selected_ontologies = json.dumps(list_ontology)
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"set-ontologies", "PUT", data=selected_ontologies, headers=HEADERS)
  
  ## 3.1 Set kappa Score
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"set-kappa-score-level"+SEP+str(kappa), "PUT")
  
  ## 3.2 Select Evidence Codes
  evidence_codes = json.dumps(["All"]) # (run "3.3 Get all available Evidence Codes" to get all options)
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"set-evidence-codes", "PUT", data=evidence_codes, headers=HEADERS)

  ## 3.3 Get all available Evidence Codes
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"get-evidence-code-info", "GET", headers=HEADERS)
  
  ## 3.4 Use GO Fusion
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"true", "PUT", headers=HEADERS)
  
  ## 3.5 Set min/max level for GO
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies"+SEP+"set-min-max-levels"+SEP+str(min_go)+SEP+str(max_go)+SEP+"false", "PUT", headers=HEADERS)
  
  ## 3.6 GO Grouping
  # Highest%20Significance, %23Genes%20%2F%20Term, %25Genes%20%2F%20Term, %25Genes%20%2F%20Term%20vs%20Cluster
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"grouping"+SEP+"true"+SEP+"Random"+SEP+leading_term_selection+SEP+"Kappa%20Score"+SEP+str(1)+SEP+str(30)+SEP+str(30), "PUT", headers=HEADERS)
  
  ## 3.7 Set Pval
  response = request_retry(CYREST_URL+CLUEGO_BASE_PATH+SEP+"ontologies/true"+SEP+str(cluego_pval), "PUT")
  
  #### Run ClueGO Analysis ####
  # Run the analysis an save log file
  analysis_name = "ClueGO Network"
  selection = "Continue analysis"                                                                         
  response = requests.get(CYREST_URL+CLUEGO_BASE_PATH+SEP+analysis_name+SEP+selection, headers=HEADERS)
  # Log File
  # log_file_name = cy_cluego_log_out
  # writeLog(response.text,log_file_name)
  try:
    response.raise_for_status()
  except:
    try:
      cluego_msg = response.json()['message']
      if "IOException: Unexpected end of ZLIB input stream" in cluego_msg:
        cluego_msg = "Unexpected end of ZLIB input stream! Please check in the 'ClueGOFiles' directory or delete them and reload ClueGO."
      eprint(f"Error: ClueGO couldn't find any pathways and exited with the following message: {cluego_msg}")
    except:
      eprint("Error: ClueGO couldn't find any pathways")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)  

  # Get network id (SUID) (CyRest function from Cytoscape)
  response = request_retry(CYREST_URL+"/v1"+SEP+"networks"+SEP+"currentNetwork", "GET", headers=HEADERS)
  current_network_suid = response.json()['data']['networkSUID']
    
  # 4.2 Get ClueGO result table
  if output_cluego:   
    response = requests.get(CYREST_URL+CLUEGO_BASE_PATH+SEP+"analysis-results"+SEP+"get-cluego-table"+SEP+str(current_network_suid))
    try:
      response.raise_for_status()
      table_file_name = output_cluego
      writeLines(response.text,table_file_name,type,uniprot_list,max_FC_len,unique_labels)
    except:
      traceback.print_exc()
      eprint("Error: No pathways found for input list")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)  
  
def get_interactions_dict(filtered_dict, search, merged_out_dict):
  ''' Add interactions found in the search database to merged_out_dict '''
  lower_filtered = [name.lower() for name in filtered_dict]
  for each_uniprot_query in merged_out_dict:    
    for name in filtered_dict:
      if 'Primary' in merged_out_dict[each_uniprot_query]:
        if name.lower() == ((merged_out_dict[each_uniprot_query]['Primary'].lower()).replace("**","")):
          if 'CommentGene' not in merged_out_dict[each_uniprot_query]:        
            interaction_list = filtered_dict[name]
            merged_out_dict[each_uniprot_query][search] += ';'.join(interaction_list)
            break
          elif not merged_out_dict[each_uniprot_query]['CommentGene']:
            interaction_list = filtered_dict[name]
            merged_out_dict[each_uniprot_query][search] += ';'.join(interaction_list)
            break            
  return(merged_out_dict)

def write_into_out(merged_out_dict, out, dup_prot_ids):
  ''' Write into interactions.csv containing for the list of input protein ids the list of its interactions, if present or comments describing reason for drop in any of the steps of the pipeline '''  
  if not out:
    return
  with open(out,'a') as csv_file:
    csv_file = open(out,'w')
    i = 0
    for each_prot in merged_out_dict:
      if 'PickedPeptide' in merged_out_dict[each_prot] or 'DroppedPeptide' in merged_out_dict[each_prot]:
        if i == 0:
          csv_file.write("ProteinID,Primary Gene,String,Genemania,Reason for Dropped Gene,Retained Peptide,Retained Site,Dropped Peptide,Dropped Site,Reason for Dropped Site \n")
        if "PickedPeptide" in merged_out_dict[each_prot]:
          ambi_picked_pep = [i for i in merged_out_dict[each_prot]["PickedPeptide"] if "**" in i]
          ambi_picked_site = [i for i in merged_out_dict[each_prot]["PickedSite"] if "**" in i]
          other_picked_pep = [i for i in merged_out_dict[each_prot]["PickedPeptide"] if "**" not in i]
          other_picked_site = [i for i in merged_out_dict[each_prot]["PickedSite"] if "**" not in i]
          picked_pep_list = ambi_picked_pep + other_picked_pep
          picked_site_list = ambi_picked_site + other_picked_site
          
        else: 
          picked_pep_list = [""]
          picked_site_list = [""]
        line = each_prot + "," + merged_out_dict[each_prot].get("Primary","") + "," + merged_out_dict[each_prot].get("String","") + "," + merged_out_dict[each_prot].get("Genemania","") + "," + merged_out_dict[each_prot].get("CommentGene","") + "," + ";".join(picked_pep_list) + ", " + ";".join(picked_site_list) + ", " + ";".join(merged_out_dict[each_prot].get("DroppedPeptide",[""])) + ", " + ";".join(merged_out_dict[each_prot].get("DroppedSite",[""])) + ", " + merged_out_dict[each_prot].get("Comment","")  + "\n"
      else:
        if i == 0:
          csv_file.write("ProteinID,Primary Gene,String,Genemania,Reason for Dropped Gene \n")
        line = each_prot + "," + merged_out_dict[each_prot].get("Primary","") + "," + merged_out_dict[each_prot].get("String","") + "," + merged_out_dict[each_prot].get("Genemania","") + "," + merged_out_dict[each_prot].get("CommentGene","") + "\n"
        
      i+=1
      csv_file.write(line)
      
    if dup_prot_ids:
      for each_prot in dup_prot_ids:
        line = each_prot + ","  + "," + "," + "," + "Duplicate Query" + "\n"
        csv_file.write(line)        
  csv_file.close()
    
def cluego_input_file(cluego_inp_file, cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file):
  ''' Perform reanalysis network based on annotation terms selected by the user. Central node is annotation term and connecting nodes are genes accosiated with annotation term '''
  top_annotations = {}
  unique_gene = []
  with open(cluego_inp_file,'r') as csv_file:
    input_file = csv.reader(csv_file, delimiter='\t')
    line_count = 0
    is_go_term = False
    is_genes_found = False
    for row in input_file:            
      if line_count == 0:
        header = row
        # Get header
        for i in range(len(row)):
          if "goterm" in row[i].lower():
            goterm = i
            is_go_term = True 
          elif "associated genes found" in row[i].lower():
            genes = i
            is_genes_found = True    
            
        if not (is_go_term and is_genes_found):
          eprint("Error: Required columns in cluego input file: GOTerm and Associated Genes Found")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
          sys.exit(1)
      
      else:
        if row:  
          if "[" in row[genes] and "]" in row[genes]:
            each_gene_list = (row[genes])[1:-1]
          each_gene_list = each_gene_list.split(', ')
        
          for each in each_gene_list:
            if each not in unique_gene:
              unique_gene.append(each)
            
          if row[goterm] not in top_annotations:
            top_annotations.update({row[goterm]:each_gene_list})
          else:
            eprint("Error: Duplicate GOTerm found in cluego input file")
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
            sys.exit(1)
          
      line_count+=1
  return(top_annotations, unique_gene)

def cy_sites_interactors_style(merged_vertex, merged_interactions, uniprot_list, max_FC_len, each_category, pval_style, type, all_prot_site_snps, uniprot_query, unique_labels,mult_mods_of_int):
  ''' Styling + visualization for the gene list interaction network & its modification sites '''
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 
  G = igraph.Graph()
  
  site_interactions = []
  sig_na = []
  fc_na = []
  get_each_site_name = []
  query_val = []
  is_snp = []
  all_fcs = {}
  all_pval = {}
  all_names = []
  all_names_without_unimod = []
  all_other_fcs = []
  all_other_pval = []
  all_length = []
  
  for each_name,each_site_gene,each_ambi_site in zip(uniprot_list['name'],uniprot_list['site'],uniprot_list['ambigious_site']):
    if each_site_gene != "NA":
      each_gene = (each_site_gene.split("-"))[1]
      each_site = (each_site_gene.split("-"))[0]
      if each_gene.lower() in merged_vertex:
        # Get corresponding prot, check if snp
        corresponding_prot = get_query_from_list(uniprot_query, [each_gene])
        if corresponding_prot:
          if list(corresponding_prot.keys())[0] in all_prot_site_snps:
            if each_site in all_prot_site_snps[list(corresponding_prot.keys())[0]]:
              is_snp.append(1.0)
            else:
              is_snp.append(0.0)
          else:
            is_snp.append(0.0)
        else:
          is_snp.append(0.0)
      
        indexOf = uniprot_list['site'].index(each_site_gene)
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          if term_FC in all_fcs:
            all_fcs[term_FC].append(uniprot_list[term_FC][indexOf])
          else:
            all_fcs.update({term_FC:[uniprot_list[term_FC][indexOf]]})
        
          if term_pval in all_pval:
            all_pval[term_pval].append(uniprot_list[term_pval][indexOf])
          else:
            all_pval.update({term_pval:[uniprot_list[term_pval][indexOf]]})
        all_names.append(each_ambi_site)
        combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
        without_unimod = re.sub(combined_pat, '', each_ambi_site)
        all_names_without_unimod.append(without_unimod)
        get_each_site_name.append(each_ambi_site)
        query_val.append("Site")
        G.add_vertex(each_site_gene)
        interaction = each_site_gene + " " + each_gene
        site_interactions.append(interaction)
        all_length.append(uniprot_list['length'][indexOf]) 
        fc_na.append(1.0)
        sig_na.append(uniprot_list['significant'][indexOf])
    else:
      if each_name.lower() in merged_vertex:
        all_names.append(each_name)
        all_names_without_unimod.append(each_name)
        query_val.append("EI")
        G.add_vertex(each_name)
        all_length.append(len(each_name))
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          if term_FC in all_fcs:
            all_fcs[term_FC].append(0.0)
          else:
            all_fcs.update({term_FC:[0.0]})
          if term_pval in all_pval:
            all_pval[term_pval].append(1.0)
          else:
            all_pval.update({term_pval:[1.0]})
        fc_na.append(0.0)
        sig_na.append(0)
      
  for each_vertex in merged_vertex:
    if each_vertex in uniprot_list['query']:
      indexOf = uniprot_list['query'].index(each_vertex)
      all_names.append(uniprot_list['ambigious_genes'][indexOf])
      all_names_without_unimod.append(uniprot_list['ambigious_genes'][indexOf])      
      fc_na.append(-1.0)
      all_other_fcs.append(-100)
      all_other_pval.append(1.0)   
      query_val.append("Gene")
      is_snp.append(0.0)
      G.add_vertex(each_vertex)  
      sig_na.append(0)
      all_length.append(len(each_vertex))
  
  count_each = 0
  for each in merged_interactions:
    each_interaction_name = each.split(" ")
    G.add_edge(each_interaction_name[0],each_interaction_name[1],name=each,interaction="pp")
    count_each +=1
  for each in site_interactions:
    each_interaction_name = each.split(" ")
    G.add_edge(each_interaction_name[0],each_interaction_name[1],name=each,interaction="ps")
    count_each +=1
    
  G.vs
  G.vs["name"] = all_names
  G.vs["substitute name"] = all_names_without_unimod
  #G.vs["SNP"] = is_snp
  domain_labels = [] 
  value_labels = []
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    if unique_labels:
      val_term_FC = unique_labels[i-1]
      val_term_pval = unique_labels[i-1] + " pval"
    else:
      val_term_FC = term_FC
      val_term_pval = term_pval
    collect_fcs = all_fcs[term_FC] + all_other_fcs
    domain_labels.append(str(val_term_FC))
    if collect_fcs:
      if not value_labels:
        value_labels = [ [str(collect_fcs[i])] for i in range(0, len(collect_fcs)) ]
      else:
        for i in range(0, len(collect_fcs)):
          value_labels[i].append(str(collect_fcs[i])) 
    G.vs[val_term_FC] = all_fcs[term_FC] + all_other_fcs
    G.vs[val_term_pval] = all_pval[term_pval] + all_other_pval
  
  G.vs["pine_query"] = query_val
  G.vs["pine_length"] = all_length  
  G.vs["pine_outline"] = sig_na
  G.vs["pine_FC_exists"] = fc_na
  degree = G.degree()
  G.vs["pine_degree"] = degree
  if max_FC_len > 1:
    tot = len(query_val)
    G.vs["pine_domain_labels"] = [domain_labels]*tot
    G.vs["pine_value_labels"] = value_labels
  cy = CyRestClient(port=CYREST_PORT)
  g_cy = cy.network.create_from_igraph(G, name="Interaction Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Initial Network Style')
  basic_settings = {
    'NODE_CUSTOMGRAPHICS_1':"org.cytoscape.BarChart",
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'NODE_SIZE':"30",
    'EDGE_WIDTH':"1",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  my_style.update_defaults(basic_settings)
  if not mult_mods_of_int:
    my_style.create_passthrough_mapping(column='substitute name', vp='NODE_LABEL', col_type='String')
  else:
    my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  node_label_size =  [
    {
      "value": "2",
      "lesser": "1",
      "equal": "10",
      "greater": "10"
    },
    {
      "value": "25",
      "lesser": "1",
      "equal": "1",
      "greater": "1"
    }
  ]
  my_style.create_continuous_mapping(column='pine_length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  columns = ['pine_length','pine_degree']

  chart={'cy_range':'[2.0,9.0]','cy_colorScheme':"Custom"}
  chart = json.dumps(chart)
  loaded_r = json.loads(chart)
  just_color = 0
  pval_sig = 0
  slash_delim = chr(92)
  
  shape_kv_pair = {
    "Function":"OCTAGON",
    "Site":"ELLIPSE",
    "Gene":"ROUND_RECTANGLE",
    "EI":"ROUND_RECTANGLE"
  }
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)

  label_color_kv_pair = {
        "Function":"#FFFFFF",
        "Gene":"#000000",
        "Site":"#000000",
        "EI":"#000000"
  }
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)
  
  width_kv_pair = {
    "Function":"210",
    "Gene":"70",
    "EI":"70",
    "Site":"60"
  }
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style, fc_na,query_val,"-1",all_names,max_FC_len, uniprot_list, unique_labels)
    if pval_style:
      pval_sig = 1
    
  elif max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list, type)
    if pval_style:
      pval_sig = 1
   
  if pval_sig == 1:
    my_style = pval(my_style)
  
  my_style = snp_bold(my_style)
  
  cy.style.apply(my_style, g_cy)
    
def cy_interactors_style(merged_vertex, merged_interactions, uniprot_list, max_FC_len, each_category, pval_style, unique_labels):
  ''' Styling + visualization for the gene list interaction network '''
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 

  G = igraph.Graph()
  for each in merged_vertex:
    G.add_vertex(each)
  
  count_each = 0
  edges = []
  for each in merged_interactions:
    each_interaction_name = each.split(" ")
    edges.append((each_interaction_name[0], each_interaction_name[1]))
    count_each +=1
  G.add_edges(edges)

  G.vs 
  G.vs["name"] = uniprot_list['name']
  G.vs["shared name"] = uniprot_list['ambigious_genes']
  
  domain_labels = []
  value_labels = []
  
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    if unique_labels:
      val_term_FC = unique_labels[i-1]
      val_term_pval = unique_labels[i-1] + " pval"
    else:
      val_term_FC = term_FC
      val_term_pval = term_pval
      
    collect_fcs = uniprot_list[term_FC] 
    domain_labels.append(str(val_term_FC))
    if collect_fcs:
      if not value_labels:
        value_labels = [ [str(collect_fcs[i])] for i in range(0, len(collect_fcs)) ]
      else:
        for i in range(0, len(collect_fcs)):
          value_labels[i].append(str(collect_fcs[i]))
    G.vs[val_term_FC] = uniprot_list[term_FC]
    G.vs[val_term_pval] = uniprot_list[term_pval]
  
  G.vs["pine_query"] = uniprot_list["query"]
  G.vs["pine_length"] = uniprot_list['length']
  
  if max_FC_len >= 1:
    G.vs["pine_outline"] = uniprot_list["significant"]
    G.vs["pine_FC_exists"] = uniprot_list["FC_exists"]
  if max_FC_len > 1:
    tot = len(uniprot_list["query"])
    G.vs["pine_domain_labels"] = [domain_labels]*tot
    G.vs["pine_value_labels"] = value_labels
    
  degree = G.degree()
  G.vs["pine_degree"] = degree  
  category_present = 0
  
  for each in each_category:
    category_present = 1
    G.vs[each] = uniprot_list[each]
  if category_present == 1:
    G.vs['pine_category_true'] = uniprot_list['category_true']
  
  cy = CyRestClient(port=CYREST_PORT)
  g_cy = cy.network.create_from_igraph(G, name="Interaction Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Initial Network Style')
  basic_settings = {
    'NODE_CUSTOMGRAPHICS_1':"org.cytoscape.BarChart",
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'NODE_SIZE':"30",
    'EDGE_WIDTH':"1",
    'NODE_SHAPE':"ROUND_RECTANGLE",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  #degree_to_label_size = StyleUtil.create_slope(min=min(degree),max=max(degree),values=(10, 20))
  node_label_size =  [
    {
      "value": "2",
      "lesser": "1",
      "equal": "10",
      "greater": "10"
    },
    {
      "value": "25",
      "lesser": "1",
      "equal": "1",
      "greater": "1"
    }
  ]
  my_style.create_continuous_mapping(column='pine_length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  columns = ['pine_length','pine_degree']

  chart={'cy_range':'[2.0,9.0]','cy_colorScheme':"Custom"}
  chart = json.dumps(chart)
  loaded_r = json.loads(chart)
  just_color = 0
  pval_sig = 0
  slash_delim = chr(92)
  
  if category_present == 1:
    my_style = get_category(my_style,uniprot_list['category_true'],"1",uniprot_list["query"],uniprot_list['name'], each_category)
    color_kv_pair = {
      "1":"#FFFFFF",
      "0":"#E2E2E2"
    }
    my_style.create_discrete_mapping(column='pine_category_true', col_type='Double', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,uniprot_list['FC_exists'],uniprot_list["query"],"1",uniprot_list['name'],max_FC_len, uniprot_list, unique_labels)
    if pval_style:
      pval_sig = 1
    
  elif (max_FC_len == 0 and category_present == 0):
    just_color = 1
    
  elif max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list, "1")
    if pval_style:
      pval_sig = 1
  
  if just_color == 1:
    my_style = color(my_style, uniprot_list)
  
  if pval_sig == 1:
    my_style = pval(my_style)
  
  cy.style.apply(my_style, g_cy)
  
def cy_category_style(merged_vertex, merged_interactions, uniprot_list, each_category):
  ''' Styling specific to category case '''
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 
  G = igraph.Graph()
  
  additional_nodes = []
  additional_len = []
  additional_query = []
  additional_each = []
  category_true = []
  length_of = []
  breadth_of = []
  all_names = []
  category_values = {}
  
  for each in merged_vertex:
    indexOf = uniprot_list['name'].index(each)
    if uniprot_list['query'][indexOf] == "NA":
      continue
    G.add_vertex(each)
    additional_query.append("Gene")
    all_names.append(each)
    val_length_of_val = len(each) #* 15
    length_of.append(val_length_of_val)
    val_breadth_of_val = 30
    for each in each_category:
      if each in category_values:
        category_values[each].append(uniprot_list[each][indexOf])
      else:
        category_values.update({each:[uniprot_list[each][indexOf]]})
    breadth_of.append(val_breadth_of_val)
    category_true.append(uniprot_list['category_true'][indexOf])
    
  for each in each_category:
    G.add_vertex(each)
    all_names.append(each)
    val_length_of = len(each)
    length_of.append(val_length_of)
    if val_length_of > 18:
      val_breadth_of_val = val_length_of/18.0*50.0
    else:
      val_breadth_of_val = 50
    breadth_of.append(val_breadth_of_val)
   
    additional_query.append("Function")
    additional_each.append(0.0)
    category_true.append(0.0)
  
  for each_node in merged_vertex:
    for each_cat in each_category:
      index = uniprot_list['name'].index(each_node)
      if uniprot_list[each_cat][index] == 1.0:
        G.add_edge(each_node,each_cat,name=each_node + " " + each_cat,interaction="with")
  G.vs
  G.vs["name"] = all_names
  
  
  for each in each_category:
    G.vs[each] = category_values.get(each,[]) + additional_each
  
  G.vs["pine_query"] = additional_query 
  G.vs["pine_length"] = length_of
  G.vs["pine_breadth"] = breadth_of
  
  degree = G.degree()
  G.vs["pine_degree"] = degree
  G.vs['pine_category_true'] = category_true
  
  cy = CyRestClient(port=CYREST_PORT)
  g_cy = cy.network.create_from_igraph(G, name="Category Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Category-Network-Style')
  basic_settings = {
    #'NODE_FILL_COLOR':"#33CCFF",
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_WIDTH':"1",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  
  shape_kv_pair = {
    "Function":"OCTAGON",
    "Gene":"ROUND_RECTANGLE",
  }
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)
  height_kv_pair = {
    "Function":"70",
    "Gene":"30",
  }
  width_kv_pair = {
    "Function":"210",
    "Gene":"60",
  }
  label_color_kv_pair = {
    "Function":"#000000",
    "Gene":"#000000",
  }
  node_label_size =  [
    {
        "value": "2",
        "lesser": "1",
        "equal": "32",
        "greater": "32"
    },
    {
        "value": "100",
        "lesser": "8",
        "equal": "8",
        "greater": "1"
    }
  ] 
  my_style.create_continuous_mapping(column='pine_length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_HEIGHT', mappings=height_kv_pair)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)

  my_style = get_category(my_style,category_true,"2",additional_query,all_names,each_category)
  
  cy.style.apply(my_style, g_cy)
  data = [
    {
    "visualPropertyDependency": "nodeSizeLocked",
    "enabled": False
    }
  ]
  response = request_retry(f"{CYREST_URL}/v1/styles/Category-Network-Style/dependencies", 'PUT', json=data)
  
def singleFC(my_style, uniprot_list, type):
  ''' Styling specific to singleFC case '''
  basic_settings = {
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'NODE_FILL_COLOR':"#E2E2E2",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  if min(uniprot_list['FC1']) != 0 and max(uniprot_list['FC1']) != 0:
    points1 =  [
      {
        "value": min(uniprot_list['FC1']),
        "lesser": "#FFFF99" ,
        "equal": "#3399FF", 
        "greater": "#3399FF"
      },
      {
        "value": 0,
        "lesser": "#FFFFFF",
        "equal": "#FFFFFF",
        "greater": "#FFFFFF"
      },
      {
        "value": max(uniprot_list['FC1']),
        "lesser": "#FF0000",
        "equal": "#FF0000", 
        "greater": "#E2E2E2" 
      }
    ]
  elif min(uniprot_list['FC1']) == 0:
    points1 =  [
      {
        "value": 0,
        "lesser": "#FFFF99",
        "equal": "#FFFFFF",
        "greater": "#FFFFFF"
      },
      {
        "value": max(uniprot_list['FC1']),
        "lesser": "#3399FF",
        "equal": "#3399FF",
        "greater": "#E2E2E2"
      }
    ]
  else:
    points1 =  [
      {
        "value": min(uniprot_list['FC1']),
        "lesser": "#FFFF99",
        "equal": "#FF0000",
        "greater": "#FF0000"
      },
      {
        "value": 0,
        "lesser": "#FFFFFF",
        "equal": "#FFFFFF",
        "greater": "#E2E2E2"
      }
    ]
  my_style.create_continuous_mapping(column='FC1',vp='NODE_FILL_COLOR',col_type='Double',points=points1)
  
  if type == "5":
    data = {
      "bypass": "true",
      "nodeList": "pine_query:EI",
      "propertyList": "Fill Color",
      "valueList": "#E2E2E2"
    }
    
  else:
    data = {
      "bypass": "true",
      "nodeList": "pine_FC_exists:0.0",
      "propertyList": "Fill Color",
      "valueList": "#E2E2E2"
    }

  response = request_retry(f"{CYREST_URL}/v1/commands/node/set properties", 'POST', json=data)
  
  data = {
    "bypass": "true",
    "nodeList": "pine_regulated:Up",
    "propertyList": "Fill Color, Label Color",
    "valueList": "#FF6633, #FFFFFF"
  }
  response = request_retry(f"{CYREST_URL}/v1/commands/node/set properties", 'POST', json=data)
    
  data = {
    "bypass": "true",
    "nodeList": "pine_regulated:Down",
    "propertyList": "Fill Color, Label Color",
    "valueList": "#0000CC, #FFFFFF"
  }
  response = request_retry(f"{CYREST_URL}/v1/commands/node/set properties", 'POST', json=data)
  
  data = {
    "bypass": "true",
    "nodeList": "pine_regulated:None",
    "propertyList": "Fill Color",
    "valueList": "#E2E2E2"
  }
  response = request_retry(f"{CYREST_URL}/v1/commands/node/set properties", 'POST', json=data)
  
  return(my_style)
  
def multipleFC(my_style,FC_exists,query,func,name,max_FC_len,uniprot_list, unique_labels):
  ''' Styling specific to multipleFC case '''
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 
  bar_columns = ""
  color_columns = ""
  min_fc = 0
  max_fc = 0
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    term_val_FC = unique_labels[i-1]
    if i == 1:
      min_fc = min(uniprot_list[term_FC])
      max_fc = max(uniprot_list[term_FC])
    else:
      if min(uniprot_list[term_FC]) < min_fc:
        min_fc = min(uniprot_list[term_FC])
      if max(uniprot_list[term_FC]) > max_fc:
        max_fc = max(uniprot_list[term_FC])
    bar_columns = bar_columns + '\"' + term_val_FC + '\"' + ","
    color_columns = color_columns + '\"' + color_code[i-1] + '\"' + ","
  bar_columns = bar_columns[:-1]
  color_columns = color_columns[:-1]
 
  value = "org.cytoscape.BarChart:{" +'\"' + "cy_range" + '\"' + ":[" + str(min_fc) + "," + str(max_fc) + "]," + '\"' + "cy_globalRange" + '\"' + ":true," + '\"' + "cy_colors" + '\"' + ":[" + color_columns + "]," + '\"' + "cy_dataColumns" + '\"' + ":[" + bar_columns + "]," + '\"' + "cy_domainLabelPosition" + '\"' + ":" + '\"' + "UP_90" + '\",' + '\"' + "cy_borderWidth" + '\"' + ':3.0,' + '\"' + "cy_separation" + '\"' + ':0.3,' + '\"' + "cy_axisWidth" + '\"' + ':5' +  "}"
  kv_pair = {
    "1":value
  }
  kv_pair_node_position = {
    "1":"S,N,c,0.00,0.00"
  }
  my_style.create_discrete_mapping(column='pine_FC_exists', col_type='Double', vp='NODE_CUSTOMGRAPHICS_1',mappings=kv_pair)
  my_style.create_discrete_mapping(column='pine_FC_exists', col_type='Double', vp='NODE_LABEL_POSITION',mappings=kv_pair_node_position)

  kv_pair_color = {}
  width_kv_pair = {}
  if func == "1":
    for i in query:
      if i == "NA":
        color_value = "#E2E2E2"
      else:
        color_value = "#FFFFFF"
      if i not in kv_pair_color:
        kv_pair_color.update({i:color_value})
    my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_FILL_COLOR',mappings=kv_pair_color)
    
  else:
    for i,does_fc_exist,each_name in zip(query,FC_exists,name):
      if i == "Function":
        color_value = "#E2E2E2"
        width_value = "3"
      elif i == "Gene" and does_fc_exist == 1.0:
        color_value = "#FFFFFF"
        width_value = "3"
      elif i == "Gene" and does_fc_exist != 1.0:
        color_value = "#FFFF99"
        width_value = "3"
      elif i == "Site":
        color_value = "#FFFFFF"
        width_value = "3"
      else:
        color_value = "#E2E2E2"
        width_value = "3"
        
      if each_name not in kv_pair_color:
        kv_pair_color.update({each_name:color_value})
        width_kv_pair.update({each_name:width_value})
        
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_FILL_COLOR',mappings=kv_pair_color)
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_BORDER_WIDTH', mappings=width_kv_pair)
  return(my_style)
  
def get_category(my_style,is_category_present,cat_val,query,name,each_category):
  ''' Styling specific to category case '''
  
  basic_settings = {
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  my_style.update_defaults(basic_settings)
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 

  bar_columns = ""
  color_columns = ""
  i = 0
  for each in each_category:
    i+=1
    bar_columns = bar_columns + '\"' + each + '\"' + ","
    color_columns = color_columns + '\"' + color_code[i-1] + '\"' + ","
  bar_columns = bar_columns[:-1]
  color_columns = color_columns[:-1]
  
  kv_pair = {}
  value = "org.cytoscape.PieChart:{" + '\"' + "cy_colors" + '\"' + ":[" + color_columns + "]," + '\"' + "cy_dataColumns" + '\"' + ":[" + bar_columns + "]}"
 
  kv_pair = {
    "1.0":value
  }
  kv_pair_node_position = {
    "1.0":"S,N,c,0.00,0.00",
    "0.0":"c,c,c,0.00,0.00"
  }
  kv_edge_width = {
    "1.0":"1",
    "0.0":"1"
  }
  kv_node_border_width = {
    "1.0":"0",
    "0.0":"3"
  }
  my_style.create_discrete_mapping(column='pine_category_true', col_type='Double', vp='NODE_CUSTOMGRAPHICS_1',mappings=kv_pair)
  my_style.create_discrete_mapping(column='pine_category_true', col_type='Double', vp='NODE_LABEL_POSITION',mappings=kv_pair_node_position)
  my_style.create_discrete_mapping(column='pine_category_true', col_type='Double', vp='EDGE_WIDTH',mappings=kv_edge_width)
  my_style.create_discrete_mapping(column='pine_category_true', col_type='Double', vp='NODE_BORDER_WIDTH', mappings=kv_node_border_width)
  
  if cat_val == "2":
    kv_pair_color = {}
    width_kv_pair = {}
    node_width = {}
    for i,each_is_category_present,each_name in zip(query,is_category_present,name):
        if i == "Function":
          color_value = "#E2E2E2"
          width_value = "3"
          node_width_value = "210"
        elif i == "Gene" and each_is_category_present == 1.0:
          color_value = "#FFFFFF"
          width_value = "0"
          node_width_value = "30"
        else:
          color_value = "#E2E2E2"
          width_value = "3"
          node_width_value = "60"
        if each_name not in kv_pair_color:
          kv_pair_color.update({each_name:color_value})
          width_kv_pair.update({each_name:width_value})
          node_width.update({each_name:node_width_value})

    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_WIDTH', mappings=node_width)      
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_FILL_COLOR',mappings=kv_pair_color)
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_BORDER_WIDTH', mappings=width_kv_pair)
  
  return(my_style)
  
def pval(my_style):
  ''' If pval significant, border node with blue color, width = 3 '''
  width_kv_pair = {
    "1":"3",
  }
  bc_kv_pair = {
    "1":"#0000FF",
  }
  my_style.create_discrete_mapping(column='pine_outline', col_type='Double', vp='NODE_BORDER_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='pine_outline', col_type='Double', vp='NODE_BORDER_PAINT', mappings=bc_kv_pair)
  return(my_style)
  
def snp_bold(my_style):
  ''' Not used. If node crepresents SNP, make label bold '''
  ff_kv_pair = {
    "1.0":"SansSerif.bold,bold,12",
    "0.0":"SansSerif.plain,plain,12",
  }
  my_style.create_discrete_mapping(column='SNP', col_type='Double', vp='NODE_LABEL_FONT_FACE', mappings=ff_kv_pair)
  return(my_style)
  
def color(my_style, uniprot_list):
  ''' Styling specific to list only (noFC) case '''
  color_kv_pair = {}
  for each_query in uniprot_list["query"]:
    if each_query == "NA":
      if each_query not in color_kv_pair:
        color_kv_pair.update({each_query:"#E2E2E2"})
    else:
      color_kv_pair.update({each_query:"#FFFFBF"})
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  return(my_style)

def cy_pathways_style(cluster, each_category, max_FC_len, pval_style, uniprot_list, type, all_prot_site_snps, uniprot_query, unique_labels,mult_mods_of_int):
  ''' Based on top annotations picked, construct function interaction network + visualization and styling '''
  G = igraph.Graph()
  color_code = ["#FF9933", "#00FFFF", "#00FF00", "#FF66FF", "#FFFF66", "#9999FF"] 
  
  cluster_list = cluster
  
  name = []
  function_only = cluster_list.keys()
  query = []
  length_of = []
  breadth_of = []
   
  merged_vertex = []
  merged_vertex_sites_only = []
  merged_vertex_unimod_only = []
  query_val_noFC = []
  count_each = 0
  all_interactions = []
  function_fc_val = {}
  is_snp = []
  up_or_down = []
  up_or_down_to_append = []
  calc_up = 0
  calc_down = 0
  total_genes = 0
  category_present = 0
  gene_list_per_interaction = []
  for each in each_category:
    category_present = 1
    break
    
  for each in cluster_list:
    if each not in merged_vertex:
      G.add_vertex(each)
      merged_vertex.append(each)
      merged_vertex_sites_only.append(each)
      merged_vertex_unimod_only.append(each)
      if each in function_only:
        calc_up = 0
        calc_down = 0
        total_genes = 0
        up_or_down_to_append = []
        gene_list_per_interaction = []
        is_snp.append(0.0)
        query.append('Function')
        if max_FC_len == 0 and not category_present:
          query_val_noFC.append('Function')
        val_length_of = len(each)
       
        length_of.append(val_length_of)
        if val_length_of > 18:
          val_breadth_of_val = val_length_of/18.0*50.0
        else:
          val_breadth_of_val = 50
        breadth_of.append(val_breadth_of_val)
        
    for each_gene in cluster_list[each]:
      if each_gene not in merged_vertex:
        G.add_vertex(each_gene)
        merged_vertex.append(each_gene)
        if type == "1" or type == "5":
          up_or_down_to_append.append("NA")
        if each_gene not in function_only:
          is_snp.append(0.0)
          indexOf = uniprot_list['name'].index(each_gene.lower())
          if uniprot_list['query'][indexOf] != "NA":
            query.append('Gene')
          else:
            query.append("EI")
          merged_vertex_sites_only.append(uniprot_list['ambigious_genes'][indexOf])
          merged_vertex_unimod_only.append(uniprot_list['ambigious_genes'][indexOf])
          if max_FC_len == 0 and not category_present:
            index_noFC = uniprot_list["name"].index(each_gene.lower())
            query_val_noFC.append((uniprot_list["query"])[index_noFC])
          val_length_of_val = len(each_gene) #* 15
          length_of.append(val_length_of_val)
          val_breadth_of_val = 30
          breadth_of.append(val_breadth_of_val)
                        
          if type == "5" or type == "6":
            indices = [i for i, x in enumerate(uniprot_list['name']) if x == each_gene.lower()]
            for each_index in indices:
              if uniprot_list['site'][each_index] != "NA":
                G.add_vertex(uniprot_list['site'][each_index])
                merged_vertex.append(uniprot_list['site'][each_index])
                merged_vertex_sites_only.append(uniprot_list['ambigious_site'][each_index])
                combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
                without_unimod = re.sub(combined_pat, '', uniprot_list['ambigious_site'][each_index])
                merged_vertex_unimod_only.append(without_unimod)
                #get corresponding prot id and check for snps
                get_corres_prot = get_query_from_list(uniprot_query, [each_gene])
                if get_corres_prot:
                  if list(get_corres_prot.keys())[0] in all_prot_site_snps:
                    if each_site in all_prot_site_snps[get_corres_prot]['Position']:
                      is_snp.append(1.0)
                    else:
                      is_snp.append(0.0)
                  else:
                    is_snp.append(0.0)
                else:
                  is_snp.append(0.0)                  
                query.append('Site')
                val_length_of_val = len(each_gene) #* 15
                length_of.append(val_length_of_val)
                val_breadth_of_val = 30
                breadth_of.append(val_breadth_of_val)
                name_edge = uniprot_list['site'][each_index] + " with " + each_gene
                G.add_edge(uniprot_list['site'][each_index],each_gene,name=name_edge)
                all_interactions.append(name_edge)
                if type == "5":
                  up_or_down_to_append.append("NA")
                  if uniprot_list['site'][each_index] not in gene_list_per_interaction:
                    gene_list_per_interaction.append(uniprot_list['site'][each_index])
                    total_genes += 1                  
                    FC_val_each_gene = (uniprot_list['FC1'])[each_index]
                    if FC_val_each_gene > 0:
                      calc_up += 1
                    elif FC_val_each_gene < 0:
                      calc_down -= 1  
      else:
        if type == "5":
          indices = [i for i, x in enumerate(uniprot_list['name']) if x == each_gene.lower()]
          gene_list_per_interaction = []
          for each_index in indices:
            if uniprot_list['site'][each_index] not in gene_list_per_interaction:
              gene_list_per_interaction.append(uniprot_list['site'][each_index])
              total_genes += 1                  
              FC_val_each_gene = (uniprot_list['FC1'])[each_index]
              if FC_val_each_gene > 0:
                calc_up += 1
              elif FC_val_each_gene < 0:
                calc_down -= 1
            
      name_edge = each + " with " + each_gene
      G.add_edge(each,each_gene,name=name_edge)
      all_interactions.append(name_edge)
      if type == "1":
        if each_gene not in gene_list_per_interaction:
          gene_list_per_interaction.append(each_gene)
          total_genes += 1
          indexOf = uniprot_list["name"].index(each_gene.lower())
          FC_val_each_gene = (uniprot_list['FC1'])[indexOf]
          if FC_val_each_gene > 0:
            calc_up += 1
          elif FC_val_each_gene < 0:
            calc_down -= 1
    if type == "1" or type == "5":            
      if (calc_up/total_genes*100) > 60: 
        up_or_down.append("Up")
      elif (abs(calc_down)/total_genes*100) > 60:
        up_or_down.append("Down")
      else: 
        up_or_down.append("None")
      up_or_down.extend(up_or_down_to_append)
  
  G.vs
  G.vs["name"] = merged_vertex
  G.vs["shared name"] = merged_vertex_sites_only
  G.vs["substitute name"] = merged_vertex_unimod_only
  is_category_present = []
  first_cat_iteration = 0
  for each in each_category:
    category_present = 1
    category = []
    k = 0
    for each_vertex_name in merged_vertex:
      if each_vertex_name.lower() in uniprot_list["name"]:
        index = uniprot_list["name"].index(each_vertex_name.lower())
        category.append((uniprot_list[each])[index])
        
        if first_cat_iteration == 0 :
          if (uniprot_list["category_true"])[index] == 1.0:
            is_category_present.append(1.0)
          else:
            is_category_present.append(0.0)
        else:
          if (uniprot_list["category_true"])[index] == 1.0:
            is_category_present[k] = 1.0
        
      else:
        category.append((''))
        if first_cat_iteration == 0:
          is_category_present.append(0.0)
      k+=1    
    first_cat_iteration +=1
    G.vs[each] = category
      
  add_term_pval = []
  FC_exists = []
  significant_val = []
  fc_exists_count = 0
  add_term_FC = []
  value_labels = []
  if type == "5" or type == "6":
    search_list_for_fc = [a.lower() for a in uniprot_list["site"]]
  else:
    search_list_for_fc = uniprot_list["name"]
  domain_labels = []
  for i in range(1,max_FC_len+1):
    k = 0
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    if unique_labels:
      val_term_FC = unique_labels[i-1]
      val_term_pval = unique_labels[i-1] + " pval"
    else:
      val_term_FC = term_FC
      val_term_pval = term_pval
    domain_labels.append(str(val_term_FC))
    if add_term_FC:
      if not value_labels:
        value_labels = [ [str(add_term_FC[i])] for i in range(0, len(add_term_FC)) ]
      else:
        for i in range(0, len(add_term_FC)):
          value_labels[i].append(str(add_term_FC[i]))       
    add_term_FC = []
    add_term_pval = []
    significant_val = []
    for each_vertex_name in merged_vertex:
      if each_vertex_name.lower() in search_list_for_fc:
        index = search_list_for_fc.index(each_vertex_name.lower())
        add_term_FC.append((uniprot_list[term_FC])[index])
        add_term_pval.append((uniprot_list[term_pval])[index])
        significant_val.append((uniprot_list["significant"])[index])
        if fc_exists_count == 0:
          if ((uniprot_list[term_FC])[index]) != 0:
            FC_exists.append(1.0)
          elif each_vertex_name in function_only:
            FC_exists.append(-1.0)
          else:
            FC_exists.append(0.0)
        else:
          if ((uniprot_list[term_FC])[index]) != 0:
            FC_exists[k] = 1.0
              
      else:
        if max_FC_len == 1 and not (type == "5" or type == "6"):
          add_term_FC.append(100.0)
       
        else:
          if each_vertex_name in function_only:
            add_term_FC.append(100.0)
          else:
            add_term_FC.append(-100.0)
          
        add_term_pval.append(1.0)
        significant_val.append(0.0)
        
        if fc_exists_count == 0 and each_vertex_name in function_only:
          FC_exists.append(-1.0)
        else:
          FC_exists.append(0.0)   
             
      k+=1
    fc_exists_count+=1

    G.vs[val_term_FC] = add_term_FC
    G.vs[val_term_pval] = add_term_pval
    G.vs["pine_outline"] = significant_val
    G.vs["pine_FC_exists"] = FC_exists
    
  if add_term_FC:
    if not value_labels:
      value_labels = [ [str(add_term_FC[i])] for i in range(0, len(add_term_FC)) ]
    else:
      for i in range(0, len(add_term_FC)):
        value_labels[i].append(str(add_term_FC[i]))  
        
  if max_FC_len > 1:
    tot = len(query)
    G.vs["pine_domain_labels"] = [domain_labels]*tot
    G.vs["pine_value_labels"] = value_labels
    
  G.vs["pine_query"] = query
  G.vs["pine_length"] = length_of
  G.vs["pine_breadth"] = breadth_of
  degree = G.degree()
  G.vs["pine_degree"] = degree
  if up_or_down:
    G.vs["pine_regulated"] = up_or_down
  if category_present == 1:
    G.vs["pine_category_true"] = is_category_present  
    
  if max_FC_len == 0 and not category_present:
    G.vs["pine_query_val"] = query_val_noFC
  
  cy = CyRestClient(port=CYREST_PORT)

  g_cy = cy.network.create_from_igraph(G, name="Ontology Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  
  my_style = cy.style.create('GAL_Style3')
  
  basic_settings = {
    'EDGE_TRANSPARENCY':"150",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_WIDTH':"1",
    'EDGE_STROKE_UNSELECTED_PAINT':"#000000"
  }
  my_style.update_defaults(basic_settings)
  if not mult_mods_of_int:
    my_style.create_passthrough_mapping(column='substitute name', vp='NODE_LABEL', col_type='String')
  else:
    my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  if type == "5" or type == "6":
    shape_kv_pair = {
      "Function":"OCTAGON",
      "Gene":"ROUND_RECTANGLE",
      "EI":"ROUND_RECTANGLE",
      "Site":"ELLIPSE"
    }
  else:
    shape_kv_pair = {
      "Function":"OCTAGON",
      "Gene":"ROUND_RECTANGLE",
      "EI":"ROUND_RECTANGLE"
    }
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)
  if type == "2" or type == "6":
    height_kv_pair = {
      "Function":"130",
      "Gene":"90",
      "EI":"90",
      "Site":"130"
    } 
  else:
    height_kv_pair = {
      "Function":"70",
      "Gene":"30",
      "EI":"30",
      "Site":"30"
    }
  if type == "5":
    width_kv_pair = {
      "Function":"210",
      "Gene":"70",
      "Site":"60",
      "EI":"70"
    }
  elif type == "2" or type == "6":
    width_kv_pair = {
      "Function":"250",
      "Gene":"120",
      "Site":"90",
      "EI":"120"
    }
  else:
    width_kv_pair = {
      "Function":"210",
      "Gene":"60",
      "Site":"60"
    }

  label_color_kv_pair = {
      "Function":"#000000",
      "Gene":"#000000",
      "Site":"#000000"
  }
  if type == "2" or type == "6":
    node_label_size =  [
      {
        "value": "2",
        "lesser": "1",
        "equal": "32",
        "greater": "32"
      },
      {
        "value": "100",
        "lesser": "8",
        "equal": "8",
        "greater": "1"
      }
    ] 
  else:
    node_label_size =  [
      {
        "value": "2",
        "lesser": "1",
        "equal": "14",
        "greater": "14"
      },
      {
        "value": "100",
        "lesser": "8",
        "equal": "8",
        "greater": "1"
      }
    ]
  my_style.create_continuous_mapping(column='pine_length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_HEIGHT', mappings=height_kv_pair)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)

  pval_sig = 0
  if max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list, type)
    if pval_style:
      pval_sig = 1
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,FC_exists,query,"2",merged_vertex, max_FC_len, uniprot_list, unique_labels)

    if pval_style:
      pval_sig = 1
    
  if category_present:
    my_style = get_category(my_style,is_category_present,"2",query,merged_vertex,each_category)
  
  if max_FC_len == 0 and not category_present:
    basic_settings = {
      'NODE_BORDER_WIDTH':"3"
    }
    my_style.update_defaults(basic_settings)
    color_kv_pair = {}
    label_color_kv_pair = { }
    for each_val in query_val_noFC:
      if each_val == "NA":
        color_kv_pair.update({each_val:"#E2E2E2"})
        label_color_kv_pair.update({each_val:"#000000"})
      elif each_val == "Function":
        color_kv_pair.update({each_val:"#333333"})
        label_color_kv_pair.update({each_val:"#FFFFFF"})
      else:
        color_kv_pair.update({each_val:"#FFFFBF"})
        label_color_kv_pair.update({each_val:"#000000"})
    
    my_style.create_discrete_mapping(column='pine_query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)
    my_style.create_discrete_mapping(column='pine_query_val', col_type='String', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  
  if pval_sig == 1:
    my_style = pval(my_style)
 
  cy.style.apply(my_style, g_cy)
  
  data = [
    {
    "visualPropertyDependency": "nodeSizeLocked",
    "enabled": False
    }
  ]
  response = request_retry(f"{CYREST_URL}/v1/styles/GAL_Style3/dependencies", 'PUT', json=data)
  
def remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, settings_file):
  ''' In case of errors, remove all output files and close Cytoscape '''
  try:
    if cy_debug:
      logging.handlers = []
      if path.exists(logging_file):
        os.remove(logging_file)
    if cy_session and path.exists(cy_session):
      os.remove(cy_session)
    if cy_out and path.exists(cy_out):
      os.remove(cy_out)
    if cy_cluego_out and path.exists(cy_cluego_out):
      os.remove(cy_cluego_out)
    if settings_file and path.exists(settings_file):
      os.remove(settings_file)
    if path_to_new_dir and path.exists(path_to_new_dir):
      os.rmdir(path_to_new_dir)
  except:
    eprint(f"Error: Failed to delete results directory at: {path_to_new_dir}")

  #Close cytoscape
  try:
    requests.get(f"{CYREST_URL}/v1/commands/command/quit")
  except:
    pass
    
def main(argv):
  '''
  Run PINE - a tool for visualizing protein-protein interactions and annotations
  Parse command line arguments
  Create an analysis session directory or load an existing session
  Start Cytoscape
  Parse input
  Search protein-protein interaction databases (STRING and Genemania)
  Obtain annotation terms for proteins using Cytoscape App, ClueGO
  Create interaction and annotation networks in Cytoscape
  Save Cytoscape and network and analysis results
  '''
  cy_in = ""
  cy_species = ""
  cy_lim = 0
  cy_score = 0.4
  cy_debug = True
  cy_map = ""
  merged_out_dict = {}
  help = False
  cy_type = ""
  cy_pval = 0
  cy_cluego_grouping = "medium"
  cy_fc_cutoff = 0
  cy_pval_cutoff = 1
  select_terms = "pathways"
  leading_term_selection = "no. of genes per term"
  cluego_reference_file = ""
  cy_cluego_inp_file = ""
  cluego_pval = 0.05
  cy_run = "both"
  cy_ptm_sites = False
  cy_fasta_file = ""
  cy_mods = ""
  leading_term_cluster = "" 
  cy_enzyme = ""
  gui_mode = False
  cy_out_dir = None
  cy_out_name = None
  cy_exe = None
  cy_ambi = False
  count_args = 0
  
  try:
    opts, args = getopt.getopt(argv, "i:s:l:t:r:m:o:j:ng:f:p:z:h:a:y:u:d:b:x:e:c:k",["in=","species=","limit=","type=","score=","mapping=","output=","significant","grouping=","fccutoff=","pvalcutoff=","visualize=","reference-path=","input_cluego=","cluego-pval=","run=","mods=","fasta-file=","enzyme=","gui","cytoscape-executable=","cytoscape-session-file=","exclude-ambiguity", "output-name="])
    for opt, arg in opts:
      count_args += 1
      if opt in ("-i","--in"):
        cy_in = arg
      elif opt in ("-s","--species"):
        cy_species = arg
      elif opt in ("-l","--limit"):
        cy_lim = arg
      elif opt in ("-t","--type"):
        cy_type = arg
      elif opt in ("-r","--score"):
        cy_score = arg
      elif opt in ("-n","--significant"):
        cy_pval = 1
      elif opt in ("-u","--run"):
        cy_run = arg
      elif opt in ("-f","--fccutoff"):
        cy_fc_cutoff = arg
      elif opt in ("-p","--pvalcutoff"):
        cy_pval_cutoff = arg
      elif opt in ("-m","--mapping"):
        cy_map = arg
      elif opt in ("-o","--output"):
        cy_out_dir = arg
      elif opt in ("-d", "--mods"):
        cy_mods = arg
      elif opt in ("-b", "--fasta-file"):
        cy_fasta_file = arg
      elif opt in ("-a","--input_cluego"):
        cy_cluego_inp_file = arg
      elif opt in ("-z","--visualize"):
        select_terms = arg
      elif opt in ("-y","--cluego-pval"):
        cluego_pval = arg
      elif opt in ("-h","--reference-path"):
        cluego_reference_file = arg
      elif opt in ("-g","--grouping"):
        cy_cluego_grouping = arg
      elif opt in ("-x","--enzyme"):
        cy_enzyme = arg
      elif opt in ("--gui",):
        gui_mode = True
      elif opt in ("-e","--cytoscape-executable"):
        cy_exe = arg
      elif opt in ("-k","--exclude-ambiguity"):
        cy_ambi = True
      elif opt in ("-j", "--output-name"):
        cy_out_name = arg
      else:
        help = True      
  except getopt.GetoptError as e:
    help = True
  
  if not count_args:
    help = True
    
  if help:
    print("PINE")
    print("---------------------------------------------------------------------------------------------")
    print("Usage:         pine.py -i input.csv -o output_dir -t input_type -s species -m cluego_map_file.gz")
    print("Argument:      -i [--in]: input file in csv format with the following headers as applicable: ProteinID, FC, pval, adj.pval, Label, Category, Peptide")
    print("Argument:      -o [--output]: path to output directory")     
    print("Argument:      -t [--type]: analysis type [Allowed: noFC, singleFC, multiFC, category, singlefc-ptm, multifc-ptm]")   
    print("Argument:      -s [--species]: species [Allowed: human, mouse, rat]")   
    print("Argument:      -x [--enzyme]: enzyme name [Allowed: Trypsin, Trypsin_p, Lys_n, Asp_n, Arg_c, Chymotrypsin, Lys_c]")   
    print("Argument:      -d [--mods]: comma separated list of modifications of interest [Example: S,T,Y or K(Unimod:1) or S[+80]]")
    print("Argument:      -b [--fastafile]: path to fasta file")
    print("Argument:      -m [--mapping]: path to cluego mapping file compressed in .gz format") 
    print("Argument:      -e [--cytoscape-executable]: the path to the Cytoscape executable")
    print("Argument(opt): -f [--fccutoff]: fold change cutoff for input [Default: abs(FC) >= 0.0]")
    print("Argument(opt): -p [--pvalcutoff]: pvalue cutoff for input [Default: pval > 1.0]")
    print("Argument(opt): -n [--significant]: outline statistically significant nodes, i.e pval>0.0")
    print("Argument(opt): -k [--exclude-ambiguity]: exclude ambigious genes and sites")  
    print("Argument(opt): -u [--run]: interaction databases [Allowed: string, genemania, both; Default: both]")
    print("Argument(opt): -r [--score]: interaction confidence score for string [Default:0.4, Range 0-1]")
    print("Argument(opt): -l [--limit]: maximum number of external interactors [Default:0, Range:0-100]")
    print("Argument(opt): -z [--visualize]: ontology type [Allowed: biological process, cellular component, molecular function, pathways, all; Default: pathways].  Pathways include REACTOME, KEGG, CLINVAR, CORUM and Wiki.")
    print("Argument(opt): -g [--grouping]: network specificity indicating general, representative and specific pathways [Allowed: global, medium, detailed; Default: medium]")
    print("Argument(opt): -y [--cluegopval]: pvalue cutoff for enrichment analysis [Default: 0.05]")
    print("Argument(opt): -h [--referencepath]: path to background reference file for enrichment")
    print("Argument(opt): -a [--inputcluego]: filtered cluego file with ontology terms of interest")
    sys.exit()
  
  if not cy_in or not cy_species or not cy_type or not cy_out_dir or not cy_exe or not cy_map:
    eprint("Error: Mandatory parameter not provided. Please provide path to cytoscape exe, path to ClueGO mapping file, input csv file, species, run type and output directory")
    sys.exit(1)

  if not os.path.isdir(cy_out_dir):
    eprint("Error: output is not a directory")
    sys.exit(1)
    
  if not path.exists(cy_in):
    eprint("Error: Path to Input file " + cy_in + " does not exist")
    sys.exit(1)
    
  if cy_fasta_file and not path.exists(cy_fasta_file):
    eprint("Error: Path to FASTA file " + cy_fasta_file + " does not exist")
    sys.exit(1)
    
  if cy_cluego_inp_file and not path.exists(cy_cluego_inp_file):
    eprint("Error: Path to ClueGO Input file " + cy_cluego_inp_file + " does not exist")
    sys.exit(1)
    
  if cluego_reference_file:
    if not path.exists(cluego_reference_file):
      eprint("Error: Path to ClueGO Reference file " + cluego_reference_file + " does not exist")
      sys.exit(1)

    try:
      with open(cluego_reference_file) as f:
        for line in f:
          pass # just check if file can be read, not a binary file
    except FileNotFoundError:
      eprint("Error: Reference file is missing")
      sys.exit(1)
    except UnicodeDecodeError:
      eprint("Error: Reference file must be a text file")
      sys.exit(1)
    
  timestamp = datetime.utcnow().replace(tzinfo=dt.timezone.utc).astimezone().replace(microsecond=0).isoformat()
  hr_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  
  path_to_new_dir = ""
  # create file names
  if cy_cluego_inp_file:
    cy_out = None
    cy_cluego_out = None
    path_to_cluego = (os.path.abspath(cy_cluego_inp_file))
    if path_to_cluego.endswith(".cluego.txt"):
      session_name = ".".join(path_to_cluego.split(".")[:-2])
    else:
      session_name = path_to_cluego
    cy_session = session_name + ".cys"
    logging_file = session_name + ".log"
    reanalyze_flag = True
    cy_settings_file = None
  else:  
    if cy_out_name is not None:
      cy_out_name = re.sub(r"[\/\\:\*\?\"<>\|]", "", cy_out_name) # remove any characters that are not alphanumeric underscore or dash
    if cy_out_name:
      path_to_new_dir = os.path.join(os.path.abspath(cy_out_dir), cy_out_name)
    else:
      path_to_new_dir = os.path.join(os.path.abspath(cy_out_dir), hr_timestamp) + "_PINE"
    if os.path.exists(path_to_new_dir):
      eprint(f"Error: Output directory already exists: {path_to_new_dir}.  Please select a new name")
      sys.exit(1)
    os.mkdir(path_to_new_dir) 
    cy_out = os.path.join(path_to_new_dir, "Interactions.csv")
    cy_cluego_out = os.path.join(path_to_new_dir, "PINE.cluego.txt")
    cy_cluego_log_out = os.path.join(path_to_new_dir, "PINE-cluego-log.txt")
    cy_session = os.path.join(path_to_new_dir, "PINE.cys")
    cy_settings_file = os.path.join(path_to_new_dir, "timestamp.json")
    logging_file = os.path.join(path_to_new_dir, "PINE.log")
    reanalyze_flag = False
    with open(cy_settings_file, "w") as f:
      json.dump({"timestamp": timestamp}, f)
    if gui_mode: # notify gui of the output directory
      print(f"COMMAND FILE-SESSION {path_to_new_dir}")
      
  if os.path.exists(cy_session):
    eprint("Session " + cy_session + " already exists. Please provide a different name")
    sys.exit(1)
  
  # set up logging
  if gui_mode:
    logging = setup_logger("PINE.log", logging_file, with_stdout=True)
  else:
    logging = setup_logger("PINE.log", logging_file)

  if cy_species.lower() == "human":
    tax_id = "9606"
    organism_name = "Homo Sapiens"
  elif cy_species.lower() == "mouse":
    tax_id = "10090"
    organism_name = "Mus Musculus"
  elif cy_species.lower() == "rat":
    tax_id = "10116"
    organism_name = "Rattus norvegicus"
  elif cy_species.lower() == "ecoli":
    tax_id = "199310"
    organism_name = "Escherichia coli"
  elif cy_species.lower() == "yeast":
    tax_id = "4932"
    organism_name = "Saccharomyces cerevisiae"
  elif cy_species.lower() == "arabidopsis thaliana":
    tax_id = "3702"
    organism_name = "Arabidopsis thaliana"
  elif cy_species.lower() == "roundworm":
    tax_id = "6239"
    organism_name = "Caenorhabditis elegans"
  elif cy_species.lower() == "zebrafish":
    tax_id = "7955"
    organism_name = "Danio rerio"
  elif cy_species.lower() == "fruit fly":
    tax_id = "7227"
    organism_name = "Drosophila melanogaster"
  elif cy_species.lower() == "bovine":
    tax_id = "9913"
    organism_name = "Bos taurus"
  elif cy_species.lower() == "chicken":
    tax_id = "9031"
    organism_name = "gallus gallus"
  elif cy_species.lower() == "pig":
    tax_id = "9823"
    organism_name = "sus scrofa"
  elif cy_species.lower() == "rabbit":
    tax_id = "9986"
    organism_name = "Oryctolagus cuniculus"
  elif cy_species.lower() == "sheep":
    tax_id = "9940"
    organism_name = "Ovis aries"
  elif cy_species.lower() == "dog":
    tax_id = "9612"
    organism_name = "Canis Lupus Familiaris"
  else:
    eprint("Error: Species not currently supported")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)

  if not ('\\' in cy_session or "/" in cy_session):
    cwd = os.getcwd()
    cy_session = cwd + "\\" + cy_session
    
  if not ('.cys' in cy_session):
    eprint("Error: Cytoscape session file must have a valid name with .cys extension")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  if not (float(cluego_pval) >=0.0 and float(cluego_pval) <=1.0):
    eprint("Error: Cluego pvalue range must be between 0 to 1")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  allowed_runs = ["string","genemania","both"]
  if cy_run.lower() not in allowed_runs:  
    eprint("Error: Run type must be one of the following: " + ','.join(allowed_runs))
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
    
  allowed_leading_term = ["highest significance", "no. of genes per term", "percent of genes per term", "percent genes per term vs cluster"]
  if leading_term_selection.lower() not in allowed_leading_term:
    eprint("Error: ClueGO leading term selection must be one of the following: " + (',').join(allowed_leading_term))
    sys.exit(1)
  elif leading_term_selection.lower() == "highest significance":
    leading_term_selection = "Highest%20Significance"
  elif leading_term_selection.lower() == "no. of genes per term":
    leading_term_selection = "%23Genes%20%2F%20Term"
  elif leading_term_selection.lower() == "percent of genes per term":
    leading_term_selection = "%25Genes%20%2F%20Term"
  elif leading_term_selection.lower() == "percent genes per term vs cluster":
    leading_term_selection = "%25Genes%20%2F%20Term%20vs%20Cluster"
    
  allowed_type = ["singlefc","multifc","nofc","category","singlefc-ptm","multifc-ptm"]
  # 1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = category; 5 = singlefc-ptm; 6 = multifc-ptm
  if cy_type.lower() not in allowed_type:
    eprint("Error: Input type must be one of the following: " + (',').join(allowed_type))
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  elif cy_type.lower() == "singlefc":
    cy_type_num = "1"
  elif cy_type.lower() == "multifc":
    cy_type_num = "2"
  elif cy_type.lower() == "nofc":
    cy_type_num = "3"
  elif cy_type.lower() == "category":
    cy_type_num = "4"
  elif cy_type.lower() == "singlefc-ptm":
    cy_type_num = "5"
  elif cy_type.lower() == "multifc-ptm":
    cy_type_num = "6"
    
  allowed_selections = ["biological process","cellular component","molecular function","pathways","all"]
  if select_terms.lower() not in allowed_selections:
    eprint("Error: The visualization type must be one of the following: " + (',').join(allowed_selections))
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  try:
    cy_fc_cutoff = float(cy_fc_cutoff)
  except:
    eprint("Error: FC cutoff must be a number") 
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)

  try:
    cy_pval_cutoff = float(cy_pval_cutoff)
  except:
    eprint("Error: PVal cutoff must be a number") 
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  if not (cy_pval_cutoff >= 0.0 and cy_pval_cutoff <= 1.0):
    eprint("Error: PVal cutoff must range between 0 and 1")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
    
  allowed_groups = ["global","medium","detailed"]
  if cy_cluego_grouping.lower() not in allowed_groups:
    eprint("Error: Cluego grouping must be one of the following:" + (',').join(allowed_groups))
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  try:  
    cy_score = float(cy_score)*1000
  except:
    eprint("Error: Invalid string score provided. Value must be between 0 to 1")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  if not (cy_score >= 0.0 and cy_score <= 1000.0):
    eprint("Error: Invalid string score provided. Value must be between 0 to 1; Confidence levels for string score- Low = 0.150, Medium = 0.400, High = 0.700, Highest = 0.900")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
  
  # Rounding off score - ex: 500.9 and above = 501, less = 500
  cy_score_ceil = math.ceil(cy_score)
  if((cy_score_ceil-cy_score) <= 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = cy_score_ceil
  elif ((cy_score_ceil-cy_score) > 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = math.floor(cy_score)
  
  if not (int(cy_lim) >= 0 and int(cy_lim) <=100):
    eprint("Error: Limit on additional interactors is 100. Please choose a number between 0 and 100")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    sys.exit(1)
    
  try:
    if cy_debug:
      logging.debug("Starting PINE Analysis...\n") 
      
    # get a local address port that this run will use
    found_open_port = False
    global CYREST_URL
    global CYREST_PORT
    for cyrest_port in CYREST_PORTS:
      if not port_in_use(cyrest_port):
        found_open_port = True
        break
    CYREST_PORT = cyrest_port
    CYREST_URL = "http://localhost:" + str(cyrest_port)

    if not found_open_port:
      eprint("Could not find an open port to start Cytoscape. Please close a previous PINE generated Cytoscape session and start the run again.")
      sys.exit(1)
    
    if not cy_cluego_inp_file:
      # Check if cluego path exists        
      if not path.exists(cy_map):
        eprint("Error: Path to ClueGO mapping file " + cy_map + " does not exist")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      else:
        if organism_name.lower() not in cy_map.lower():
          eprint("Error: Species mismatch.  Species parameter provided is " + organism_name + " which does not match species contained in path to ClueGO mapping file is " + cy_map)
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
          sys.exit(1)
    
      if not (".gz" in cy_map and "gene2accession" in cy_map):
        eprint("Error: ClueGO mapping file must refer to the species gene2accession .gz file")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)

    database_dict = {} 
    mods_list = []
    if cy_type_num == "5" or cy_type_num == "6":
      if not cy_fasta_file or not cy_mods or not cy_enzyme:
        eprint("Error: Fasta file, Enzyme and List of Modifications are mandatory for site analysis")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)        
      try:
        sys.stdout.write("Send"+ str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")) + "\n")
        sys.stdout.flush()
        database_dict = db_handling(cy_fasta_file)
        sys.stdout.write("Retrieved"+ str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")) + "\n")
        sys.stdout.flush()
        if len(database_dict) == 0:
          eprint("Error: Fasta file is empty or not in fasta format")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
          sys.exit(1)
      except FileNotFoundError:
        eprint("Error: Fasta file is missing")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      except UnicodeDecodeError:
        eprint("Error: Fasta file must be in fasta format")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      except PineError as e:
        eprint(f"Error: {e}")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      mods_list = cy_mods.split(",") 
      allowed_enzyme = ['trypsin', 'trypsin_p', 'lys_n', 'asp_n', 'arg_c', 'chymotrypsin', 'lys_c']
      if cy_enzyme.lower() not in allowed_enzyme:
        eprint("Error: Enzyme must be one of the following: " + ','.join(allowed_enzyme))
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)

    #Read input and obtain protid 
    if cy_debug:
      logging.debug("Step 1: Start processing the input protein list at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
    unique_each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict,initial_length, site_info_dict, ambigious_sites, unique_labels, dup_prot_ids, mult_mods_of_int = preprocessing(cy_in, cy_type_num, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out, database_dict, mods_list, cy_fasta_file, cy_enzyme, path_to_new_dir, logging_file, cy_fc_cutoff, cy_pval_cutoff, cy_ambi, cy_settings_file)
    
    # FC and Pval cutoff
    if (cy_type_num == "1" or cy_type_num == "2") and not (cy_fc_cutoff == 0.0 and cy_pval_cutoff == 1.0):
      unique_each_protein_list, prot_list, merged_out_dict = inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict)  
    
    if cy_type_num == "5" or cy_type_num == "6":
      for each_prot_id in merged_out_dict:
        comment_drop = False
        if not 'PickedPeptide' in merged_out_dict[each_prot_id]:
          comment_drop = True
        elif not merged_out_dict[each_prot_id]['PickedPeptide']:
          comment_drop = True
        if comment_drop:
          if 'CommentGene' in merged_out_dict[each_prot_id]:
            merged_out_dict[each_prot_id]['CommentGene'] += 'All peptides dropped;'
          else:
            merged_out_dict[each_prot_id].update({'CommentGene':'All peptides dropped;'})
    else:
      logging.debug(f"Remaining query: {len(unique_each_protein_list)} unique protein IDs")
            
    # Limit query inpt number = 1500
    if len(unique_each_protein_list) > 1500:
      eprint("Error: The query input is too big. Currently supporting upto 1500 query protein ids")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)
    
    if len(unique_each_protein_list) == 0:
      eprint("Error: No query protein ids found. Please check input, settings or filters")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1) 
    
    # open cytoscape
    subprocess.Popen([cy_exe, "-R", str(CYREST_PORT)])
    if gui_mode:
      print(f"COMMAND CYREST-PORT {str(CYREST_PORT)}")
    try:
      request_retry(f"{CYREST_URL}/v1/version", "GET")
    except CytoscapeError:
      eprint("Error: Cytoscape not responding. Please start the run again")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)

    # Check Cytoscape version
    request = request_retry(f'{CYREST_URL}/v1/version', 'GET')
    cy_version = request.json()
    check_cy_ver = re.match('^([0-9]{1,}\.[0-9]{1,})', cy_version['cytoscapeVersion'])
    if float(check_cy_ver.group(1)) < 3.7:
      eprint("Error: Cytoscape version must be 3.7.0 and above")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      sys.exit(1)

    # Start a new session 
    request_retry(f'{CYREST_URL}/v1/commands/session/new', 'GET')
    
    # Apps installed
    request = request_retry(f'{CYREST_URL}/v1/commands/apps/list installed', 'POST')
    apps_installed = request.json()
    app_reactome = False
    app_genemania = False
    app_cluego = False
    ver_reactome = ""
    ver_genemania = ""
    ver_cluego = ""
    
    for each_app in apps_installed['data']:
      app_name = each_app['appName']
      app_version = each_app['version']
      if app_name == "GeneMANIA":
        app_genemania = True
        ver_genemania = app_version
      elif app_name == "ClueGO":
        app_cluego = True
        ver_cluego = app_version
      elif app_name == "ReactomeFIPlugIn":
        app_reactome = True
        ver_reactome = app_version
    
    if cy_run.lower() == "genemania" or cy_run.lower() == "both":
      check_genemania_ver = re.match('^([0-9]{1,}\.[0-9]{1,})', ver_genemania)
      if float(check_genemania_ver.group(1)) < 3.5:
        eprint("Error: Cytoscape app GeneMANIA v3.5.0 or above not installed or not responding properly")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
        
    if not cy_cluego_inp_file:  
      check_cluego_ver = re.match('^([0-9]{1,}\.[0-9]{1,})', ver_cluego)
      if float(check_cluego_ver.group(1)) < 2.5:
        eprint("Error: Cytoscape app ClueGO v2.5.0 or above not installed or not responding properly")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)

      if ver_cluego not in cy_map:
        match = re.search(r"v(\d{1}\.\d{1}\.\d{1})", cy_map)
        eprint("Error: ClueGO version installed in Cytoscape is " + ver_cluego + " but ClueGO version selected in Setup page is " + match.group(1) + ". Please select ClueGO version " + ver_cluego + " in Setup page or reinstall ClueGO version " + match.group(1) + " in Cytoscape.")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      
      if cluego_reference_file and (ver_cluego == "2.5.0" or ver_cluego == "2.5.1" or ver_cluego == "2.5.2" or ver_cluego == "2.5.3" or ver_cluego == "2.5.4"):
        eprint("Error: Using ClueGO custom reference file needs version 2.5.5 or above of ClueGO")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
    
    if cy_run.lower() == "genemania" or cy_run.lower() == "both":
      body = dict(offline=True)
      response = request_retry(f"{CYREST_URL}/v1/commands/genemania/organisms", 'POST', json=body)
      genemania_bool = False
      for each in response.json()['data']['organisms']:
        if tax_id == str(each['taxonomyId']):
          genemania_bool = True
      if not genemania_bool:
        eprint("Error: Please install " + cy_species + " dataset in Genemania")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      
    #Uniprot API call to get primary gene, synonym
    if cy_debug:
      logging.debug("\nStep 2: Start the uniprot api call at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
      logging.debug("Uniprot query: " + str(len(unique_each_protein_list)))
        
    uniprot_query,each_primgene_list,merged_out_dict,ambigious_genes = uniprot_api_call(unique_each_protein_list, prot_list, cy_type_num, cy_debug, logging, merged_out_dict, organism_name, cy_session, cy_out, cy_cluego_out, cy_cluego_inp_file, path_to_new_dir, logging_file, site_info_dict, cy_ambi, cy_settings_file)
    
    all_prot_site_snps = {}
         
    drop_dupeprimgene_prot = {}
    unique_each_primgene_list = each_primgene_list
    genes_before_initial_drop = unique_each_primgene_list
    if len(each_primgene_list) != len(set(each_primgene_list)):
      if cy_ambi:
        warning = []
        count_dup_drops = 0
        dupe_gene_list = [item for item, count in collections.Counter(each_primgene_list).items() if count > 1]
        unique_each_primgene_list = [x for x in each_primgene_list if x.lower() not in [name.lower() for name in dupe_gene_list]]
        for each_dup_gene in dupe_gene_list:
          drop_dupeprimgene_prot_dict = get_query_from_list(uniprot_query, [each_dup_gene])
          drop_dupeprimgene_prot = list(drop_dupeprimgene_prot_dict.keys())
          count_dup_drops += len(drop_dupeprimgene_prot)
          for each_dropped_prot in drop_dupeprimgene_prot:
            if 'CommentGene' in merged_out_dict:
              merged_out_dict[each_dropped_prot]['CommentGene'] += "Duplicate primary gene;"
            else:
              merged_out_dict[each_dropped_prot].update({'CommentGene':'Duplicate primary gene;'})
          warning.append(each_dup_gene + "(" + ",".join(drop_dupeprimgene_prot) + ")")
        if cy_debug:
          logging.debug("AMBIGUITY WARNING - Uniprot duplicate primary gene mapping: " + str(count_dup_drops))
          logging.warning("Dropping queries: " + ','.join(warning))
      else:
        unique_each_primgene_list = list(set(each_primgene_list))
        warning = []
        count_dup_drops = 0
        all_dup_genes = [item for item, count in collections.Counter(each_primgene_list).items() if count > 1]
        for each_dup_gene in all_dup_genes:
          ambigious_genes.append(each_dup_gene)
          each_dup_prot_dict = get_query_from_list(uniprot_query, [each_dup_gene])
          each_dup_prot = list(each_dup_prot_dict.keys())
          retained_prot = ""
          all_dropped_prot = []
          bool1 = 0
          for get_each_dup_prot in each_dup_prot:
            if uniprot_query[get_each_dup_prot]['Reviewed'].lower() == "reviewed" and bool1 == 0:
              retained_prot = get_each_dup_prot
              bool1 = 1
            else:
             all_dropped_prot.append(get_each_dup_prot)
          if not retained_prot:          
            all_dropped_prot = each_dup_prot[1:]
            retained_prot = each_dup_prot[0]
            
          count_dup_drops += len(all_dropped_prot)
          merged_out_dict[retained_prot]['Primary'] = merged_out_dict[retained_prot]['Primary'] + "**"
          for each_dropped_prot in all_dropped_prot:
            del uniprot_query[each_dropped_prot]
            if 'CommentGene' in merged_out_dict:
              merged_out_dict[each_dropped_prot]['CommentGene'] += "Duplicate primary gene;"
            else:
              merged_out_dict[each_dropped_prot].update({'CommentGene':'Duplicate primary gene;'})
          warning.append(each_dup_gene + "(" + ",".join([retained_prot]+all_dropped_prot) + ")")
        if cy_debug:
          logging.debug("AMBIGUITY WARNING - Uniprot duplicate primary gene mapping: " + str(count_dup_drops))
          logging.warning("Dropping all but first query: " + ','.join(warning))            
        
    if cy_cluego_inp_file:
      leading_term_cluster, unique_each_primgene_list = cluego_input_file(cy_cluego_inp_file, cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      cy_lim = 0
      if cy_debug:
        logging.debug("Limiting query to re-analyze terms: " + str(len(unique_each_primgene_list)) )
  
    if not unique_each_primgene_list:
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      eprint("Error: No query for String and Genemania")
      sys.exit(1)
    
    string_filtered_dict = {}
    genemania_filtered_dict = {}
    
    if cy_run.lower() == "string" or cy_run.lower() == "both":
      if cy_debug:
        logging.debug("\nStep 3: String mapping started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
        logging.debug("String query: " + str(len(unique_each_primgene_list)))
 
      #Send gene list to string, get mapping & interactions
      string_interaction, string_unique_vertex, string_mapping, merged_out_dict = create_string_cytoscape(uniprot_query,unique_each_primgene_list, tax_id, cy_lim, cy_score, cy_debug, logging, merged_out_dict, genes_before_initial_drop, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
      #Categorize interaction nodes as primary gene, secondary gene, external synonym and external interactor
      string_category = categorize_gene(string_unique_vertex, string_mapping, uniprot_query)
      #Get a filtered dictionary of interactions(primary-primary; primary-external interactor; external interactor-external interactor)
      string_filtered_dict, merged_out_dict = get_search_dicts(string_interaction, string_category, logging, cy_debug, uniprot_query, string_mapping, merged_out_dict, "String")
      # Update merged_out_dict
      merged_out_dict = get_interactions_dict(string_filtered_dict, 'String', merged_out_dict)

    if cy_run.lower() == "genemania" or cy_run.lower() == "both":
      if cy_debug:
        logging.debug("\nStep 4: Genemania mapping started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
        logging.debug("Genemania query: " + str(len(unique_each_primgene_list)))
    
      #Send gene list to string, get mapping & interactions
      genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict = create_genemania_interactions(uniprot_query,unique_each_primgene_list,tax_id,cy_lim,"10",cy_debug,logging, merged_out_dict, cy_session, cy_out, cy_cluego_out, genes_before_initial_drop, path_to_new_dir, logging_file, cy_settings_file)
      #Categorize interaction nodes as primary gene, secondary gene, external synonym and external interactor
      genemania_category = categorize_gene(genemania_unique_vertex, genemania_mapping, uniprot_query)
      #Get a filtered dictionary of interactions(primary-primary; primary-external interactor; external interactor-external interactor)
      genemania_filtered_dict, merged_out_dict = get_search_dicts(genemania_interaction, genemania_category, logging, cy_debug, uniprot_query, genemania_mapping, merged_out_dict, "Genemania")
      # Update merged_out_dict
      merged_out_dict = get_interactions_dict(genemania_filtered_dict, 'Genemania', merged_out_dict)
    
    #Merge String and Genemania interactions
    if cy_debug:
      logging.debug("\nStep 5: Merge interactions started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
    
    interaction_skip = False
    if not string_filtered_dict and not genemania_filtered_dict:
      if not cy_cluego_inp_file: 
        eprint("Error: No interactions found in String and Genemania")
        remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
        sys.exit(1)
      else:
        interaction_skip = True
   
    unique_merged_interactions = []
    unique_nodes = []
    #Get a list of merged interactions (ex: [node1 node2, node3 node4] etc) and a list of unique nodes (ex: [node1, node2, node3, node4])
    if (cy_run.lower() == "string" or cy_run.lower() == "both") and not interaction_skip:
      unique_merged_interactions,unique_nodes = get_merged_interactions(string_filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, cy_type_num)
    if (cy_run.lower() == "genemania" or cy_run.lower() == "both") and not interaction_skip:
      unique_merged_interactions,unique_nodes = get_merged_interactions(genemania_filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, cy_type_num)
    
    for each_merged_dict_node in merged_out_dict:
      if 'Primary' in merged_out_dict[each_merged_dict_node] and merged_out_dict[each_merged_dict_node]['Primary']:
        if ((merged_out_dict[each_merged_dict_node]['Primary'].lower()).replace("**","")) not in unique_nodes:
          if 'CommentGene' in merged_out_dict[each_merged_dict_node]:
            if not merged_out_dict[each_merged_dict_node]['CommentGene']:
              merged_out_dict[each_merged_dict_node]['CommentGene'] += 'No interactions found;'
          else:
            merged_out_dict[each_merged_dict_node].update({'CommentGene':'No interactions found;'})
          
    if cy_debug:
      logging.debug("Total merged query nodes: " + str(len([i for i in unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])))
      logging.debug("Total merged interactions: " + str(len(unique_merged_interactions)))
    
    #Get uniprot query, primary gene and FC values together
    uniprot_list = {}

    if not cy_cluego_inp_file:
      for each_node in unique_nodes:
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num, site_info_dict, ambigious_sites, ambigious_genes)    
    else:
      lower_unique_each_primgene_list = [x.lower() for x in unique_each_primgene_list]
      for each_node in lower_unique_each_primgene_list:
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num, site_info_dict, ambigious_sites, ambigious_genes)
    
    #Interactors styling   
    if not interaction_skip: 
      if not (cy_type_num == "5" or cy_type_num == "6"):         
        cy_interactors_style(unique_nodes, unique_merged_interactions, uniprot_list, max_FC_len, each_category, cy_pval, unique_labels)
        
      else:    
        cy_sites_interactors_style(unique_nodes, unique_merged_interactions, uniprot_list, max_FC_len, each_category, cy_pval, cy_type_num, all_prot_site_snps, uniprot_query, unique_labels,mult_mods_of_int)
    
    #Category styling   
    if cy_type_num == "4":    
      cy_category_style(unique_nodes, unique_merged_interactions, uniprot_list, each_category)
    
    coverage = 0.0
    
    if not cy_cluego_inp_file:
      if cy_debug:
        logging.debug("\nStep 6: ClueGO started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
  
      #Start ClueGO
      response = request_retry(f'{CYREST_URL}/v1/apps/cluego/start-up-cluego', 'POST')
      
      if cy_debug:
        # Number of ClueGO query + EI = x + y
        logging.debug("ClueGO query + External Interactor: " + str(len([i for i in unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])) + " + " + str(len([i for i in unique_nodes if i.lower() not in [y.lower() for y in unique_each_primgene_list] ])))  

      filtered_unique_nodes, merged_out_dict = cluego_filtering(unique_nodes, cy_map, uniprot_query, cy_debug, logging, merged_out_dict, unique_each_primgene_list, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)

      if cy_debug:
        #Number of ClueGO query + EI = x + y
        logging.debug("Total ClueGO query + External Interactor: " + str(len([i for i in filtered_unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])) + " + " + str(len([i for i in filtered_unique_nodes if i.lower() not in [y.lower() for y in unique_each_primgene_list] ])))      
    
      final_length = len([i for i in filtered_unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])
      coverage = final_length/initial_length *100
      cluego_run(organism_name,cy_cluego_out,filtered_unique_nodes,cy_cluego_grouping,select_terms, leading_term_selection,cluego_reference_file,cluego_pval, cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file, cy_cluego_log_out, cy_type_num, uniprot_list, max_FC_len, unique_labels)
    
    if leading_term_cluster:
      cy_pathways_style(leading_term_cluster, each_category, max_FC_len, cy_pval, uniprot_list, cy_type_num, all_prot_site_snps, uniprot_query, unique_labels,mult_mods_of_int)
    
    if cy_debug:
      if not cy_cluego_inp_file:
        logging.debug("\nQuery coverage = " + str(round(coverage, 2)) + "%")
      logging.debug("\nRun completed successfully")
     
    #Write into outfile
    write_into_out(merged_out_dict, cy_out, dup_prot_ids)
    request_retry(f"{CYREST_URL}/v1/session?file=" + urllib.parse.quote_plus(cy_session), 'POST')

  except CytoscapeError as e:
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    eprint("Error: Cytoscape not responding. Please start the run again")
    sys.exit(1)

  except Exception as e:
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out, path_to_new_dir, logging_file, cy_settings_file)
    cytoscape_not_open_msg = "No connection could be made because the target machine actively refused it"
    cytoscape_not_open_msg2 = "Remote end closed connection without response"
    cytoscape_not_responding_msg = "Expecting value: line 1 column 1 (char 0)"
    if cytoscape_not_open_msg in str(e) or cytoscape_not_open_msg2 in str(e):
      eprint("Error: Cytoscape not responding. Please start the run again")
      sys.exit(1)
    elif cytoscape_not_responding_msg in str(e):
      eprint("Error: Cytoscape not responding. Please start the run again")
      sys.exit(1)
    else:
      traceback.print_exc()
      eprint("Fatal error")
      sys.exit(1)
      
if __name__ == "__main__":
  main(sys.argv[1:])
