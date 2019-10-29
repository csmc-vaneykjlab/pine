#usage: python get_interactions_table_4.py -i Mouse_rev_unrev_input.csv -t nofc -s rat -m "C:\Users\SundararamN\ClueGOConfiguration\v2.5.5\ClueGOSourceFiles\Organism_Rattus norvegicus\Rattus norvegicus.gene2accession_2019.05.20.txt.gz" -d
#Goal: Given list of protein IDs (with/without their corresponding FC and pval), construct 1) an interaction network among all proteins in list 2) construct a pathway network of proteins in the list 
try:
  from py2cytoscape import cyrest
except ImportError:
  eprint("Error: Please install module py2cytoscape. [Installation: pip install py2cytoscape]")
from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
try:
  import requests
except ImportError:
  eprint("ImportError: Please make sure you are using Python3")
from urllib.parse import quote
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
import sys
import time
import re
import logging
import gzip
import math
import igraph
import traceback
import collections
import subprocess

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

def setup_logger(name, log_file, level=logging.DEBUG, with_stdout=False):
  handler = logging.FileHandler(log_file,mode='w')
  logger = logging.getLogger(name)
  logger.setLevel(level)
  logger.addHandler(handler)
  if with_stdout:
    logger.addHandler(logging.StreamHandler(sys.stdout))
  return logger

def db_handling(db_file):
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
          seq += line.strip("\n\r")
            
      if proteinID:
        if '|' in proteinID:
          proteinID = proteinID.split('|')[1]
      
        if proteinID not in database_dict:
          database_dict[proteinID] = [seq]
        else:
          database_dict[proteinID].append(seq)
     
        proteinID = line.strip(">\n").split(" ")[0]
  return(database_dict)
          
def find_mod(seq): #this is for find all the modification position in the peptide
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

def find(seq,database):#put \ before [ and ] this is for find the peptide (which remove the modification) in the database
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
  
def preprocessing(inp, type, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out, database_dict, include_list, db_file, enzyme):
  is_prot_col = False
  is_FC = False
  is_pval = False
  is_cat = False
  is_multFC = False
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
  to_return_unique_protids_length = 0
  ambigious_sites = {}
  pep_to_mod = {}
  with open(inp,'r') as csv_file:
    '''
    Read input and collect columns based on type of analysis chosen
    Types:1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = category; 5=singlefc-ptm; 6=multifc-ptm
    '''
    input_file = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in input_file:            
      if line_count == 0:
        header = row
        # Get header
        for i in range(len(row)):                                                                               
          if "proteinid" in row[i].lower():
            protein = i
            is_prot_col = True
          if type == "2" or type == "1":
            if "fc" in row[i].lower():
              FC=i
              is_FC = True  
            elif "adj.pval" in row[i].lower():
              pval = i
              is_pval = True
            elif "pval" in row[i].lower() and not is_pval:
              pval=i
              is_pval = True  
          if type == "4":
            if "category" in row[i].lower():
              cat = i
              is_cat = True
          if type == "2" or type == "6":
            is_multFC = True
            if "label" in row[i].lower():
              label = i  
              is_label_col = True
          if type == "5" or type == "6":
            if "peptide" in row[i].lower():
              peptide_col = i
              is_pep_col = True
            elif "fc" in row[i].lower():
              FC = i
              is_FC = True          
        
        if not is_prot_col:
          eprint("Error: Column 'ProteinID' not found")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit(1)
      else: 
        if '/' in row[protein]:
          split_val = row[protein].split('/')
          if len(split_val) > 2:
            eprint("Error: Multiple protein IDs not considered")
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
            sys.exit(1)
          elif len(split_val) == 2:
            if not bool(re.match('^[0-9]$', split_val[0])):
              eprint("Error: Multiple protein IDs not considered")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
              sys.exit(1)
        
        if '|' in row[protein]:
          row[protein] = row[protein].split('|')[1]
          
        if not bool(re.match('^[A-Za-z0-9\-]+$', row[protein])):
          err_msg = "Error: Invalid proteinID: " + row[protein]
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit(1)
        
        if type == "1" or type == "2" or type == "5" or type == "6":
          skip = False
          if not is_pval:
            get_pval = 1.0
          else:
            try:
              float(row[pval])
              get_pval = row[pval]
            except ValueError:
              get_pval = 1.0
              
          if is_FC:
            try:
              float(row[FC])
              get_fc_val = row[FC]
            except ValueError:
              get_fc_val = 0.0
        
        if type == "5" or type == "6":
          all_mods_for_prot = []
          if include_list:
            modInSeq_dict = {}
            if is_pep_col:
              peptide = row[peptide_col]
            else:
              print("Error: Column 'Peptide' not found")
              sys.exit(1)
            protein_list_id = row[protein]
            modInSeq_all_dict = find_mod(peptide)
            combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
            peptide_sub = re.sub(combined_pat, '', peptide)  #remove modification from peptide sequence
            peptide_sub = peptide_sub.strip('"')		           
            if protein_list_id in database_dict:
              for each_seq_in_db_dict in database_dict[protein_list_id]:
                seqInDatabase_list = find(peptide_sub,each_seq_in_db_dict)
                seqInDatabase = seqInDatabase_list[0]
                break
            else:
              eprint("Error: ProteinID '" + str(protein_list_id) + "' not found in " + str(db_file) + ". Check fasta database file provided.")  
              sys.exit(1)            
            for k in include_list:
              for key,value in modInSeq_all_dict.items():                   
                combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
                key = re.sub('\+','',key)                               
                k = re.sub('\+','',k)
                key_match = re.search(combined_pat,key)
                k_match = re.search(combined_pat,k)
                if not (key_match and k_match):
                  key1 = re.sub(combined_pat, '', key)
                  k1 = re.sub(combined_pat,'',k)               
                if k1.lower() in key1.lower(): 
                  for each_val in value:
                    val = int(each_val)+int(seqInDatabase)+1
                    if k1 in modInSeq_dict:            
                      modInSeq_dict[k1].append(val)                          
                    else:
                      modInSeq_dict[k1] = [val]
                    all_mods_for_prot.append(key)
            sites = ""
            sites_list = []
            for key,value in modInSeq_dict.items():
              for each_val in value:
                sites += key + str(each_val) + "/"
            sites = sites.strip("/")
            if sites: 
              if type == "5":
                max_FC_len = 1
                if is_prot_col and is_FC:
                  if protein_list_id not in site_info_dict:
                    site_info_dict[protein_list_id] = {}
                    site_info_dict[protein_list_id][sites] = {}
                    site_info_dict[protein_list_id][sites].update({peptide:[[float(get_fc_val)],[get_pval],all_mods_for_prot]})
                  else:
                    if sites not in site_info_dict[protein_list_id]:
                      site_info_dict[protein_list_id][sites] = {}
                      site_info_dict[protein_list_id][sites].update({peptide:[[float(get_fc_val)],[get_pval],all_mods_for_prot]})
                    else:
                      if not peptide in site_info_dict[protein_list_id][sites]:
                        site_info_dict[protein_list_id][sites].update({peptide:[[float(get_fc_val)],[get_pval],all_mods_for_prot]})
                      
                else:
                  eprint("Error: Required columns- ProteinID, Peptide, FC")
                  remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                  sys.exit(1)
              
              if type == "6":              
                if is_prot_col and is_FC and is_label_col:
                  if protein_list_id not in site_info_dict:
                    site_info_dict[protein_list_id] = {}
                    site_info_dict[protein_list_id][sites] = {}
                    site_info_dict[protein_list_id][sites].update({peptide: [[float(get_fc_val)],[get_pval],all_mods_for_prot,[row[label]]] })
                  else:
                    if sites not in site_info_dict[protein_list_id]:
                       site_info_dict[protein_list_id][sites] = {}
                       site_info_dict[protein_list_id][sites].update({peptide: [[float(get_fc_val)],[get_pval],all_mods_for_prot,[row[label]]] })
                    else:
                      if peptide not in site_info_dict[protein_list_id][sites]:
                        site_info_dict[protein_list_id][sites].update({peptide: [[float(get_fc_val)],[get_pval],all_mods_for_prot,[row[label]]] })
                      else:
                        if row[label] not in site_info_dict[protein_list_id][sites][peptide]:
                          site_info_dict[protein_list_id][sites][peptide][0].append(float(get_fc_val))
                          site_info_dict[protein_list_id][sites][peptide][1].append(float(get_pval))
                          site_info_dict[protein_list_id][sites][peptide][3].append(row[label])
                        #else dup protid, sites, peptide, label
                        
                  if row[label] not in unique_labels:
                    unique_labels.append(row[label])                
                else:
                  eprint("Error: Required columns- ProteinID, Peptide, FC, Label")
                  remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                  sys.exit(1)
                  
        if row[protein] not in each_protein_list:
          each_protein_list.append(row[protein])
          if type == "1" or type == "2" or type == "3":          
            if type == "1":
              if is_prot_col and is_FC:           
                prot_list.update({row[protein]:[[float(get_fc_val)],[get_pval]]})
                max_FC_len = 1
              else:
                eprint("Error: Required columns- ProteinID, FC")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit(1)
                
            elif type == "2":
              if is_prot_col and is_FC and is_label_col:
                prot_list[row[protein]] = {}
                prot_list[row[protein]].update({row[label]:[float(get_fc_val),get_pval]})
                if row[label] not in unique_labels:
                  unique_labels.append(row[label])                
              else:
                eprint("Error: Required columns- ProteinID, FC, Label")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit(1)
                
            elif type == "3":
              if is_prot_col:
                continue
              else:
                eprint("Error: Required columns- ProteinID")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit(1)
            
          elif type == "4":
            if (not is_prot_col) or (not is_cat):
              eprint("Error: Required columns- ProteinID, Category")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
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

  if len(unique_labels) > 10:
    eprint("Error: Number of unique labels should not exceed 10")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit(1)
    
  if len(each_category) > 10:
    eprint("Error: Number of categories should not exceed 10")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit(1)
  
  if cy_out:
    csv_file = open(cy_out,'w')
    csv_file.write("ProteinID,Primary Gene,String,Genemania,Comment,\n")

  if cy_debug:
    logging.debug("Initial query: " + str(len(each_protein_list)))
  initial_length = len(each_protein_list)
  to_return_unique_protids_length = len(set(each_protein_list))

  if type == "1": 
    dropping_repeats = []    
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
      if all_dropped:       
        logging.debug("Duplicate query: " + str((initial_length)-len(each_protein_list)))
        logging.warning("WARNING - Dropping queries: " + ','.join(all_dropped))
    
    for each_dupe_query in all_dropped:
      line = each_dupe_query + ",,,,," + "Duplicate query;\n"
      if cy_out:
        csv_file.write(line)  
      
  elif type == "3":
    if repeat_prot_ids:
      unique_each_protein_list = list(set(each_protein_list))    
      if cy_debug:
        logging.debug("Duplicate query: " + str(len(each_protein_list)-len(unique_each_protein_list)))
        logging.warning("WARNING - Dropping queries: " + ','.join(repeat_prot_ids)) 
      for each_dupe_query in repeat_prot_ids:
        line = each_dupe_query + ",,,," + "Duplicate query;\n"
        if cy_out:
          csv_file.write(line) 
      each_protein_list = unique_each_protein_list      

  elif type == "4":
    unique_each_protein_list = list(set(each_protein_list))
    if retain_prot_ids:
      if cy_debug:
        logging.debug("Duplicate query: " + str(len(retain_prot_ids)))
        logging.warning("WARNING - Dropping queries: " + ','.join(retain_prot_ids)) 
      for each_dupe_query in retain_prot_ids:
        line = each_dupe_query + ",,,," + "Duplicate query;\n"
        if cy_out:
          csv_file.write(line) 
    each_protein_list = unique_each_protein_list

  elif type == "2":
      unique_each_protein_list = list(set(each_protein_list))
      count_dropped = 0
      additional_dropped = []
      if retain_prot_ids:
        for each_dupe_query in retain_prot_ids:
          line = each_dupe_query + ",,,," + "Duplicate query;\n"

          if cy_out:
            csv_file.write(line)
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
        list_of_duplicates = additional_dropped + list(repeat_prot_ids_2.keys()) + retain_prot_ids
        list_of_duplicates = sorted(list_of_duplicates)
        if list_of_duplicates:
          logging.debug("Duplicate query: " + str(len(list_of_duplicates)))
          logging.warning("WARNING - Dropping queries: " + ','.join(list_of_duplicates)) 
      
      for each_dupe_query in list_of_duplicates:
        line = each_dupe_query + ",,,," + "Duplicate query;\n"
        if cy_out:
          csv_file.write(line)
    
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
        keys = list(site_info_dict[each_protid][each_site].keys())
        if len(keys) > 1:
          all_mods = [v[2] for k,v in site_info_dict[each_protid][each_site].items()]
          unique_all_mods = [list(x) for x in set(tuple(x) for x in all_mods)]
          if len(unique_all_mods) > 1:
            non_unique = True
            for each_key_pep in site_info_dict[each_protid][each_site]:
              new_each_site = ""
              for each_modsite in site_info_dict[each_protid][each_site][each_key_pep][2]:
                new_each_site_num = re.sub("[^0-9]", "", each_modsite)
                match = re.match(r"([A-Za-z]+)([0-9]+)", each_site)
                if match:
                  if new_each_site:
                    new_each_site += "/"
                  new_each_site += match.group(1) + "{" + new_each_site_num + "}" + match.group(2)
              new_each_site_info = site_info_dict[each_protid][each_site][each_key_pep]
              
              if each_protid not in site_info_dict_rearrange:
                site_info_dict_rearrange[each_protid] = {}
                site_info_dict_rearrange[each_protid].update({new_each_site:new_each_site_info})
              else:
                site_info_dict_rearrange[each_protid].update({new_each_site:new_each_site_info})  
         
          else:              
            each_site_info, dropped_pep = ptm_scoring(site_info_dict[each_protid][each_site], enzyme, include_list)
            for each_dropped_pep in dropped_pep:
              count_dropped+=1
              if each_protid in all_dropped_pep:
                all_dropped_pep[each_protid].append(each_dropped_pep)
                all_dropped_warning += ", " + each_dropped_pep
                ambigious_sites[each_protid].append(each_site)                
              else:
                all_dropped_pep.update({each_protid:[each_dropped_pep]})
                all_dropped_warning += each_protid + "(" + each_dropped_pep
                ambigious_sites.update({each_protid:[each_site]})
            if all_dropped_warning:
              all_dropped_warning += ")"
            
        else:
          each_site_info = site_info_dict[each_protid][each_site][keys[0]]
        if not non_unique:
          if each_protid not in site_info_dict_rearrange:
            site_info_dict_rearrange[each_protid] = {}
            site_info_dict_rearrange[each_protid].update({each_site:each_site_info})
          else:
            site_info_dict_rearrange[each_protid].update({each_site:each_site_info})
    site_info_dict = site_info_dict_rearrange     
    if cy_debug and all_dropped_pep:
      logging.debug("Lower score query: " + str(count_dropped))
      logging.warning("WARNING - Dropping queries: " + all_dropped_warning)
           
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
    logging.debug("Unique remaining queries: " + str(len(each_protein_list)))
  if cy_out:    
    csv_file.close() 
  
  return(each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict, to_return_unique_protids_length, site_info_dict, ambigious_sites)

def ptm_scoring(site_dict, enzyme, include_list):
  enzyme_info = {'trypsin':{'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP', 'RP']}, 
                 'trypsin_p':{'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : []}, 
                 'lys_n':{'terminus' : 'N' , 'cleave' : ['K'], 'exceptions' : []}, 
                 'asp_n':{'terminus' : 'N' , 'cleave' : ['B','D'], 'exceptions' : []}, 
                 'arg_c':{'terminus' : 'C' , 'cleave' : ['R'], 'exceptions' : ['RP']}, 
                 'chymotrypsin':{'terminus' : 'C' , 'cleave' : ['F','Y','W','L'], 'exceptions' : ['FP','YP','WP','LP']}, 
                 'lys_c' : {'terminus' : 'C' , 'cleave' : ['K'], 'exceptions' : ['KP']}}
                 
  all_peptides = list(site_dict.keys())
  top_score = [0] * len(all_peptides)
 
  for each_peptide in all_peptides: 
    # Other Mods
    total_mods = 0
    om = 0
    all_mods_dict = find_mod(each_peptide)
    for k in include_list:
      for key,value in all_mods_dict.items():  
        combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))
        key = re.sub('\+','',key)                               
        k = re.sub('\+','',k)
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
    # count miscleavages
    for each_cleave in enzyme_info[enzyme]['cleave']:
      miscleave += each_peptide.count(each_cleave)
    # Exclude exceptions
    for each_except in enzyme_info[enzyme]['exceptions']:
      exception += each_peptide.count(each_except)
    index = all_peptides.index(each_peptide)
    
    #Score update
    top_score[index] += miscleave - exception + om 
  # Pick lowest score
  min_score = min(top_score)
  indices = [i for i, x in enumerate(top_score) if x == min_score]
  largest_avg_FC = 0
  index = 0
  for each_index in indices:
    each_pick_pep = all_peptides[each_index]
    n = 0
    avg_FC = 0     
    for each_site_pos in site_dict[each_pick_pep][0]:
      avg_FC += abs(each_site_pos)
      n+= 1
    avg_FC = avg_FC/n
    if abs(avg_FC) >= largest_avg_FC:     
      largest_avg_FC = abs(each_site_pos)
      index = each_index
  pick_pep = all_peptides[index]
  dropped_pep = [x for x in all_peptides if not x == pick_pep]
  
  return(site_dict[pick_pep],dropped_pep)  
    
def inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict):
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
      merged_out_dict[each_prot] = {}
      merged_out_dict[each_prot].update({'Primary':'', 'Comment':'FC/Pval cutoff not met' , 'String':'', 'Genemania':'', 'ClueGO':''})
      
  if cy_debug:
    logging.debug("FC and PVal cutoff not met: " + str(len(queries_dropped)))
    logging.warning("WARNING - Dropping queries: " + ','.join(queries_dropped))
    
  return(unique_each_protein_list, prot_list, merged_out_dict)
 
def uniprot_api_call(each_protein_list, prot_list, type, cy_debug, logging, merged_out_dict, species, cy_session, cy_out, cy_cluego_out):
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

  params = {
  'from':'ACC+ID',
  'to':'ACC',
  'format':'tab',
  'query':query_term,
  'columns': 'id,genes(PREFERRED),genes(ALTERNATIVE),organism,feature(NATURAL VARIANT),last-modified'
  }
  data = urllib.parse.urlencode(params).encode("utf-8")
  request = urllib2.Request(url, data)
  response = urllib2.urlopen(request)
  page = response.read(200000)
  decode = page.decode("utf-8")
  list1=decode.split('\n')
  list1 = list1[1:]
  for each_list1 in list1:
    remaining_isoforms = []
    if each_list1 != '':
      uniprot_list = each_list1.split('\t')
      uniprot_protid = uniprot_list[0]
      prot = uniprot_list[6]
      split_prot_list = prot.split(',')
      if len(split_prot_list) > 0:
        remaining_isoforms = split_prot_list[1:]
        all_isoforms.extend(remaining_isoforms)
      each_prot = split_prot_list[0]
      primary_gene = uniprot_list[1]
      merged_out_dict[each_prot] = {}
      comment_merged = ""
      # Pick first gene in case of multiple primary genes for single uniprot ids 
      if ";" in primary_gene:
        prot_with_mult_primgene.append(each_prot + "(" + primary_gene + ") ")
        #comment_merged = "Multiple primary genes;"
        primary_gene = primary_gene.split(";")[0]
        if primary_gene not in ambigious_gene:
          ambigious_gene.append(primary_gene)
      if remaining_isoforms:
        if primary_gene not in ambigious_gene:
          ambigious_gene.append(primary_gene) 
      synonym_gene = uniprot_list[2]
      synonym_gene = synonym_gene.split(" ")
      organism_name = uniprot_list[3]
      uniprot_query[each_prot] = {}
      uniprot_query[each_prot].update({'Uniprot':uniprot_protid,'Primary':primary_gene,'Synonym':synonym_gene,'Organism':organism_name,'Natural_variant':uniprot_list[4],'Date_modified':uniprot_list[5]})
      # Do not add empty primary gene to list
      if primary_gene: # Duplicates in primary gene  
        each_primgene_list.append(primary_gene)
      else:
        no_primgene_val.append(each_prot)
        comment_merged = "Primary gene unavailable;"
          
      if uniprot_protid:
        each_protid_list.append(uniprot_protid)
        
      merged_out_dict[each_prot].update({'Primary':primary_gene, 'Comment':comment_merged, 'String':'', 'Genemania':'', 'ClueGO':''})

  for each_prot_in_input in each_protein_list:
    if each_prot_in_input not in uniprot_query:
      no_uniprot_val.append(each_prot_in_input)
      uniprot_query[each_prot_in_input] = {}
      uniprot_query[each_prot_in_input].update({'Uniprot':"NA",'Primary':"NA",'Synonym':"NA",'Organism':"NA"})
      merged_out_dict[each_prot_in_input] = {}
      comment_merged = "Uniprot query not mapped;" 
      merged_out_dict[each_prot_in_input].update({'Primary':'', 'Comment':comment_merged, 'String':'', 'Genemania':'', 'ClueGO':''})
  
  organisms = [dict['Organism'] for dict in uniprot_query.values() if dict['Organism'] != "NA"]
  unique_organisms = list(set(organisms))

  if len(unique_organisms) > 1:
    eprint("Error: Protein list is of more than 1 organism: " + ','.join(unique_organisms))
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit(1)
  
  if not species.lower() in unique_organisms[0].lower():
    eprint("Error: Species mismatch. Species parameter provided is " + species + " and species of protein list is " + unique_organisms[0])
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit(1)
   
  if cy_debug:
    if prot_with_mult_primgene:
      logging.warning("WARNING - Uniprot Multiple primary genes: " + ','.join(prot_with_mult_primgene))
  
    if no_uniprot_val:
      logging.debug("Uniprot query not mapped: " + str(len(no_uniprot_val)))
      logging.warning("WARNING - Dropping queries: " + ','.join(no_uniprot_val))

    if no_primgene_val:
      logging.debug("Uniprot Primary gene unavailable: " + str(len(no_primgene_val)))
      logging.warning("WARNING - Dropping queries: " + ','.join(no_primgene_val))
  
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

def get_query_from_list(uniprot_query, list):
  dropped_list = {}
  for key in uniprot_query:
    for each_dupe_gene in list:
      if uniprot_query[key]['Primary'].lower() == each_dupe_gene.lower():
        dropped_list.update({key:each_dupe_gene})
        break
  return(dropped_list)

def check_dup_preferred_gene(mapping, interaction, unique_vertex, search, uniprot_query, cy_debug, logging, merged_out_dict):
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
        comment_merged = str(search) + "- not mapped;"
        merged_out_dict[each_get_dupe_ids]['Comment'] += comment_merged
        warning.append(each_get_dupe_ids + "(" + get_dupe_query_proids[each_get_dupe_ids] + ") ") 
    
    if len(get_dupe_query_proids_2) != 0:    
      for each_get_dupe_ids in get_dupe_query_proids_2:
        comment_merged = str(search) + "- not mapped;"
        merged_out_dict[each_get_dupe_ids]['Comment'] += comment_merged
        warning.append(each_get_dupe_ids + "(" + get_dupe_query_proids_2[each_get_dupe_ids] + ") ")  
      
  else:
    return(interaction, unique_vertex, dupe_preferred_list, merged_out_dict, warning)
    
  return(dropped_interaction_list, unique_each_preferred_list, dupe_preferred_list, merged_out_dict, warning)
  
def overwrite_values(overwrite_preferred, interaction, unique_vertex, mapping, merged_out_dict):
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
  
def create_string_cytoscape(uniprot_query,each_inp_list, species, limit, score, cy_debug, logging, merged_out_dict, genes_before_initial_drop):
  string_db_out = {}
  string_mapping = {}
  string_list_input = "\n".join(each_inp_list)
  
  data = {
  'identifiers': string_list_input,
  'species': species,
  'echo_query': 1,
  "limit" : 1
  }
  
  response = requests.post('https://string-db.org/api/tsv-no-header/get_string_ids', data=data)
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
  string_db_out.update(pd.read_table(StringIO(response.text), header=None))
  string1 = list(string_db_out[2].values.flatten())
  string2 = list(string_db_out[3].values.flatten())
  
  string_interaction,string_unique_vertex = remove_duplicates_within(string1,string2)
 
  no_mapping = [i for i in each_inp_list if i.lower() not in [x.lower() for x in string_mapping]]
  no_mapping_prots =  get_query_from_list(uniprot_query, no_mapping)

  no_interactions = [i for i in list(string_mapping.keys()) if string_mapping[i].lower() not in [x.lower() for x in string_unique_vertex] ]
  no_interactions_prots = get_query_from_list(uniprot_query, no_interactions)
  # Drop duplicate Preferred genes
  string_interaction,string_unique_vertex,string_dupe_preferred,merged_out_dict,dup_pref_warning = check_dup_preferred_gene(string_mapping,string_interaction,string_unique_vertex,"String", uniprot_query, cy_debug, logging, merged_out_dict)
  
  if overwrite_preferred:
    string_interaction,string_unique_vertex,string_mapping,merged_out_dict = overwrite_values(overwrite_preferred, string_interaction, string_unique_vertex, string_mapping, merged_out_dict)
  
  if len(genes_before_initial_drop) != len(each_inp_list):
    string_interaction,string_unique_vertex = drop_ei_if_query(string_interaction,string_unique_vertex,genes_before_initial_drop,each_inp_list)
    
  if cy_debug: 
    if len(no_mapping) != 0:
      logging.debug("String queries not mapped: " + str(len(no_mapping) + len(dup_pref_warning)))
      warning = []
      for each in no_mapping_prots:
        merged_out_dict[each]['Comment'] += "String- not mapped;"
        warning.append(each + "(" + no_mapping_prots[each] + ") ")
      if dup_pref_warning:
        logging.debug("WARNING - Dropping queries: " + ','.join(warning) + "," + ','.join(dup_pref_warning))
      else:
        logging.debug("WARNING - Dropping queries: " + ','.join(warning))
    
    if len(no_interactions) != 0:
      logging.debug("String queries with no interactions: " + str(len(no_interactions)))
      warning = []
      for each in no_interactions_prots:
        merged_out_dict[each]['Comment'] += "String- no interaction;"
        warning.append(each + "(" + no_interactions_prots[each] + ") ")
      logging.debug("WARNING - Dropping queries: " + ','.join(warning))

  return(string_interaction, string_unique_vertex, string_mapping, merged_out_dict)
  
def remove_duplicates_within(node1,node2):
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

def create_genemania_interactions(uniprot_query,each_inp_list,species,limit,att_limit,cy_debug,logging,merged_out_dict,cy_session, cy_out, cy_cluego_out, genes_before_initial_drop):
  '''
  Genemania interactions for primary genes
  Remove duplicates
  '''
  genemania_interaction = []
  genemania_unique_vertex = []
  genemania_mapping = {}
  genemania_node = {}
  overwrite_preferred = {}
  
  try:
    join_genes = '|'.join(each_inp_list)
    #body = dict(attrLimit=str(att_limit), geneLimit=str(limit), genes=join_genes, organism=species)
    body = dict(attrLimit=str(att_limit), geneLimit=str(limit), genes=join_genes, organism=species, offline=True)
    get_genemania = requests.post('http://localhost:1234/v1/commands/genemania/search', json=body)
    uploaded_list = get_genemania.json()
    current_network_suid = str(uploaded_list['data']['network'])
    request = 'http://localhost:1234/v1/networks/' + str(uploaded_list['data']['network']) + '/tables/defaultedge'
    resp = requests.get(request, json=body)
    edge_info =resp.json()
    request = 'http://localhost:1234/v1/networks/' + str(uploaded_list['data']['network']) +'/tables/defaultnode'
    resp = requests.get(request, json=body)
    node_info = resp.json()
    requests.delete("http://localhost:1234/v1/networks/" + str(current_network_suid))
    
  except Exception as e:
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)  
    try:
      requests.get("http://localhost:1234/v1/commands/command/quit")
      eprint("Error: Genemania timed out- please try again later")
      sys.exit(1)
    except:
      eprint("Error: Cytoscape must be open")
      sys.exit(1)
      
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
      logging.debug("Genemania queries not mapped: " + str(len(no_mapping)))
      warning = []
      for each in no_mapping_prots:
        warning.append(each + "(" + no_mapping_prots[each] + ") ")
        merged_out_dict[each]['Comment'] += "Genemania- not mapped;"
      logging.debug("WARNING - Dropping queries: " + ','.join(warning))
    
    if len(no_interactions) != 0:
      logging.debug("Genemania queries with no interactions: " + str(len(no_interactions)))
      warning = []
      for each in no_interactions_prots:
        warning.append(each + "(" + no_interactions_prots[each] + ") ")
        merged_out_dict[each]['Comment'] += "Genemania- no interaction;"
      logging.debug("WARNING - Dropping queries: " + ','.join(warning))
  
  return(genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict)

def categorize_gene(unique_vertex, mapping, uniprot_query):
    #primary gene, synonym gene, external synonym, external interactor
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
        
        # To be changed, drop any secondary gene or external synonym and change below statement
        # Query = node; for primary gene  & query = ""; for external interactor
        count_of_occurences = list(mapping.values()).count(node)
        query = list(mapping.keys())[list(mapping.values()).index(node)]
    
      else:
        # If gene name not in string mapping table -> external interactor
        node_category = "External Interactor"
        categories.update({node:node_category})

    return(categories)

def get_search_dicts(interaction, categories, logging, cy_debug, uniprot_query, mapping, merged_out_dict, search):
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
      logging.debug(str(search) + " dropped interaction category query nodes: " + str(len(filtered_dropped_nodes)))
      warning = []
      for each in filtered_dropped_nodes_prot:
        warning.append(each + "(" + filtered_dropped_nodes_prot[each] + "->" + get_dropped_query_primgenes[filtered_dropped_nodes_prot[each]] + ") ")
        merged_out_dict[each]['Comment'] += str(search) + "- dropped interaction category;"
      logging.debug("WARNING - Dropping queries: " + ','.join(warning))
    logging.debug(str(search) + " Total query nodes: " + str(len(primary_nodes)))
    logging.debug(str(search) + " Total interactions: " + str(count_retained_interactions))
    
  return(search_dict, merged_out_dict)  
  
def get_everything_together(each,uniprot_query, uniprot_list, max_FC_len, each_category, type, site_dict, ambigious_sites, ambigious_genes):
  '''
  For genes in merged interaction list, append other info: FC, pval, category etc
  '''  
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
          if each_pval and float(each_pval) < 0.05:
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
          if each_pval and float(each_pval) < 0.05:
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
          uniprot_list.update({term_FC:[0.0], term_pval:['NA']})
        uniprot_list.update({'query':['NA'],'significant':['NA'],'FC_exists':[0]})
        not_in_list +=1
      else:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list[term_FC].append(0.0)
          uniprot_list[term_pval].append('NA')
        uniprot_list['query'].append('NA')
        uniprot_list['significant'].append('NA')
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
          uniprot_list['length'].append(len(each))
          
        if not_in_list == 0:
          count_FC_uniprot = 1
          for each_FC,each_pval in zip(site_dict[prot_val][each_site][0],site_dict[prot_val][each_site][1]):
            if each_pval and float(each_pval) < 0.05:
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
            if each_pval and float(each_pval) < 0.05:
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
          uniprot_list.update({term_FC:[0.0], term_pval:['NA']})
        uniprot_list.update({'query':['NA'],'significant':['NA'],'FC_exists':[0],'site':['NA'],'ambigious_site':['NA']})
        not_in_list +=1
      else:
        for i in range(1,max_FC_len+1):
          term_FC = 'FC' + str(i)
          term_pval = 'pval' + str(i)
          uniprot_list[term_FC].append(0.0)
          uniprot_list[term_pval].append('NA')
        uniprot_list['query'].append('NA')
        uniprot_list['significant'].append('NA')
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
      
def get_merged_interactions(filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, type):
  for each_node1 in filtered_dict:
    for each_node2 in filtered_dict[each_node1]:
      interaction = each_node1 + " " + each_node2
      backwards = each_node2 + " " + each_node1
      if (interaction.lower() not in unique_merged_interactions) and (backwards.lower() not in unique_merged_interactions):
        unique_merged_interactions.append(interaction.lower())
      
      if each_node2.lower() not in unique_nodes:
        unique_nodes.append(each_node2.lower())
        
    if each_node1.lower() not in unique_nodes:
      unique_nodes.append(each_node1.lower())
      
  return(unique_merged_interactions, unique_nodes)        

def remove_list_duplicates(list):
  dupe_gene_list = []
  set_list = []
  for x in list:
    if x.lower() in [set_name.lower() for set_name in set_list]:
      if x.lower() not in [dupe_name.lower() for dupe_name in dupe_gene_list]:
        dupe_gene_list.append(x)
    else:
      set_list.append(x)
  return(dupe_gene_list)

def cluego_filtering(unique_nodes, cluego_mapping_file, uniprot_query, cy_debug, logging, merged_out_dict, unique_each_primgene_list, cy_session, cy_out, cy_cluego_out):
  # Categorize nodes as primary or other
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
          eprint("Error: Required columns in ClueGO file: SymbolID and UniqueID#EntrezGeneID")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
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
      merged_out_dict[each]['Comment'] += "Cluego- not found;"
    if cy_debug and (warning1 or warning2):
      logging.debug("ClueGO query + EI unavailable: " + str(len(not_found_query)) + " + " + str(len(not_found_ei)))
      if warning1:
        logging.warning("WARNING - Dropped queries: " + ','.join(warning1))
      if warning2:
        logging.warning("WARNING - Dropped EI: " + ','.join(warning2))
 
  # Drop queries with duplicate preferred gene list
  for each_acceptable in acceptable_genes_list:
    each_preferred_gene.append(acceptable_genes_list[each_acceptable]['GenePreferredName'])
  
  get_dupe_query_primgenes = {}
  '''
  if len(each_preferred_gene) != len(set(each_preferred_gene)):
    dupe_gene_list = remove_list_duplicates(each_preferred_gene)
    unique_each_preferred_list = [x for x in each_preferred_gene if x.lower() not in [name.lower() for name in dupe_gene_list]]
    get_dupe_query_primgenes = {}
    for each_dupe in dupe_gene_list:
      for k,v in acceptable_genes_list.items():
        if v['GenePreferredName'] == each_dupe:
          get_dupe_query_primgenes.update({k:v['GenePreferredName']})
    get_dupe_query_proids = get_query_from_list(uniprot_query,list(get_dupe_query_primgenes.keys()))
    
    sorted_x = sorted(get_dupe_query_proids.items(), key=lambda kv: kv[1])
    get_dupe_query_proids = collections.OrderedDict(sorted_x)    
    warning = []
    for each_in_list in get_dupe_query_proids:
      merged_out_dict[each_in_list]['Comment'] += "ClueGO- duplicate preferred gene;"
      warning.append(each_in_list + "(" + get_dupe_query_proids[each_in_list] + "->" + get_dupe_query_primgenes[get_dupe_query_proids[each_in_list]] + ") ")
    #Not found
    if cy_debug:
      logging.debug("ClueGO duplicate preferred gene: " + str(len(get_dupe_query_proids)))
      logging.warning("WARNING - Dropping queries: " + ','.join(warning))
  '''
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
    merged_out_dict[each_in_list]['Comment'] += "ClueGO non-primary query"
    warning.append(each_in_list + "(" + get_drop_query_proids[each_in_list] + "->" + drop_non_primary[get_drop_query_proids[each_in_list]] + ") ")
    
  if cy_debug:
    if len(not_found_query) != 0 or len(not_found_ei) != 0:
      logging.debug("ClueGO non-primary query + EI genes: " + str(len(not_found_query)) + " + " + str(len(not_found_ei)))
      if warning:
        logging.warning("WARNING - Dropping queries: " + ','.join(warning))
      if not_found_ei:
        logging.warning("WARNING - Dropping EI: " + ','.join(not_found_ei))
      #logging.warning("WARNING - Dropping queries: " + ','.join(get_drop_query_proids) + " with non-primary cluego genes " + ','.join(drop_non_primary))  
  
  filtered_preferred_unique_list = [acceptable_genes_list[i]['GenePreferredName'] for i in filtered_unique_list if i.lower() not in [x.lower() for x in drop_non_primary]]
  
  return(filtered_preferred_unique_list,merged_out_dict)

SEP = "/"
PORT_NUMBER = "1234"
HOST_ADDRESS = "localhost"
HEADERS = {'Content-Type': 'application/json'}

CYTOSCAPE_BASE_URL = "http://"+HOST_ADDRESS+":"+PORT_NUMBER+SEP+"v1"
CLUEGO_BASE_URL = CYTOSCAPE_BASE_URL+SEP+"apps"+SEP+"cluego"+SEP+"cluego-manager"

def writeLines(lines,out_file):
    file = open(out_file,'w')
    for line in lines:
        file.write(line)
    file.close()

def writeBin(raw,out_file):
    file = open(out_file,'wb')
    file.write(raw)
    file.close()
    
def cluego_run(organism_name,output_cluego,merged_vertex,group,select_terms, leading_term_selection, reference_file,cluego_pval):
  '''
  Run Cluego for interactor gene list
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
  response = requests.put(CLUEGO_BASE_URL+SEP+"organisms"+SEP+"set-organism"+SEP+str(organism_name), headers=HEADERS)
  
  ## Use custom reference file
  if reference_file:
    response = requests.put(CLUEGO_BASE_URL+SEP+"stats/Enrichment%2FDepletion%20(Two-sided%20hypergeometric%20test)/Bonferroni%20step%20down/false/false/true/"+SEP+reference_file)
  
  # Set the number of Clusters
  number = 1
  max_input_panel_number = number
  response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"max-input-panel"+SEP+str(max_input_panel_number))
  
  cluster = 1
  gene_list = json.dumps(merged_vertex)
  response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"upload-ids-list"+SEP+quote(str(cluster)), data=gene_list, headers=HEADERS)
  
  # 2.5 Set analysis properties for a Cluster
  input_panel_index = 1
  node_shape = "Ellipse" # ("Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle","V")
  cluster_color = "#ff0000" # The color in hex, e.g. #F3A455
  
  no_restrictions = False # "True" for no restricions in number and percentage per term
  response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"set-analysis-properties"+SEP+str(input_panel_index)+SEP+node_shape+SEP+quote(cluster_color)+SEP+str(min_number_of_genes_per_term)+SEP+str(min_percentage_of_genes_mapped)+SEP+str(no_restrictions), headers=HEADERS)
  
  #Select visual style
  visual_style = "ShowClusterDifference"
  response = requests.put(CLUEGO_BASE_URL+SEP+"cluster"+SEP+"select-visual-style"+SEP+visual_style, headers=HEADERS)
  
  ## 3.1 Get all available Ontologies
  response = requests.get(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"get-ontology-info", headers=HEADERS)
  ontology_info = response.json()
 
  i = 0
  list_ontology = []
  for each_ontology in ontology_info:
    get_each = ontology_info[each_ontology]
    m = re.search(r"name=([A-Za-z\-]+)\,",get_each)
    ontologies = m.group(1)
    if select_terms.lower() == "biological process" or select_terms.lower() == "all":
      if "BiologicalProcess" in ontologies:
        list_ontology.append(str(i)+";"+"Ellipse")
    if select_terms.lower() == "subcellular location" or select_terms.lower() == "all":
      if "CellularComponent" in ontologies:
        list_ontology.append(str(i)+";"+"Ellipse")
    if select_terms.lower() == "molecular function" or select_terms.lower() == "all":
      if "MolecularFunction" in ontologies:
        list_ontology.append(str(i)+";"+"Ellipse")
    if select_terms.lower() == "pathways" or select_terms.lower() == "all":
      if "Human-diseases" in ontologies or "KEGG" in ontologies or "Pathways" in ontologies or "WikiPathways" in ontologies:
        list_ontology.append(str(i)+";"+"Ellipse")
    i+=1

  ####Select Ontologies
  #selected_ontologies = json.dumps(["0;Ellipse","1;Triangle","2;Rectangle","3;Ellipse","4;Triangle","5;Rectangle","6;Ellipse","7;Triangle","8;Rectangle","9;Ellipse","10;Triangle","11;Rectangle","12;Ellipse"]) # (run "3.1 Get all available Ontologies" to get all options)
  selected_ontologies = json.dumps(list_ontology)
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-ontologies", data=selected_ontologies, headers=HEADERS)
  
  ## 3.1 Set kappa Score
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-kappa-score-level"+SEP+str(kappa))
  
  ## 3.2 Select Evidence Codes
  evidence_codes = json.dumps(["All"]) # (run "3.3 Get all available Evidence Codes" to get all options)
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-evidence-codes", data=evidence_codes, headers=HEADERS)

  ## 3.3 Get all available Evidence Codes
  response = requests.get(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"get-evidence-code-info", headers=HEADERS)
  
  ## 3.4 Use GO Fusion
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"true", headers=HEADERS)
  
  ## 3.5 Set min/max level for GO
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies"+SEP+"set-min-max-levels"+SEP+str(min_go)+SEP+str(max_go)+SEP+"false", headers=HEADERS)
  
  ## 3.6 GO Grouping
  # Highest%20Significance, %23Genes%20%2F%20Term, %25Genes%20%2F%20Term, %25Genes%20%2F%20Term%20vs%20Cluster
  response = requests.put(CLUEGO_BASE_URL+SEP+"grouping"+SEP+"true"+SEP+"Random"+SEP+leading_term_selection+SEP+"Kappa%20Score"+SEP+str(1)+SEP+str(30)+SEP+str(30), headers=HEADERS)
  
  ## 3.7 Set Pval
  response = requests.put(CLUEGO_BASE_URL+SEP+"ontologies/true"+SEP+str(cluego_pval))
  
  #### Run ClueGO Analysis ####
  # Run the analysis an save log file
  analysis_name = "ClueGO Network"
  selection = "Continue analysis"                                                                         
  response = requests.get(CLUEGO_BASE_URL+SEP+analysis_name+SEP+selection, headers=HEADERS)
  #log_file_name = "ClueGO-log.txt"
  #writeLines(response.text,log_file_name)

  # Get network id (SUID) (CyRest function from Cytoscape)
  response = requests.get(CYTOSCAPE_BASE_URL+SEP+"networks"+SEP+"currentNetwork", headers=HEADERS)
  current_network_suid = response.json()['data']['networkSUID']
  
  # 4.2 Get ClueGO result table
  if output_cluego:
    response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-cluego-table"+SEP+str(current_network_suid))
    table_file_name = output_cluego
    writeLines(response.text,table_file_name)
  
def get_interactions_dict(filtered_dict, search, merged_out_dict):
  lower_filtered = [name.lower() for name in filtered_dict]
  for each_uniprot_query in merged_out_dict:    
    for name in filtered_dict:
      if name.lower() == merged_out_dict[each_uniprot_query]['Primary'].lower():
        interaction_list = filtered_dict[name]
        merged_out_dict[each_uniprot_query][search] += ';'.join(interaction_list)
        break
        
  return(merged_out_dict)

def write_into_out(merged_out_dict, out):
  if not out:
    return
  with open(out,'a') as csv_file:
    for each_prot in merged_out_dict:
      line = each_prot + "," + merged_out_dict[each_prot]['Primary'] + "," + merged_out_dict[each_prot]['String'] + "," + merged_out_dict[each_prot]['Genemania'] + "," + merged_out_dict[each_prot]['Comment'] + "\n"
      csv_file.write(line)
  csv_file.close()
    
def cluego_input_file(cluego_inp_file, cy_debug, logging, cy_session, cy_out, cy_cluego_out):
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
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit(1)
      
      else:
        if row:  
          if "[" in row[genes] and "]" in row[genes]:
            each_gene_list = (row[genes])[1:-1]
          if ", " in row[genes]:
            each_gene_list = each_gene_list.split(', ')
        
          for each in each_gene_list:
            if each not in unique_gene:
              unique_gene.append(each)
            
          if row[goterm] not in top_annotations:
            top_annotations.update({row[goterm]:each_gene_list})
          else:
            eprint("Error: Duplicate GOID found in cluego input file")
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
            sys.exit(1)
          
      line_count+=1
  return(top_annotations, unique_gene)

def cy_sites_interactors_style(merged_vertex, merged_interactions, uniprot_list, max_FC_len, each_category, pval_style, type):
  '''
  Styling + visualization for the entire gene list interaction network & its sites
  '''

  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  G = igraph.Graph()
  
  site_interactions = []
  site_length = []
  sig_na = []
  fc_na = []
  fc_val = []
  pval_val = []
  val_fc_gene = {}
  get_each_site_name = []
  include_gene_data = []
  fc_merged_vertex = {}
  query_val = []
  for each_site_gene,each_ambi_site in zip(uniprot_list['site'],uniprot_list['ambigious_site']):
    each_gene = (each_site_gene.split("-"))[1]
    each_site = (each_site_gene.split("-"))[0]
    if each_gene.lower() in merged_vertex:
      get_each_site_name.append(each_ambi_site)
      query_val.append("Site")
      G.add_vertex(each_site_gene)
      interaction = each_site_gene + " " + each_gene
      site_interactions.append(interaction)
      site_length.append(len(each_site))        
  for each_vertex in merged_vertex:
    if each_vertex in uniprot_list['name']:
      indexOf = uniprot_list['name'].index(each_vertex.lower())
      include_gene_data.append(uniprot_list['ambigious_genes'][indexOf])      
      fc_na.append(-1)
      fc_val.append(0.0)
      pval_val.append(-1)
    else:
      fc_na.append(0)    
    query_val.append("Gene")
    G.add_vertex(each_vertex)  
    sig_na.append(0)
    
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
  G.vs["name"] = get_each_site_name + include_gene_data
  G.vs["length"] = site_length + uniprot_list['length']
  
  G.vs["significant"] = uniprot_list["significant"] + sig_na
  G.vs["FC_exists"] = uniprot_list["FC_exists"] + fc_na
  
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    G.vs[term_FC] = uniprot_list[term_FC] + fc_val
    G.vs[term_pval] = uniprot_list[term_pval] + pval_val
  
  G.vs["query"] = query_val
  
  degree = G.degree()
  G.vs["degree"] = degree
  
  cy = CyRestClient()
  g_cy = cy.network.create_from_igraph(G, name="Interaction Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Initial Network Style')
  basic_settings = {
    'NODE_CUSTOMGRAPHICS_1':"org.cytoscape.BarChart",
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999",
    'NODE_SIZE':"30",
    'EDGE_WIDTH':"2"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  degree_to_label_size = StyleUtil.create_slope(min=min(degree),max=max(degree),values=(10, 20))
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
  my_style.create_continuous_mapping(column='length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  columns = ['length','degree']

  chart={'cy_range':'[2.0,9.0]','cy_colorScheme':"Custom"}
  chart = json.dumps(chart)
  loaded_r = json.loads(chart)
  just_color = 0
  pval_sig = 0
  slash_delim = chr(92)
  
  shape_kv_pair = {
    "Function":"OCTAGON",
    "Site":"ELLIPSE",
    "Gene":"ROUND_RECTANGLE"
  }
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)
  if type == "5":
    label_color_kv_pair = {
        "Function":"#FFFFFF",
        "Gene":"#FFFFFF",
        "Site":"#000000"
    }
  else:
    label_color_kv_pair = {
        "Function":"#FFFFFF",
        "Gene":"#000000",
        "Site":"#000000"
    }
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)
  
  width_kv_pair = {
    "Function":"210",
    "Gene":"70",
    "Site":"60"
  }
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,uniprot_list["FC_exists"] + fc_na,query_val,"-1",get_each_site_name + include_gene_data,max_FC_len, uniprot_list)
    if pval_style:
      pval_sig = 1
    
  elif max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list)
    if pval_style:
      pval_sig = 1
   
  if pval_sig == 1:
    my_style = pval(my_style)
  
  cy.style.apply(my_style, g_cy)
  return(fc_merged_vertex)
  
def cy_interactors_style(merged_vertex, merged_interactions, uniprot_list, max_FC_len, each_category, pval_style):
  '''
  Styling + visualization for the entire gene list interaction network
  '''

  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  G = igraph.Graph()
  
  for each in merged_vertex:
    G.add_vertex(each)
  
  count_each = 0
  for each in merged_interactions:
    each_interaction_name = each.split(" ")
    G.add_edge(each_interaction_name[0],each_interaction_name[1],name=each,interaction="pp")
    count_each +=1

  G.vs
 
  G.vs["name"] = uniprot_list['name']
  G.vs["shared name"] = uniprot_list['ambigious_genes']
  G.vs["length"] = uniprot_list['length']
  
  if max_FC_len >= 1:
    G.vs["significant"] = uniprot_list["significant"]
    G.vs["FC_exists"] = uniprot_list["FC_exists"]
  
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    G.vs[term_FC] = uniprot_list[term_FC]
    G.vs[term_pval] = uniprot_list[term_pval]
  
  G.vs["query"] = uniprot_list["query"]
  
  category_present = 0
  
  for each in each_category:
    category_present = 1
    G.vs[each] = uniprot_list[each]
  if category_present == 1:
    G.vs['category_true'] = uniprot_list['category_true']

  degree = G.degree()
  G.vs["degree"] = degree
  
  cy = CyRestClient()
  g_cy = cy.network.create_from_igraph(G, name="Interaction Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Initial Network Style')
  basic_settings = {
    'NODE_CUSTOMGRAPHICS_1':"org.cytoscape.BarChart",
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999",
    'NODE_SIZE':"30",
    'EDGE_WIDTH':"2"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  degree_to_label_size = StyleUtil.create_slope(min=min(degree),max=max(degree),values=(10, 20))
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
  my_style.create_continuous_mapping(column='length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  columns = ['length','degree']

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
    my_style.create_discrete_mapping(column='category_true', col_type='Double', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,uniprot_list['FC_exists'],uniprot_list["query"],"1",uniprot_list['name'],max_FC_len, uniprot_list)
    if pval_style:
      pval_sig = 1
    
  elif (max_FC_len == 0 and category_present == 0):
    just_color = 1
    
  elif max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list)
    if pval_style:
      pval_sig = 1
  
  if just_color == 1:
    my_style = color(my_style, uniprot_list)
  
  if pval_sig == 1:
    my_style = pval(my_style)
  
  cy.style.apply(my_style, g_cy)

def cy_category_style(merged_vertex, merged_interactions, uniprot_list, each_category):
  '''
  Styling + visualization for the cateory network
  '''
  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  G = igraph.Graph()
  
  additional_nodes = []
  additional_len = []
  additional_query = []
  additional_each = []
  additional_category_true = []
  length_of = []
  breadth_of = []
  
  for each in merged_vertex:
    G.add_vertex(each)
    additional_query.append("Gene")
    val_length_of_val = len(each) #* 15
    length_of.append(val_length_of_val)
    val_breadth_of_val = 30
    breadth_of.append(val_breadth_of_val)
    
  for each in each_category:
    G.add_vertex(each)
    additional_nodes.append(each)
    val_length_of = len(each)
    length_of.append(val_length_of)
    if val_length_of > 18:
      val_breadth_of_val = val_length_of/18.0*50.0
    else:
      val_breadth_of_val = 50
    breadth_of.append(val_breadth_of_val)
   
    additional_query.append("Function")
    additional_each.append(0.0)
    additional_category_true.append(0.0)
  
  for each_node in merged_vertex:
    for each_cat in each_category:
      index = uniprot_list['name'].index(each_node)
      if uniprot_list[each_cat][index] == 1.0:
        G.add_edge(each_node,each_cat,name=each_node + " " + each_cat,interaction="with")
  G.vs
  G.vs["name"] = uniprot_list['name'] + additional_nodes
  G.vs["length"] = length_of
  G.vs["breadth"] = breadth_of
  G.vs["query"] = additional_query
  
  for each in each_category:
    G.vs[each] = uniprot_list[each] + additional_each
    
  G.vs['category_true'] = uniprot_list['category_true'] + additional_category_true

  degree = G.degree()
  G.vs["degree"] = degree
  
  cy = CyRestClient()
  g_cy = cy.network.create_from_igraph(G, name="Category Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  my_style = cy.style.create('Category-Network-Style')
  basic_settings = {
    #'NODE_FILL_COLOR':"#33CCFF",
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_WIDTH':"2"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  
  shape_kv_pair = {
    "Function":"OCTAGON",
    "Gene":"ROUND_RECTANGLE",
  }
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)
  height_kv_pair = {
    "Function":"70",
    "Gene":"30",
  }
  width_kv_pair = {
    "Function":"210",
    "Gene":"60",
  }
  label_color_kv_pair = {
    "Function":"#FFFFFF",
    "Gene":"#000000",
  }
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
  my_style.create_continuous_mapping(column='length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_HEIGHT', mappings=height_kv_pair)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)

  my_style = get_category(my_style,uniprot_list['category_true']+additional_category_true,"2",additional_query,uniprot_list['name']+additional_nodes,each_category)
  
  cy.style.apply(my_style, g_cy)
  data = [
    {
    "visualPropertyDependency": "nodeSizeLocked",
    "enabled": False
    }
  ]
  response = requests.put("http://localhost:1234/v1/styles/Category-Network-Style/dependencies", json=data)
  
def singleFC(my_style, uniprot_list):
  '''
  Styling specific to singleFC case
  '''
  basic_settings = {
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999"
  }
  if min(uniprot_list['FC1']) != 0 and max(uniprot_list['FC1']) != 0:
    points1 =  [
      {
        "value": min(uniprot_list['FC1']),
        "lesser": "#0000CC" ,
        "equal": "#3399FF", 
        "greater": "#3399FF"
      },
      {
        "value": 0,
        "lesser": "#E2E2E2",
        "equal": "#E2E2E2",
        "greater": "#E2E2E2"
      },
      {
        "value": max(uniprot_list['FC1']),
        "lesser": "#FF0000",
        "equal": "#FF0000", 
        "greater": "#FF6633" 
      }
    ]
  elif min(uniprot_list['FC1']) == 0:
    points1 =  [
      {
        "value": 0,
        "lesser": "#E2E2E2",
        "equal": "#E2E2E2",
        "greater": "#E2E2E2"
      },
      {
        "value": max(uniprot_list['FC1']),
        "lesser": "#3399FF",
        "equal": "#3399FF",
        "greater": "#0000CC"
      }
    ]
  else:
    points1 =  [
      {
        "value": min(uniprot_list['FC1']),
        "lesser": "#FF6633",
        "equal": "#FF0000",
        "greater": "#FF0000"
      },
      {
        "value": 0,
        "lesser": "#E2E2E2",
        "equal": "#E2E2E2",
        "greater": "#E2E2E2"
      }
    ]
  my_style.create_continuous_mapping(column='FC1',vp='NODE_FILL_COLOR',col_type='Double',points=points1)
  return(my_style)
  
def multipleFC(my_style,FC_exists,query,func,name,max_FC_len,uniprot_list):
  '''
  Styling specific to multipleFC case
  '''
  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  bar_columns = ""
  color_columns = ""
  min_fc = 0
  max_fc = 0
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    if i == 1:
      min_fc = min(uniprot_list[term_FC])
      max_fc = max(uniprot_list[term_FC])
    else:
      if min(uniprot_list[term_FC]) < min_fc:
        min_fc = min(uniprot_list[term_FC])
      if max(uniprot_list[term_FC]) > max_fc:
        max_fc = max(uniprot_list[term_FC])
    bar_columns = bar_columns + '\"' + term_FC + '\"' + ","
    color_columns = color_columns + '\"' + color_code[i-1] + '\"' + ","
  bar_columns = bar_columns[:-1]
  color_columns = color_columns[:-1]
 
  value = "org.cytoscape.BarChart:{" +'\"' + "cy_range" + '\"' + ":[" + str(min_fc) + "," + str(max_fc) + "]," + '\"' + "cy_globalRange" + '\"' + ":true," + '\"' + "cy_colors" + '\"' + ":[" + color_columns + "]," + '\"' + "cy_dataColumns" + '\"' + ":[" + bar_columns + "]}"

  kv_pair = {
    "1":value
  }
  kv_pair_node_position = {
    "1":"S,N,c,0.00,0.00"
  }
  my_style.create_discrete_mapping(column='FC_exists', col_type='Double', vp='NODE_CUSTOMGRAPHICS_1',mappings=kv_pair)
  my_style.create_discrete_mapping(column='FC_exists', col_type='Double', vp='NODE_LABEL_POSITION',mappings=kv_pair_node_position)

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
    my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_FILL_COLOR',mappings=kv_pair_color)
  
  else:
    for i,does_fc_exist,each_name in zip(query,FC_exists,name):
      if i == "Function":
        color_value = "#333333"
        width_value = "1"
      elif i == "Gene" and does_fc_exist == 1.0:
        color_value = "#FFFFFF"
        width_value = "1"
      elif i == "Gene" and does_fc_exist != 1.0:
        color_value = "#FFFF99"
        width_value = "1"
      elif i == "Site":
        color_value = "#FFFFFF"
        width_value = "1"
      else:
        color_value = "#E2E2E2"
        width_value = "1"
        
      if each_name not in kv_pair_color:
        kv_pair_color.update({each_name:color_value})
        width_kv_pair.update({each_name:width_value})
        
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_FILL_COLOR',mappings=kv_pair_color)
    my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_BORDER_WIDTH', mappings=width_kv_pair)
  return(my_style)
  
def get_category(my_style,is_category_present,cat_val,query,name,each_category):
  '''
  Styling specific to category case
  '''
  basic_settings = {
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999"
  }
  my_style.update_defaults(basic_settings)
  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 

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
    "1.0":"2",
    "0.0":"2"
  }
  
  my_style.create_discrete_mapping(column='category_true', col_type='Double', vp='NODE_CUSTOMGRAPHICS_1',mappings=kv_pair)
  my_style.create_discrete_mapping(column='category_true', col_type='Double', vp='NODE_LABEL_POSITION',mappings=kv_pair_node_position)
  my_style.create_discrete_mapping(column='category_true', col_type='Double', vp='EDGE_WIDTH',mappings=kv_edge_width)

  if cat_val == "2":
    kv_pair_color = {}
    width_kv_pair = {}
    node_width = {}
    for i,each_is_category_present,each_name in zip(query,is_category_present,name):
        if i == "Function":
          color_value = "#333333"
          width_value = "1"
          node_width_value = "210"
        elif i == "Gene" and each_is_category_present == 1.0:
          color_value = "#FFFFFF"
          width_value = "0"
          node_width_value = "30"
        else:
          color_value = "#E2E2E2"
          width_value = "1"
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
  '''
  If pval significant, border node with red color
  '''
  #"SansSerif.plain,plain,12"

  width_kv_pair = {
    "1":"2",
  }
  bc_kv_pair = {
    "1":"#0000FF",
  }
  my_style.create_discrete_mapping(column='significant', col_type='String', vp='NODE_BORDER_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='significant', col_type='String', vp='NODE_BORDER_PAINT', mappings=bc_kv_pair)
  return(my_style)
  
def color(my_style, uniprot_list):
  '''
  Styling specific to list only (noFC) case
  '''
  color_kv_pair = {}
  for each_query in uniprot_list["query"]:
    if each_query == "NA":
      if each_query not in color_kv_pair:
        color_kv_pair.update({each_query:"#E2E2E2"})
    else:
      color_kv_pair.update({each_query:"#FFFFBF"})
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  return(my_style)

def cy_pathways_style(cluster, each_category, max_FC_len, pval_style, uniprot_list, type, fc_merged_vertex):
  '''
  Based on top clusters picked, construct function interaction network + visualization and styling
  '''
  G = igraph.Graph()
  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  
  cluster_list = cluster
  
  name = []
  function_only = cluster_list.keys()
  query = []
  length_of = []
  breadth_of = []
   
  merged_vertex = []
  merged_vertex_sites_only = []
  query_val_noFC = []
  count_each = 0
  all_interactions = []
  function_fc_val = {}

  category_present = 0
  for each in each_category:
    category_present = 1
    break
      
  for each in cluster_list:
    if each not in merged_vertex:
      G.add_vertex(each)
      merged_vertex.append(each)
      merged_vertex_sites_only.append(each)
      if each in function_only:
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
        if each_gene in function_only:
          merged_vertex_sites_only.append(each_gene)
          query.append('Function')
          if max_FC_len == 0 and not category_present:
            query_val_noFC.append('Function')
          val_length_of = len(each_gene)
         
          length_of.append(val_length_of)
          if val_length_of > 18:
            val_breadth_of_val = val_length_of/18.0*50.0
          else:
            val_breadth_of_val = 50
          breadth_of.append(val_breadth_of_val)
         
        else:
          query.append('Gene')
          indexOf = uniprot_list['name'].index(each_gene.lower())
          merged_vertex_sites_only.append(uniprot_list['ambigious_genes'][indexOf])
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
              G.add_vertex(uniprot_list['site'][each_index])
              merged_vertex.append(uniprot_list['site'][each_index])
              merged_vertex_sites_only.append(uniprot_list['ambigious_site'][each_index])
              query.append('Site')
              val_length_of_val = len(each_gene) #* 15
              length_of.append(val_length_of_val)
              val_breadth_of_val = 30
              breadth_of.append(val_breadth_of_val)
              name_edge = uniprot_list['site'][each_index] + " with " + each_gene
              G.add_edge(uniprot_list['site'][each_index],each_gene,name=name_edge)
              all_interactions.append(name_edge)
              
      if max_FC_len == 1 and not (type == "5" or type == "6"):
        if each_gene.lower() in uniprot_list["name"]:
          index = uniprot_list["name"].index(each_gene.lower())
          FC_val_each_gene = (uniprot_list['FC1'])[index]
          
      name_edge = each + " with " + each_gene
      G.add_edge(each,each_gene,name=name_edge)
      all_interactions.append(name_edge)
      
  G.vs
  G.vs["query"] = query
  G.vs["name"] = merged_vertex
  G.vs["shared name"] = merged_vertex_sites_only
  G.vs["length"] = length_of
  G.vs["breadth"] = breadth_of
  degree = G.degree()
  G.vs["degree"] = degree
  
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
    
  if category_present == 1:
    G.vs["category_true"] = is_category_present  
  
  add_term_pval = []
  FC_exists = []
  fc_exists_count = 0

  if type == "5" or type == "6":
    search_list_for_fc = [a.lower() for a in uniprot_list["site"]]
  else:
    search_list_for_fc = uniprot_list["name"]
  for i in range(1,max_FC_len+1):
    k = 0
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
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
          else:
            FC_exists.append(0.0)
        else:
          if ((uniprot_list[term_FC])[index]) != 0:
            FC_exists[k] = 1.0
              
      else:
        if max_FC_len == 1 and not (type == "5" or type == "6"):
          add_term_FC.append(100.0)
       
        else:
          add_term_FC.append(0.0)
          
        add_term_pval.append("NA")
        significant_val.append("NA")
        
        if fc_exists_count == 0:
          FC_exists.append(0.0)    
             
      k+=1
    fc_exists_count+=1

    G.vs[term_FC] = add_term_FC
    G.vs[term_pval] = add_term_pval
    G.vs["significant"] = significant_val
    G.vs["FC_exists"] = FC_exists
    
  if max_FC_len == 0 and not category_present:
    G.vs["query_val"] = query_val_noFC
  
  cy = CyRestClient()

  g_cy = cy.network.create_from_igraph(G, name="Ontology Network")
  
  cy.layout.apply(name='cose', network=g_cy)
  
  my_style = cy.style.create('GAL_Style3')
  
  basic_settings = {
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_WIDTH':"2"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  
  shape_kv_pair = {
    "Function":"OCTAGON",
    "Gene":"ROUND_RECTANGLE",
    "Site":"ELLIPSE"
  }
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_SHAPE', mappings=shape_kv_pair)
  height_kv_pair = {
    "Function":"70",
    "Gene":"30",
    "Site":"30"
  }
  if type == "5" or type == "6":
    width_kv_pair = {
      "Function":"210",
      "Gene":"70",
      "Site":"60"
    }
  else:
    width_kv_pair = {
      "Function":"210",
      "Gene":"60",
      "Site":"60"
    }
  if type == "5":
    label_color_kv_pair = {
      "Function":"#FFFFFF",
      "Gene":"#FFFFFF",
      "Site":"#000000"
    }
  else:
    label_color_kv_pair = {
      "Function":"#FFFFFF",
      "Gene":"#000000",
      "Site":"#000000"
    }
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
  my_style.create_continuous_mapping(column='length',vp='NODE_LABEL_FONT_SIZE',col_type='Double',points=node_label_size)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_HEIGHT', mappings=height_kv_pair)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_WIDTH', mappings=width_kv_pair)
  my_style.create_discrete_mapping(column='query', col_type='String', vp='NODE_LABEL_COLOR', mappings=label_color_kv_pair)

  pval_sig = 0
  if max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list)
    if pval_style:
      pval_sig = 1
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,FC_exists,query,"2",merged_vertex, max_FC_len, uniprot_list)

    if pval_style:
      pval_sig = 1
    
  if category_present:
    my_style = get_category(my_style,is_category_present,"2",query,merged_vertex,each_category)
  
  if max_FC_len == 0 and not category_present:
    color_kv_pair = {}
    for each_val in query_val_noFC:
      if each_val == "NA":
        color_kv_pair.update({each_val:"#E2E2E2"})
      elif each_val == "Function":
        color_kv_pair.update({each_val:"#333333"})
      else:
        color_kv_pair.update({each_val:"#FFFFBF"})

    my_style.create_discrete_mapping(column='query_val', col_type='String', vp='NODE_FILL_COLOR', mappings=color_kv_pair)
  
  if pval_sig == 1:
    my_style = pval(my_style)
 
  cy.style.apply(my_style, g_cy)
  
  data = [
    {
    "visualPropertyDependency": "nodeSizeLocked",
    "enabled": False
    }
  ]
  response = requests.put("http://localhost:1234/v1/styles/GAL_Style3/dependencies", json=data)
  
def remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out):
  if cy_debug:
    logging.handlers = []
    if path.exists("Cytoscape.log"):
      os.remove("Cytoscape.log")
  if cy_session and path.exists(cy_session):
    os.remove(cy_session)
  if cy_out and path.exists(cy_out):
    os.remove(cy_out)
  if cy_cluego_out and path.exists(cy_cluego_out):
    os.remove(cy_cluego_out)
    
def main(argv):
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
  cy_exe = None
  
  try:
    opts, args = getopt.getopt(argv, "i:s:l:t:r:m:o:ng:f:p:z:h:a:y:u:d:b:x:e:c:",["in=","species=","limit=","type=","score=","mapping=","output=","significant","grouping=","fccutoff=","pvalcutoff=","visualize=","reference-path=","input_cluego=","cluego-pval=","run=","mods=","fasta-file=","enzyme=","gui","cytoscape-executable=","cytoscape-session-file="])
    for opt, arg in opts:
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
      else:
        help = True      
  except getopt.GetoptError as e:
    help = True
    
  if not cy_in or not cy_species or not cy_type or not cy_out_dir or not cy_exe:
    help = True
    
  if help:
    print("PINE")
    print("---------------------------------------------------------------------------------------------")
    print("Usage:         cytoscape_api.py -i input.csv -o output_dir -c cluego_out.txt -t input_type -s species -m cluego_map_file.gz")
    print("Argument:      -i [--in]: input file in csv format with the following headers as applicable: ")
    print("Argument:      -o [--output]: path to output directory")     
    print("Argument:      -t [--type]: analysis type [Allowed: noFC, singleFC, multiFC, category, singlefc-ptm, multifc-ptm]")   
    print("Argument:      -s [--species]: species [Allowed: human, mouse, rat]")   
    print("Argument:      -x [--enzyme]: enzyme name [Allowed:]")   
    print("Argument:      -d [--mods]: comma separated list of modifications of interest")
    print("Argument:      -b [--fastafile]: path to fasta file")
    print("Argument:      -m [--mapping]: path to cluego mapping file compressed in .gz format") 
    print("Argument:      -e [--cytoscape-executable]: the path to the Cytoscape executable")
    print("Argument(opt): -f [--fccutoff]: fold change cutoff for input [Default: abs(FC) >= 0.0]")
    print("Argument(opt): -p [--pvalcutoff]: pvalue cutoff for input [Default: pval > 1.0]")
    print("Argument(opt): -n [--significant]: outline statistically significant nodes")    
    print("Argument(opt): -u [--run]: interaction databases [Allowed: string, genemania, both; Default: both]")
    print("Argument(opt): -r [--score]: interaction confidence score for string [Default:0.4, Range 0-1]")
    print("Argument(opt): -l [--limit]: maximum number of external interactors [Default:0, Range:0-100]")
    print("Argument(opt): -z [--visualize]: ontology type [Allowed: biological process, subcellular location, molecular function, pathways, all; Default: pathways]")
    print("Argument(opt): -g [--grouping]: network specificity indicating general, representative and specific pathways [Allowed: global, medium, detailed; Default: medium]")
    print("Argument(opt): -y [--cluegopval]: pvalue cutoff for enrichment analysis [Default: 0.05]")
    print("Argument(opt): -h [--referencepath]: path to background reference file for enrichment")
    print("Argument(opt): -a [--inputcluego]: filtered cluego file with ontology terms of interest")
    sys.exit()

  if not os.path.isdir(cy_out_dir):
    eprint("Error: output is not a directory")
    sys.exit(1)
  
  timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
  
  # create file names
  if cy_cluego_inp_file:
    cy_out = None
    cy_cluego_out = None
    path_to_cluego = (os.path.abspath(cy_cluego_inp_file))
    session_name = path_to_cluego.split(".")[0]
    cy_session = session_name + ".cys"
    logging_file = session_name + ".log"
    reanalyze_flag = True
  else:  
    path_to_new_dir = os.path.join(os.path.abspath(cy_out_dir), timestamp) + "_PINE"
    os.mkdir(path_to_new_dir) 
    cy_out = os.path.join(path_to_new_dir, "Interactions.csv")
    cy_cluego_out = os.path.join(path_to_new_dir, "PINE.cluego.txt")
    cy_session = os.path.join(path_to_new_dir, "PINE.cys")
    logging_file = os.path.join(path_to_new_dir, "PINE.log")
    reanalyze_flag = False
    if gui_mode: # notify gui of the output directory
      print(f"COMMAND FILE-SESSION {path_to_new_dir}")
      
  if os.path.exists(cy_session):
    eprint("Session " + cy_session + " already exists. Please provide a different name")
    sys.exit(1)
  
  # set up logging
  if gui_mode:
    logging = setup_logger("PINE", logging_file, with_stdout=True)
  else:
    logging = setup_logger("PINE", logging_file)

  if cy_species.lower() == "human":
    tax_id = "9606"
    organism_name = "Homo Sapiens"
  elif cy_species.lower() == "mouse":
    tax_id = "10090"
    organism_name = "Mus Musculus"
  elif cy_species.lower() == "rat":
    tax_id = "10116"
    organism_name = "Rattus norvegicus"
  else:
    eprint("Error: Species currently supported are human, mouse, rat")
    sys.exit(1)

  if not ('\\' in cy_session or "/" in cy_session):
    cwd = os.getcwd()
    cy_session = cwd + "\\" + cy_session
    
  if not ('.cys' in cy_session):
    eprint("Error: Cytoscape session file must have a valid name with .cys extension")
    sys.exit(1)
  
  if not (float(cluego_pval) >=0.0 and float(cluego_pval) <=1.0):
    eprint("Error: Cluego pvalue range must be between 0 to 1")
    sys.exit(1)
  
  allowed_runs = ["string","genemania","both"]
  if cy_run.lower() not in allowed_runs:  
    eprint("Error: String run type must be one of the following: " + ','.join(allowed_runs))
    sys.exit(1)
    
  # Highest%20Significance, %23Genes%20%2F%20Term, %25Genes%20%2F%20Term, %25Genes%20%2F%20Term%20vs%20Cluster
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
  # 1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = category
  if cy_type.lower() not in allowed_type:
    eprint("Error: Input type must be one of the following: " + (',').join(allowed_type))
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
  
  if (cy_ptm_sites and cy_type.lower() == "nofc") or (cy_ptm_sites and cy_type.lower() == "category"):
    eprint("Error: PTM visualization only supported along with singlefc or multifc analysis")
    sys.exit(1)
    
  allowed_selections = ["biological process","subcellular location","molecular function","pathways","all"]
  if select_terms.lower() not in allowed_selections:
    eprint("Error: The visualization type must be one of the following: " + (',').join(allowed_selections))
    sys.exit(1)
  
  try:
    cy_fc_cutoff = float(cy_fc_cutoff)
  except:
    eprint("Error: FC cutoff must be a number") 
    sys.exit(1)

  try:
    cy_pval_cutoff = float(cy_pval_cutoff)
  except:
    eprint("Error: PVal cutoff must be a number") 
    sys.exit(1)
  
  if not (cy_pval_cutoff >= 0.0 and cy_pval_cutoff <= 1.0):
    eprint("Error: PVal cutoff must range between 0 and 1")
    sys.exit(1)
    
  allowed_groups = ["global","medium","detailed"]
  if cy_cluego_grouping.lower() not in allowed_groups:
    eprint("Error: Cluego grouping must be one of the following:" + (',').join(allowed_groups))
    sys.exit(1)
    
  cy_score = float(cy_score)*1000
  if not (cy_score >= 0.0 and cy_score <= 1000.0):
    eprint("Error: Invalid string score provided. Value must be between 0 to 1; Confidence levels for string score- Low = 0.150, Medium = 0.400, High = 0.700, Highest = 0.900")
    sys.exit(1)
  
  # Rounding off score - ex: 500.9 and above = 501, less = 500
  cy_score_ceil = math.ceil(cy_score)
  if((cy_score_ceil-cy_score) <= 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = cy_score_ceil
  elif ((cy_score_ceil-cy_score) > 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = math.floor(cy_score)
  
  if not (int(cy_lim) >= 0 and int(cy_lim) <=100):
    eprint("Error: Limit on additional interactors is 100. Please choose a number between 0 and 100")
    sys.exit(1)
    
  try:
    if cy_debug:
      logging.debug("Starting PINE Analysis...\n\n") 
    try:
      r = requests.get("http://localhost:1234/v1/version")
      path_to_docs = os.path.expanduser("~\Documents")
      session_filename = os.path.join(path_to_docs, timestamp + "-exited-session.cys") # save current session
      r = requests.get("http://localhost:1234/v1/commands/session/save%20as?file=" + urllib.parse.quote(session_filename, safe=""))
      r = requests.get("http://localhost:1234/v1/commands/command/quit")
      logging.warning("Cytoscape was already open with an existing session.  Saved existing session to: " + session_filename)
      wait_counter = 0
      while wait_counter < 120: # give 2 minutes to exit cytoscape
        try:
          requests.get("http://localhost:1234/v1/version")
        except:
          break
        time.sleep(5)
        wait_counter += 5
    except:
      pass

    # open cytoscape
    subprocess.Popen([cy_exe])
    wait_counter = 0
    while wait_counter < 120: # give 2 minutes max for cytoscape to open
      try:
        r = requests.get("http://localhost:1234/v1/version")
        test = r.json()
        break
      except:
        time.sleep(5)
        wait_counter += 5

    # Check Cytoscape version
    request = requests.get('http://localhost:1234/v1/version')
    cy_version = request.json()
    if not bool(re.match('^3.7', cy_version['cytoscapeVersion'])):
      eprint("Error: Cytoscape version must be 3.7.0 and above")
      sys.exit(1)
    
    # Start a new session 
    requests.post('http://localhost:1234/v1/commands/session/new')
    
    # Apps installed
    request = requests.post('http://localhost:1234/v1/commands/apps/list installed')
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
    
    if not app_genemania or not ver_genemania == "3.5.1":
      eprint("Error: Cytoscape app GeneMANIA v3.5.1 not installed or not responding properly")
      sys.exit(1)
      
    if not app_cluego or not (ver_cluego == "2.5.4" or ver_cluego == "2.5.5"):
      eprint("Error: Cytoscape app ClueGO v2.5.4/v2.5.5 not installed or not responding properly")
      sys.exit(1)
    
    if cluego_reference_file and ver_cluego != "2.5.5":
      eprint("Error: Using ClueGO custom reference file needs version 2.5.5 of ClueGO")
      sys.exit(1)
    
    body = dict(offline=True)
    response = requests.post("http://localhost:1234/v1/commands/genemania/organisms", json=body)
    genemania_bool = False
    for each in response.json()['data']['organisms']:
      if tax_id == str(each['taxonomyId']):
        genemania_bool = True
    if not genemania_bool:
      eprint("Error: Please install " + cy_species + " dataset in Genemania")
      sys.exit(1)
    
    if not cy_cluego_inp_file:
      # Check if cluego path exists
      if not cy_map:
        eprint("Error: ClueGO mapping file path must be provided")
        sys.exit(1)
        
      if not path.exists(cy_map):
        eprint("Error: Path to ClueGO mapping file " + cy_map + " does not exist")
        sys.exit(1)
      else:
        if ver_cluego not in cy_map:
          eprint("Error: ClueGO version mismatch. Version installed is " + ver_cluego + " and version contained in path to ClueGO mapping file " + cy_map)
          sys.exit(1)
        if organism_name.lower() not in cy_map.lower():
          eprint("Error: Species mismatch.  Species parameter provided is " + organism_name + " and species contained in path to ClueGO mapping file is " + cy_map)
          sys.exit(1)
    
      if not (".gz" in cy_map or "gene2accession" in cy_map):
        eprint("Error: ClueGO mapping file must refer to the species gene2accession gz file")
        sys.exit(1)
    
    database_dict = {} 
    mods_list = []
    if cy_type_num == "5" or cy_type_num == "6":
      if not cy_fasta_file or not cy_mods:
        print("Error: Fasta file and List of Modifications are mandatory for site analysis")
        sys.exit(1)        
      database_dict = db_handling(cy_fasta_file)
      mods_list = cy_mods.split(",") 
       
    #Read input and obtain protid 
    if cy_debug:
      logging.debug("Step 1: Start processing the input protein list at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
    unique_each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict,initial_length, site_info_dict, ambigious_sites = preprocessing(cy_in, cy_type_num, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out, database_dict, mods_list, cy_fasta_file, cy_enzyme)
    
    # FC and Pval cutoff
    if not (cy_fc_cutoff == 0.0 and cy_pval_cutoff == 1.0):
      unique_each_protein_list, prot_list, merged_out_dict = inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict)
    
    # Limit query inpt number = 1500
    if len(unique_each_protein_list) > 1500:
      eprint("Error: The query input is too big. Currently supporting upto 1500 query protein ids")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
      sys.exit(1)
     
    #Uniprot API call to get primary gene, synonym
    if cy_debug:
      logging.debug("\nStep 2: Start the uniprot api call at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
      logging.debug("Uniprot query: " + str(len(unique_each_protein_list)))
        
    uniprot_query,each_primgene_list,merged_out_dict,ambigious_genes = uniprot_api_call(unique_each_protein_list, prot_list, cy_type_num, cy_debug, logging, merged_out_dict, organism_name, cy_session, cy_out, cy_cluego_out)
        
    if cy_cluego_inp_file:
      if cy_debug:
        logging.debug("\nStep 2.5: Start processing the input ClueGO list at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
      leading_term_cluster, each_primgene_list = cluego_input_file(cy_cluego_inp_file, cy_debug, logging, cy_session, cy_out, cy_cluego_out)
      cy_lim = 0
      if cy_debug:
        logging.debug("ClueGO query: " + str(len(each_primgene_list)) )
        
    drop_dupeprimgene_prot = {}
    unique_each_primgene_list = each_primgene_list
    genes_before_initial_drop = unique_each_primgene_list
    if len(each_primgene_list) != len(set(each_primgene_list)):
      dupe_gene_list = remove_list_duplicates(each_primgene_list)
      unique_each_primgene_list = [x for x in each_primgene_list if x.lower() not in [name.lower() for name in dupe_gene_list]]
      drop_dupeprimgene_prot = get_query_from_list(uniprot_query, dupe_gene_list)
      warning = []
      sorted_x = sorted(drop_dupeprimgene_prot.items(), key=lambda kv: kv[1])
      drop_dupeprimgene_prot = collections.OrderedDict(sorted_x)
      for each_dropped_prot in drop_dupeprimgene_prot:
        merged_out_dict[each_dropped_prot]['Comment'] += "Duplicate primary gene;"
        warning.append(each_dropped_prot + "(" + drop_dupeprimgene_prot[each_dropped_prot] + ")")
      if cy_debug:
        logging.debug("Uniprot duplicate primary gene: " + str(len(drop_dupeprimgene_prot)))
        logging.warning("WARNING - Dropping queries: " + ','.join(warning))
  
    if not unique_each_primgene_list:
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
      eprint("Error: No query for String and Genemania")
      sys.exit(1)
    
    string_filtered_dict = {}
    genemania_filtered_dict = {}
    
    if cy_run.lower() == "string" or cy_run.lower() == "both":
      if cy_debug:
        logging.debug("\nStep 3: String mapping started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
        logging.debug("String query: " + str(len(unique_each_primgene_list)))
 
      #Send gene list to string, get mapping & interactions
      string_interaction, string_unique_vertex, string_mapping, merged_out_dict = create_string_cytoscape(uniprot_query,unique_each_primgene_list, tax_id, cy_lim, cy_score, cy_debug, logging, merged_out_dict, genes_before_initial_drop)
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
      genemania_interaction, genemania_unique_vertex, genemania_mapping, merged_out_dict = create_genemania_interactions(uniprot_query,unique_each_primgene_list,tax_id,cy_lim,"10",cy_debug,logging, merged_out_dict, cy_session, cy_out, cy_cluego_out, genes_before_initial_drop)
      #Categorize interaction nodes as primary gene, secondary gene, external synonym and external interactor
      genemania_category = categorize_gene(genemania_unique_vertex, genemania_mapping, uniprot_query)
      #Get a filtered dictionary of interactions(primary-primary; primary-external interactor; external interactor-external interactor)
      genemania_filtered_dict, merged_out_dict = get_search_dicts(genemania_interaction, genemania_category, logging, cy_debug, uniprot_query, genemania_mapping, merged_out_dict, "Genemania")
      # Update merged_out_dict
      merged_out_dict = get_interactions_dict(genemania_filtered_dict, 'Genemania', merged_out_dict)
           
    #Merge String and Genemania interactions
    if cy_debug:
      logging.debug("\nStep 5: Merge String and Genemania interactions started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
  
    if not string_filtered_dict and not genemania_filtered_dict:
      eprint("Error: No interactions found in String and Genemania")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
      sys.exit(1)
   
    unique_merged_interactions = []
    unique_nodes = []
    #Get a list of merged interactions (ex: [node1 node2, node3 node4] etc) and a list of unique nodes (ex: [node1, node2, node3, node4])
    if cy_run.lower() == "string" or cy_run.lower() == "both":
      unique_merged_interactions,unique_nodes = get_merged_interactions(string_filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, cy_type_num)
    if cy_run.lower() == "genemania" or cy_run.lower() == "both":
      unique_merged_interactions,unique_nodes = get_merged_interactions(genemania_filtered_dict, unique_merged_interactions, unique_nodes, max_FC_len, each_category, uniprot_query, cy_type_num)
    
    if cy_debug:
      logging.debug("Total merged query nodes: " + str(len([i for i in unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])))
      logging.debug("Total merged interactions: " + str(len(unique_merged_interactions)))
    
    # Get uniprot query, primary gene and FC values together
    uniprot_list = {}
    if not cy_cluego_inp_file:
      for each_node in unique_nodes:
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num, site_info_dict, ambigious_sites, ambigious_genes)    
    else:
      lower_unique_each_primgene_list = [x.lower() for x in unique_each_primgene_list]
      for each_node in lower_unique_each_primgene_list:
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num, site_info_dict, ambigious_sites, ambigious_genes)
    
    # Interactors styling 
    fc_merged_vertex = {}
    if not (cy_type_num == "5" or cy_type_num == "6"):         
      cy_interactors_style(unique_nodes, unique_merged_interactions, uniprot_list, max_FC_len, each_category, cy_pval)
    else:     
      fc_merged_vertex = cy_sites_interactors_style(unique_nodes, unique_merged_interactions, uniprot_list, max_FC_len, each_category, cy_pval, cy_type_num)
    
    # Category styling   
    if cy_type_num == "4":    
      cy_category_style(unique_nodes, unique_merged_interactions, uniprot_list, each_category)
    
    if not cy_cluego_inp_file:
      if cy_debug:
        logging.debug("\nStep 6: ClueGO started at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
  
      #Start ClueGO  
      requests.post('http://localhost:1234/v1/apps/cluego/start-up-cluego')

      if cy_debug:
        # Number of ClueGO query + EI = x + y
        logging.debug("ClueGO query + EI: " + str(len([i for i in unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])) + " + " + str(len([i for i in unique_nodes if i.lower() not in [y.lower() for y in unique_each_primgene_list] ])))  

      filtered_unique_nodes, merged_out_dict = cluego_filtering(unique_nodes, cy_map, uniprot_query, cy_debug, logging, merged_out_dict, unique_each_primgene_list, cy_session, cy_out, cy_cluego_out)
    
      if cy_debug:
        # Number of ClueGO query + EI = x + y
        logging.debug("Total ClueGO query + EI: " + str(len([i for i in filtered_unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])) + " + " + str(len([i for i in filtered_unique_nodes if i.lower() not in [y.lower() for y in unique_each_primgene_list] ])))      
    
      final_length = len([i for i in filtered_unique_nodes if i.lower() in [x.lower() for x in unique_each_primgene_list] ])
      coverage = final_length/initial_length *100
      if cy_debug:
        logging.debug("\nQuery coverage = " + str(round(coverage, 2)) + "%")
      cluego_run(organism_name,cy_cluego_out,filtered_unique_nodes,cy_cluego_grouping,select_terms, leading_term_selection,cluego_reference_file,cluego_pval)
    
    if leading_term_cluster:
      cy_pathways_style(leading_term_cluster, each_category, max_FC_len, cy_pval, uniprot_list, cy_type_num, fc_merged_vertex)
    
    if cy_debug:
      logging.debug("\n\nRun completed sunccessfully")
      
    ## Write into outfile
    write_into_out(merged_out_dict, cy_out)
    requests.post("http://localhost:1234/v1/session?file=" + cy_session)
    requests.get("http://localhost:1234/v1/commands/command/quit")
    
  except Exception as e:
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    cytoscape_not_open_msg = "No connection could be made because the target machine actively refused it"
    cytoscape_not_open_msg2 = "Remote end closed connection without response"
    cytoscape_not_responding_msg = "Expecting value: line 1 column 1 (char 0)"
    if cytoscape_not_open_msg in str(e) or cytoscape_not_open_msg2 in str(e):
      eprint("Error: Cytoscape must be open")
      sys.exit(1)
    elif cytoscape_not_responding_msg in str(e):
      traceback.print_exc()
      eprint("Error: Cytoscape not responding. Please restart and wait for it to fully load")
      sys.exit(1)
    else:
      traceback.print_exc()
      sys.exit(1)
      
if __name__ == "__main__":
  main(sys.argv[1:])