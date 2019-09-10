#usage: python get_interactions_table_4.py -i Mouse_rev_unrev_input.csv -t nofc -s rat -m "C:\Users\SundararamN\ClueGOConfiguration\v2.5.5\ClueGOSourceFiles\Organism_Rattus norvegicus\Rattus norvegicus.gene2accession_2019.05.20.txt.gz" -d
#Goal: Given list of protein IDs (with/without their corresponding FC and pval), construct 1) an interaction network among all proteins in list 2) construct a pathway network of proteins in the list 
try:
  from py2cytoscape import cyrest
except ImportError:
  print("Error: Please install module py2cytoscape. [Installation: pip install py2cytoscape]")
from py2cytoscape.data.cynetwork import CyNetwork
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.style import StyleUtil
try:
  import requests
except ImportError:
  print("ImportError: Please make sure you are using Python3")
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
import tkinter as tk
from tkinter import filedialog as fd

def setup_logger(name, log_file, level=logging.DEBUG):
  handler = logging.FileHandler(log_file,mode='w')
  logger = logging.getLogger(name)
  logger.setLevel(level)
  logger.addHandler(handler)
  return logger

def preprocessing(inp, type, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out):
  is_prot_col = False
  is_FC = False
  is_pval = False
  is_cat = False
  is_multFC = False
  is_label_col = False
  each_protein_list = []
  prot_list = {}
  max_FC_len = 0
  unique_labels = []
  each_category = []
  initial_length = 0
  repeat_prot_ids = []
  retain_prot_ids = []
  repeat_prot_ids_2 = {}
  to_return_unique_protids_length = 0
  
  with open(inp,'r') as csv_file:
    '''
    Read input and collect columns based on type of analysis chosen
    Types:1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = category
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
          if type == "2":
            is_multFC = True
            if "label" in row[i].lower():
              label = i  
              is_label_col = True    
        
        if not is_prot_col:
          print("Error: Column 'ProteinID' not found")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit()          
      else: 
        if '/' in row[protein]:
          split_val = row[protein].split('/')
          if len(split_val) > 2:
            print("Error: Multiple protein IDs not considered")
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
            sys.exit()
          elif len(split_val) == 2:
            if not bool(re.match('^[0-9]$', split_val[0])):
              print("Error: Multiple protein IDs not considered")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
              sys.exit()
        
        if '|' in row[protein]:
          row[protein] = row[protein].split('|')[1]
          
        if not bool(re.match('^[A-Za-z0-9\-]+$', row[protein])):
          err_msg = "Error: Invalid proteinID: " + row[protein]
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit()
        
        if type == "1" or type == "2":
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
              
        if row[protein] not in each_protein_list:
          each_protein_list.append(row[protein])
          if type == "1" or type == "2" or type == "3":          
            if type == "1":
              if is_prot_col and is_FC:           
                prot_list.update({row[protein]:[[float(get_fc_val)],[get_pval]]})
                max_FC_len = 1
              else:
                print("Error: Required columns- ProteinID, FC")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit()
                
            elif type == "2":
              if is_prot_col and is_FC and is_label_col:
                prot_list[row[protein]] = {}
                prot_list[row[protein]].update({row[label]:[float(get_fc_val),get_pval]})
                if row[label] not in unique_labels:
                  unique_labels.append(row[label])                
              else:
                print("Error: Required columns- ProteinID, FC, Label")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit()
                
            elif type == "3":
              if is_prot_col:
                continue
              else:
                print("Error: Required columns- ProteinID")
                remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
                sys.exit()
            
          elif type == "4":
            if (not is_prot_col) or (not is_cat):
              print("Error: Required columns- ProteinID, Category")
              remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
              sys.exit()
              
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
  
  if len(unique_labels) > 6:
    print("Error: Number of unique labels should not exceed 5")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit()
    
  if len(each_category) > 6:
    print("Error: Number of categories should not exceed 5")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit()
  
  with open(cy_out,'w') as csv_file:
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
      if retain_prot_ids:
        each_protein_list = list(set(unique_each_protein_list))
        print(each_protein_list)
        
      all_dropped = dropping_repeats + retain_prot_ids
      all_dropped = sorted(all_dropped)
      if cy_debug:
        logging.debug("Duplicate query: " + str((initial_length)-len(each_protein_list)))
        logging.warning("WARNING - Dropping queries: " + ','.join(all_dropped))
      for each_dupe_query in all_dropped:
        line = each_dupe_query + ",,,,," + "Duplicate query;\n"
        csv_file.write(line) 
        
    elif type == "3":
      if repeat_prot_ids:
        unique_each_protein_list = list(set(each_protein_list))    
        if cy_debug:
          logging.debug("Duplicate query: " + str(len(each_protein_list)-len(unique_each_protein_list)))
          logging.warning("WARNING - Dropping queries: " + ','.join(repeat_prot_ids)) 
        for each_dupe_query in repeat_prot_ids:
          line = each_dupe_query + ",,,," + "Duplicate query;\n"
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
          csv_file.write(line) 
      each_protein_list = unique_each_protein_list
        
    elif type == "2":
        unique_each_protein_list = list(set(each_protein_list))
        count_dropped = 0
        additional_dropped = []
        if retain_prot_ids:
          for each_dupe_query in retain_prot_ids:
            line = each_dupe_query + ",,,," + "Duplicate query;\n"
            csv_file.write(line)
        if repeat_prot_ids_2:         
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
        
        list_of_duplicates = additional_dropped + list(repeat_prot_ids_2.keys()) + retain_prot_ids
        list_of_duplicates = sorted(list_of_duplicates)
        if cy_debug:
          logging.debug("Duplicate query: " + str(len(list_of_duplicates)))
          logging.warning("WARNING - Dropping queries: " + ','.join(list_of_duplicates)) 
        
        for each_dupe_query in list_of_duplicates:
          line = each_dupe_query + ",,,," + "Duplicate query;\n"
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
    
    if cy_debug:
      logging.debug("Unique remaining queries: " + str(len(each_protein_list)))      
  csv_file.close()
 
  # Limit query inpt number = 1100
  if len(each_protein_list) > 1100:
    print("Error: The query input is too big. Currently supporting upto 1000 query protein ids")
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit()  
  
  return(each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict, to_return_unique_protids_length)

def inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict):
  queries_dropped = []
  for each_prot in list(prot_list):
    delete_each_prot = False
    for each_fc_val,each_pval in zip(prot_list[each_prot][0],prot_list[each_prot][1]):
      if not (abs(each_fc_val) >= abs(cy_fc_cutoff) and each_pval <= cy_pval_cutoff):
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

  params = {
  'from':'ACC+ID',
  'to':'ACC',
  'format':'tab',
  'query':query_term,
  'columns': 'id,genes(PREFERRED),genes(ALTERNATIVE),organism'
  }
  data = urllib.parse.urlencode(params).encode("utf-8")
  request = urllib2.Request(url, data)
  response = urllib2.urlopen(request)
  page = response.read(200000)
  decode = page.decode("utf-8")
  list1=decode.split('\n')
  list1 = list1[1:]
  for each_list1 in list1:
    if each_list1 != '':
      uniprot_list = each_list1.split('\t')
      uniprot_protid = uniprot_list[0]
      prot = uniprot_list[4]
      split_prot_list = prot.split(',')
      for each_prot in split_prot_list:
        primary_gene = uniprot_list[1]
        merged_out_dict[each_prot] = {}
        comment_merged = ""
        # Pick first gene in case of multiple primary genes for single uniprot ids 
        if ";" in primary_gene:
          prot_with_mult_primgene.append(each_prot + "(" + primary_gene + ") ")
          #comment_merged = "Multiple primary genes;"
          primary_gene = primary_gene.split(";")[0]
        synonym_gene = uniprot_list[2]
        synonym_gene = synonym_gene.split(" ")
        organism_name = uniprot_list[3]
        uniprot_query[each_prot] = {}
        uniprot_query[each_prot].update({'Uniprot':uniprot_protid,'Primary':primary_gene,'Synonym':synonym_gene,'Organism':organism_name})
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
    print("Error: Protein list is of more than 1 organism: " + unique_organisms)
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit()
  
  if not species.lower() in unique_organisms[0].lower():
    print("Error: Species mismatch. Species parameter provided is " + species + " and species of protein list is " + unique_organisms[0])
    remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
    sys.exit()
   
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
  
  return(uniprot_query,each_primgene_list,merged_out_dict)

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
      print("Error: Genemania timed out- please try again later")
      sys.exit()
    except:
      print("Error: Cytoscape must be open")
      print(e)
      sys.exit()
      
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
  
def get_everything_together(each,uniprot_query, uniprot_list, max_FC_len, each_category, type):
  '''
  For genes in merged interaction list, append other info: FC, pval, category etc
  '''  
  not_in_list = 1
  if not uniprot_list:
    uniprot_list.update({'name':[each]})
    length_each = len(each)
    uniprot_list.update({'length':[length_each]})
    not_in_list = 0
  else:
    uniprot_list['name'].append(each)
    length_each = len(each)
    uniprot_list['length'].append(length_each)
  
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
          print("Error: Required columns in ClueGO file: SymbolID and UniqueID#EntrezGeneID")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit()
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
      min_number_of_genes_per_term = 5     
    else:
      min_number_of_genes_per_term = 20
    min_percentage_of_genes_mapped = 0
    min_go = 1
    max_go = 6
    kappa = 0.5
  elif group.lower() == "medium":
    if len(merged_vertex) <= 500:
      min_number_of_genes_per_term = 5
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
  response = requests.get(CLUEGO_BASE_URL+SEP+"analysis-results"+SEP+"get-cluego-table"+SEP+str(current_network_suid))
  table_file_name = output_cluego
  writeLines(response.text,table_file_name)

def cluego_pick_clusters(cy_cluego_out):
  cluster = {}
  leading_term_cluster = {}
  
  total_terms = 0
  with open(cy_cluego_out,'r') as input_file:
    in_file = csv.reader(input_file, delimiter='\t')
    line_count = 0
    for row in in_file:
      if line_count == 0:
        header = row
        # Get header
        for i in range(len(row)):
          if row[i].lower() == "leadinggoterm":
            leading = i
          elif row[i].lower() == "goterm":
            goterm = i
          elif row[i].lower() == "associated genes found":
            genes = i
          elif row[i].lower() == "gogroups":
            gogroup = i
      else:
        each_gene_list = ((row[genes].replace("[","")).replace("]","")).replace(" ","")
        each_gene_list = each_gene_list.split(',')    
        list_gogroup = ((row[gogroup].replace("[","")).replace("]","")).replace(" ","")
        list_gogroup = list_gogroup.split(',')
        if row[gogroup] not in cluster:
          cluster[row[gogroup]] = {}
          cluster[row[gogroup]].update({"NoTerms":1})         
        else:
            cluster[row[gogroup]]["NoTerms"] += 1
        if row[leading].lower() == "true":        
          cluster[row[gogroup]].update({"GOTerm":row[goterm],"Genes":each_gene_list})        
        total_terms +=1        
      line_count += 1  
  
  total_no_genes = sorted(list(set([cluster[x]['NoTerms'] for x in cluster ])), reverse = True)
  top_5 = total_no_genes[:5]
  
  for each_cluster in cluster:
    if cluster[each_cluster]['NoTerms'] in top_5:
      already_present = False
      for name in leading_term_cluster:
        if cluster[each_cluster]["GOTerm"].lower() == name.lower():
          already_present = True
          if not (len(leading_term_cluster[name]) > len(cluster[each_cluster]["Genes"])):
            leading_term_cluster[name] = cluster[each_cluster]["Genes"]
            
      if not already_present:  
        leading_term_cluster.update({cluster[each_cluster]["GOTerm"]:cluster[each_cluster]["Genes"]})
  return(leading_term_cluster)
  
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
          print("Error: Required columns in cluego input file: GOTerm and Associated Genes Found")
          remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
          sys.exit()            
      
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
            print("Error: Duplicate GOID found in cluego input file")
            remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
            sys.exit()
          
      line_count+=1
  return(top_annotations, unique_gene)

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
    my_style = multipleFC(my_style,uniprot_list['FC_exists'],uniprot_list["query"],"1",uniprot_list['name'],max_FC_len)
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
        "lesser": "#FF6633",
        "equal": "#FF0000",
        "greater": "#FF0000"
      },
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
  
def multipleFC(my_style,FC_exists,query,func,name,max_FC_len):
  '''
  Styling specific to multipleFC case
  '''
  color_code = ["#3366FF", "#33FFFF", "#FF6600", "#FFFF66", "#FF0000", "#006666", "#33FF33", "#FFCCCC", "#3300FF", "#CCCCFF"] 
  bar_columns = ""
  color_columns = ""
  for i in range(1,max_FC_len+1):
    term_FC = 'FC' + str(i)
    bar_columns = bar_columns + '\"' + term_FC + '\"' + ","
    color_columns = color_columns + '\"' + color_code[i-1] + '\"' + ","
  bar_columns = bar_columns[:-1]
  color_columns = color_columns[:-1]
  
  value = "org.cytoscape.BarChart:{" + '\"' + "cy_colors" + '\"' + ":[" + color_columns + "]," + '\"' + "cy_dataColumns" + '\"' + ":[" + bar_columns + "]}"
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

def cy_pathways_style(cluster, each_category, max_FC_len, pval_style, uniprot_list):
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
  
  pos_or_neg_func_value = {}
  merged_vertex = []
  query_val_noFC = []
  count_each = 0
    
  category_present = 0
  for each in each_category:
    category_present = 1
    break
      
  for each in cluster_list:
    if each not in merged_vertex:
      G.add_vertex(each)
      merged_vertex.append(each)
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
    
    val_pos_or_neg = 0
    for each_gene in cluster_list[each]:
      if each_gene not in merged_vertex:
        G.add_vertex(each_gene)
        merged_vertex.append(each_gene)
        if each_gene in function_only:
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
          if max_FC_len == 0 and not category_present:
            index_noFC = uniprot_list["name"].index(each_gene.lower())
            query_val_noFC.append((uniprot_list["query"])[index_noFC])
          val_length_of_val = len(each_gene) #* 15
          length_of.append(val_length_of_val)
          val_breadth_of_val = 30
          breadth_of.append(val_breadth_of_val)
          if max_FC_len == 1:
            if each_gene in uniprot_list["name"]:
              index = uniprot_list["name"].index(each_gene.lower())
              FC_val_each_gene = (uniprot_list['FC1'])[index]
              if FC_val_each_gene > 0:
                val_pos_or_neg += 1
              elif FC_val_each_gene < 0:
                val_pos_or_neg -= 1
            
      name_edge = each + " with " + each_gene
      G.add_edge(each,each_gene,name=name_edge)

    pos_or_neg_func_value.update({each:val_pos_or_neg})
   
  G.vs
  G.vs["query"] = query
  G.vs["name"] = merged_vertex
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
  for i in range(1,max_FC_len+1):
    k = 0
    term_FC = 'FC' + str(i)
    term_pval = 'pval' + str(i)
    add_term_FC = []
    add_term_pval = []
    significant_val = []
    for each_vertex_name in merged_vertex:
      if each_vertex_name.lower() in uniprot_list["name"]:
        index = uniprot_list["name"].index(each_vertex_name.lower())
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
        if max_FC_len == 1:
          add_term_value = pos_or_neg_func_value[each_vertex_name]
          if add_term_value >= 0:
            add_term_value = 100
          elif add_term_value < 0:
            add_term_value = -100
          add_term_FC.append(add_term_value)
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
    #'NODE_FILL_COLOR':"#33CCFF",
    'EDGE_TRANSPARENCY':"15",
    'NODE_BORDER_PAINT':"#999999",
    'EDGE_WIDTH':"2"
  }
  my_style.update_defaults(basic_settings)
  my_style.create_passthrough_mapping(column='shared name', vp='NODE_LABEL', col_type='String')
  
  shape_kv_pair = {
    "Function":"ELLIPSE",
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

  pval_sig = 0
  if max_FC_len == 1:
    my_style = singleFC(my_style, uniprot_list)
    if pval_style:
      pval_sig = 1
  
  if max_FC_len > 1:
    my_style = multipleFC(my_style,FC_exists,query,"2",merged_vertex, max_FC_len)

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
    
def main_(argv):
  Gui().start()

class Gui:
  def __init__(self):
    self.root = tk.Tk()
    
    self.cy_in = tk.StringVar()
    self.cy_out = tk.StringVar()
    self.cy_map = tk.StringVar()
    self.select_terms = tk.StringVar()
    self.cy_type = tk.StringVar()
    self.cy_species = tk.StringVar()
    self.cy_score = tk.StringVar()
    self.cy_cluego_grouping = tk.StringVar()
    self.cy_debug = tk.BooleanVar()
    self.cluego_pval = tk.StringVar()
    
    # no interface yet
    self.cy_lim = 0
    self.logging = ""
    self.merged_out_dict = {}
    self.help = False
    self.cy_pval = 0
    self.cy_fc_cutoff = 0
    self.cy_pval_cutoff = 1
    self.leading_term_selection = "no. of genes per term"
    self.cluego_reference_file = ""
    self.cy_cluego_inp_file = ""
    self.cy_session = "PINE_" + str(datetime.now().strftime("%d%b%Y_%H%M%S")) + ".cys"
    self.cy_run = "both"
    
  def start(self):   
    frame = tk.Frame(self.root, relief=tk.SUNKEN)
    frame.grid(row=0, column=0, padx=5, pady=5)
    tk.Button(frame, text="Choose an input file", command=self.get_input_filename).pack()
    tk.Label(frame, text="Input file:").pack()
    tk.Label(frame, textvariable=self.cy_in).pack()
    
    frame = tk.Frame(self.root, relief=tk.SUNKEN)
    frame.grid(row=0,column=1, padx=5, pady=5)
    tk.Button(frame, text="Choose an output directory", command=self.get_output_directory).pack()
    tk.Label(frame, text="Output directory:").pack()
    tk.Label(frame, textvariable=self.cy_out).pack()
    
    frame = tk.Frame(self.root, bd=1)
    frame.grid(row=0, column=2, padx=5, pady=5)
    tk.Button(frame, text="Choose the ClueGO mapping file", command=self.get_cluego_mapping_filename).pack()
    tk.Label(frame, text="ClueGO mapping file:").pack()
    tk.Label(frame, textvariable=self.cy_map).pack()
    
    frame = tk.Frame(self.root, bd=1)
    frame.grid(row=1, column=0)
    allowed_select_terms = ["biological process", "subcellular location", "molecular function", "pathways", "all"]
    self.select_terms.set(allowed_select_terms[3])
    tk.Label(frame, text="Select term:").pack()
    tk.OptionMenu(frame, self.select_terms, *allowed_select_terms).pack()
    
    allowed_types = ["noFC", "singleFC", "multiFC", "category"]
    self.cy_type.set(allowed_types[2])
    tk.Label(frame, text="Select type:").pack()
    tk.OptionMenu(frame, self.cy_type, *allowed_types).pack()
    
    allowed_species = ["human", "mouse", "rat"]
    self.cy_species.set(allowed_types[1])
    tk.Label(frame, text="Select species:").pack()
    tk.OptionMenu(frame, self.cy_species, *allowed_species).pack()
    
    frame = tk.Frame(self.root, bd=1)
    frame.grid(row=1, column=1)
    tk.Label(frame, text="Enter an interaction confidence score").pack()
    tk.Entry(frame, textvariable=self.cy_score).pack()
    
    allowed_cluego_grouping = ["global", "medium", "detailed"]
    self.cy_cluego_grouping.set(allowed_cluego_grouping[1])
    tk.Label(frame, text="Select ClueGO grouping:").pack()
    tk.OptionMenu(frame, self.cy_cluego_grouping, *allowed_cluego_grouping).pack()
    
    tk.Checkbutton(frame, text="Debug mode", variable=self.cy_debug).pack()
    
    tk.Label(frame, text="Enter a ClueGO p-value").pack()
    tk.Entry(frame, textvariable=self.cluego_pval).pack()
    
    frame = tk.Frame(self.root)
    frame.grid(row=2, column=0)
    tk.Button(frame, text="Save session", command=self.save_session).pack()
    tk.Button(frame, text="Load session", command=self.load_session).pack()
    
    frame = tk.Frame(self.root, bd=1)
    frame.grid(row=2, column=1)
    tk.Button(frame, text="Run", command=self.run).pack()
    
    self.root.mainloop()
    
  def get_input_filename(self):
    filename = fd.askopenfilename()
    self.cy_in.set(filename)
    
  def get_output_directory(self):
    filename = fd.askdirectory()
    self.cy_out.set(filename)
    
  def get_cluego_mapping_filename(self):
    filename = fd.askopenfilename()
    self.cy_map.set(filename)
    
  def save_session(self):
    filename = fd.asksaveasfilename(filetypes=[("Session files", "*.session")])
    if not filename:
      return
    if not filename.endswith(".session"):
      filename += ".session"
    
    session_data = {
      "cy_in": self.cy_in.get(),
      "cy_out": self.cy_out.get(),
      "cy_map": self.cy_map.get(),
      "select_terms": self.select_terms.get(),
      "cy_type": self.cy_type.get(),
      "cy_species": self.cy_species.get(),
      "cy_score": self.cy_score.get(),
      "cy_cluego_grouping": self.cy_cluego_grouping.get(),
      "cy_debug": self.cy_debug.get(),
      "cluego_pval": self.cluego_pval.get(),
    }    
    
    with open(filename, "w") as f:
      f.write(json.dumps(session_data))
      
  def load_session(self):
    filename = fd.askopenfilename(filetypes=[("Session files", "*.session")])
    if not filename:
      return
    
    with open(filename) as f:
      session_data = json.loads(f.read())
      self.cy_in.set(session_data["cy_in"])
      self.cy_out.set(session_data["cy_out"])
      self.cy_map.set(session_data["cy_map"])
      self.select_terms.set(session_data["select_terms"])
      self.cy_type.set(session_data["cy_type"])
      self.cy_species.set(session_data["cy_species"])
      self.cy_score.set(session_data["cy_score"])
      self.cy_cluego_grouping.set(session_data["cy_cluego_grouping"])
      self.cy_debug.set(session_data["cy_debug"])
      self.cluego_pval.set(session_data["cluego_pval"])
    
  def run(self):
    if self.cy_debug.get():
      self.logging = setup_logger("cytoscape","Cytoscape.log")
    run(
      self.cy_in.get(),
      self.cy_species.get(),
      self.cy_lim,
      self.cy_score.get(),
      self.cy_debug.get(),
      self.logging,
      self.cy_map.get(),
      self.cy_out.get() + "/Merged_Input.csv", # cy_out
      self.cy_out.get() + "/ClueGo_Output.txt", # cy_cluego_out
      self.merged_out_dict,
      self.help,
      self.cy_type.get(),
      self.cy_pval,
      self.cy_cluego_grouping.get(),
      self.cy_fc_cutoff,
      self.cy_pval_cutoff,
      self.select_terms.get(),
      self.leading_term_selection,
      self.cluego_reference_file,
      self.cy_cluego_inp_file,
      self.cy_out.get() + "/PINE_" + str(datetime.now().strftime("%d%b%Y_%H%M%S")) + ".cys", # cy_session
      self.cluego_pval.get(),
      self.cy_run
    )
    
def main(argv):
  cy_in = ""
  cy_species = ""
  cy_lim = 0
  cy_score = 0.4
  cy_debug = False
  logging = ""
  cy_map = ""
  cy_out = "Merged_Input.csv"
  cy_cluego_out = "ClueGO_Output.txt"
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
  cy_session = "PINE_" + str(datetime.now().strftime("%d%b%Y_%H%M%S")) + ".cys"
  cluego_pval = 0.05
  cy_run = "both"
  
  try:
    opts, args = getopt.getopt(argv, "i:s:l:t:r:dm:o:c:ng:f:p:z:e:h:a:v:y:u:",["--in=","--species=","--limit=","--type=","--score=","--debug","--mapping=","--output=","--output_cluego=","--significant","--grouping=","--fccutoff=","--pvalcutoff=","--visualize=","leading-term=","--reference-path=","--input_cluego=","--save-session=","--cluego-pval=","--run="])
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
        cy_out = arg
      elif opt in ("-c","--output_cluego"):
        cy_cluego_out = arg
      elif opt in ("-a","--input_cluego"):
        cy_cluego_inp_file = arg
      elif opt in ("-e","--leading-term"):
        leading_term_selection = arg
      elif opt in ("-z","--visualize"):
        select_terms = arg
      elif opt in ("-y","--cluego-pval"):
        cluego_pval = arg
      elif opt in ("-h","--reference-path"):
        cluego_reference_file = arg
      elif opt in ("-g","--grouping"):
        cy_cluego_grouping = arg
      elif opt in ("-v","--save-session"):
        cy_session = arg
      elif opt in ("-d","--debug"):
        cy_debug = True
        logging = setup_logger("cytoscape","Cytoscape.log")
      else:
        help = True      
  except getopt.GetoptError:
    help = True
    
  if not cy_in or not cy_species or not cy_type:
    help = True
    
  if help:
    print("Cytoscape Visualization")
    print("---------------------------------------------------------------------------------------------")
    print("Usage:         cytoscape_api.py -i input.csv -o merged.csv -c cluego_out.txt -t input_type -s species -m cluego_map_file.gz")
    print("Argument:      -i [--in]: input file in csv format")
    print("Argument:      -o [--output]: output merged interactions information file")     
    print("Argument:      -t [--type]: input type [Allowed: noFC, singleFC, multiFC, category]") 
    print("Argument:      -s [--species]: species [Supported: human, mouse, rat]")  
    print("Argument:      -m [--mapping]: path to cluego mapping file compressed in .gz format")
    print("Argument:      -v [--save-session]: path to save cytoscape session")
    print("Argument(opt): -l [--limit]: number of interactors [Default:0, Range:0-100]")
    print("Argument(opt): -r [--score]: interaction confidence score [Default:0.4, Range 0-1]")
    print("Argument(opt): -c [--output_cluego]: output cluego file in txt format")
    print("Argument(opt): -a [--input_cluego]: cluego input file in txt format")
    print("Argument(opt): -g [--grouping]: cluego grouping [Allowed: global, medium, detailed; Default: medium]")
    print("Argument(opt): -z [--visualize]: cluego analysis type [Allowed: biological process, subcellular location, molecular function, pathways, all; Default: pathways]")
    print("Argument(opt): -e [--leading-term]: cluego leading term selection criteria [Allowed: highest significance, no. of genes per term, percent of genes per term, percent genes per term vs cluster; Default: no. of genes per term]")
    print("Argument(opt): -h [--reference-path]: cluego custom reference path")
    print("Argument(opt): -f [--fccutoff]: cutoff for fold change [Default: abs(FC) >= 0.0]")
    print("Argument(opt): -p [--pvalcutoff]: cutoff for p value [Default: pval > 1.0]")
    print("Argument(opt): -n [--significant]: outline statistically significant nodes")  
    print("Argument(opt): -d [--debug]: generate logging file called cytoscape.log")
    print("Argument(opt): -h [--help]: display help")
    sys.exit()
    
  run(cy_in, cy_species, cy_lim, cy_score, cy_debug, logging, cy_map, cy_out, cy_cluego_out, merged_out_dict, help, cy_type, cy_pval, cy_cluego_grouping, cy_fc_cutoff, cy_pval_cutoff, select_terms, leading_term_selection, cluego_reference_file, cy_cluego_inp_file, cy_session, cluego_pval, cy_run)

def run(cy_in, cy_species, cy_lim, cy_score, cy_debug, logging, cy_map, cy_out, cy_cluego_out, merged_out_dict, help, cy_type, cy_pval, cy_cluego_grouping, cy_fc_cutoff, cy_pval_cutoff, select_terms, leading_term_selection, cluego_reference_file, cy_cluego_inp_file, cy_session, cluego_pval, cy_run):
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
    print("Error: Species currently supported are human, mouse, rat")
    sys.exit()

  if not ('\\' in cy_session or "/" in cy_session):
    cwd = os.getcwd()
    cy_session = cwd + "\\" + cy_session
    
  if not ('.cys' in cy_session):
    print("Error: Cytoscape session file must have a valid name with .cys extension")
    sys.exit()
  
  if not (float(cluego_pval) >=0.0 and float(cluego_pval) <=1.0):
    print("Error: Cluego pvalue range must be between 0 to 1")
    sys.exit()
  
  allowed_runs = ["string","genemania","both"]
  if cy_run.lower() not in allowed_runs:  
    print("Error: String run type must be one of the following: " + ','.join(allowed_runs))
    sys.exit()
    
  # Highest%20Significance, %23Genes%20%2F%20Term, %25Genes%20%2F%20Term, %25Genes%20%2F%20Term%20vs%20Cluster
  allowed_leading_term = ["highest significance", "no. of genes per term", "percent of genes per term", "percent genes per term vs cluster"]
  if leading_term_selection.lower() not in allowed_leading_term:
    print("Error: ClueGO leading term selection must be one of the following: " + (',').join(allowed_leading_term))
    sys.exit()
  elif leading_term_selection.lower() == "highest significance":
    leading_term_selection = "Highest%20Significance"
  elif leading_term_selection.lower() == "no. of genes per term":
    leading_term_selection = "%23Genes%20%2F%20Term"
  elif leading_term_selection.lower() == "percent of genes per term":
    leading_term_selection = "%25Genes%20%2F%20Term"
  elif leading_term_selection.lower() == "percent genes per term vs cluster":
    leading_term_selection = "%25Genes%20%2F%20Term%20vs%20Cluster"
    
  allowed_type = ["singlefc","multifc","nofc","category"]
  # 1 = SingleFC; 2 = MultiFC; 3 = List only; 4 = category
  if cy_type.lower() not in allowed_type:
    print("Error: Input type must be one of the following: " + (',').join(allowed_type))
    sys.exit()
  elif cy_type.lower() == "singlefc":
    cy_type_num = "1"
  elif cy_type.lower() == "multifc":
    cy_type_num = "2"
  elif cy_type.lower() == "nofc":
    cy_type_num = "3"
  else:
    cy_type_num = "4"
    
  allowed_selections = ["biological process","subcellular location","molecular function","pathways","all"]
  if select_terms.lower() not in allowed_selections:
    print("Error: The visualization type must be one of the following: " + (',').join(allowed_selections))
    sys.exit()
  
  try:
    cy_fc_cutoff = float(cy_fc_cutoff)
  except:
    print("Error: FC cutoff must be a number") 
    sys.exit()

  try:
    cy_pval_cutoff = float(cy_pval_cutoff)
  except:
    print("Error: PVal cutoff must be a number") 
    sys.exit()  
  
  if not (cy_pval_cutoff >= 0.0 and cy_pval_cutoff <= 1.0):
    print("Error: PVal cutoff must range between 0 and 1")
    sys.exit()
    
  allowed_groups = ["global","medium","detailed"]
  if cy_cluego_grouping.lower() not in allowed_groups:
    print("Error: Cluego grouping must be one of the following:" + (',').join(allowed_groups))
    sys.exit()
    
  cy_score = float(cy_score)*1000
  if not (cy_score >= 0.0 and cy_score <= 1000.0):
    print("Error: Invalid string score provided. Value must be between 0 to 1; Confidence levels for string score- Low = 0.150, Medium = 0.400, High = 0.700, Highest = 0.900")
    sys.exit()
  
  # Rounding off score - ex: 500.9 and above = 501, less = 500
  cy_score_ceil = math.ceil(cy_score)
  if((cy_score_ceil-cy_score) <= 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = cy_score_ceil
  elif ((cy_score_ceil-cy_score) > 0.1 and (cy_score_ceil-cy_score) > 0.0):
    cy_score = math.floor(cy_score)
  
  if not (int(cy_lim) >= 0 and int(cy_lim) <=100):
    print("Error: Limit on additional interactors is 100. Please choose a number between 0 and 100")
    sys.exit()
    
  try:
    # Check Cytoscape version
    request = requests.get('http://localhost:1234/v1/version')
    cy_version = request.json()
    if not bool(re.match('^3.7', cy_version['cytoscapeVersion'])):
      print("Error: Cytoscape version must be 3.7.0 and above")
      sys.exit()
    
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
      print("Error: Cytoscape app GeneMANIA v3.5.1 not installed or not responding properly")
      sys.exit()
      
    if not app_cluego or not (ver_cluego == "2.5.4" or ver_cluego == "2.5.5"):
      print("Error: Cytoscape app ClueGO v2.5.4/v2.5.5 not installed or not responding properly")
      sys.exit()
    
    if cluego_reference_file and ver_cluego != "2.5.5":
      print("Error: Using ClueGO custom reference file needs version 2.5.5 of ClueGO")
      sys.exit()
    
    body = dict(offline=True)
    response = requests.post("http://localhost:1234/v1/commands/genemania/organisms", json=body)
    genemania_bool = False
    for each in response.json()['data']['organisms']:
      if tax_id == str(each['taxonomyId']):
        genemania_bool = True
    if not genemania_bool:
      print("Error: Please install " + cy_species + " dataset in Genemania")
      sys.exit()
    
    if not cy_cluego_inp_file:
      # Check if cluego path exists
      if not cy_map:
        print("Error: ClueGO mapping file path must be provided")
        sys.exit()
        
      if not path.exists(cy_map):
        print("Error: Path to ClueGO mapping file " + cy_map + " does not exist")
        sys.exit()
      else:
        if ver_cluego not in cy_map:
          print("Error: ClueGO version mismatch. Version installed is " + ver_cluego + " and version contained in path to ClueGO mapping file " + cy_map)
          sys.exit()
        if organism_name.lower() not in cy_map.lower():
          print("Error: Species mismatch.  Species parameter provided is " + organism_name + " and species contained in path to ClueGO mapping file is " + cy_map)
          sys.exit()
    
      if not (".gz" in cy_map or "gene2accession" in cy_map):
        print("Error: ClueGO mapping file must refer to the species gene2accession gz file")
        sys.exit()
    
    #Read input and obtain protid 
    if cy_debug:
      logging.debug("Step 1: Start processing the input protein list at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
    unique_each_protein_list, prot_list, max_FC_len, each_category, merged_out_dict,initial_length = preprocessing(cy_in, cy_type_num, cy_debug, logging, merged_out_dict, cy_out, cy_session, cy_cluego_out)
    
    # FC and Pval cutoff
    if not (cy_fc_cutoff == 0.0 and cy_pval_cutoff == 1.0):
      unique_each_protein_list, prot_list, merged_out_dict = inp_cutoff(cy_fc_cutoff, cy_pval_cutoff, unique_each_protein_list, prot_list, cy_debug, logging, merged_out_dict)
    
    #Uniprot API call to get primary gene, synonym
    if cy_debug:
      logging.debug("\nStep 2: Start the uniprot api call at " + str(datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")))
      logging.debug("Uniprot query: " + str(len(unique_each_protein_list)))
        
    uniprot_query,each_primgene_list,merged_out_dict = uniprot_api_call(unique_each_protein_list, prot_list, cy_type_num, cy_debug, logging, merged_out_dict, organism_name, cy_session, cy_out, cy_cluego_out)
        
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
      print("Error: No query for String and Genemania")
      sys.exit()
    
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
      print("Error: No interactions found in String and Genemania")
      remove_out(cy_debug, logging, cy_session, cy_out, cy_cluego_out)
      sys.exit()
   
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
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num)    
    else:
      lower_unique_each_primgene_list = [x.lower() for x in unique_each_primgene_list]
      for each_node in lower_unique_each_primgene_list:
        uniprot_list = get_everything_together(each_node, uniprot_query, uniprot_list, max_FC_len, each_category, cy_type_num)
    
    # Interactors styling    
    cy_interactors_style(unique_nodes, unique_merged_interactions, uniprot_list, max_FC_len, each_category, cy_pval)
    
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
      leading_term_cluster = cluego_pick_clusters(cy_cluego_out)
    
    if leading_term_cluster:
      cy_pathways_style(leading_term_cluster, each_category, max_FC_len, cy_pval, uniprot_list)
    
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
      print("Error: Cytoscape must be open")
      traceback.print_exc()
    elif cytoscape_not_responding_msg in str(e):
      print("Error: Cytoscape not responding. Please restart and wait for it to fully load")
    else:
      traceback.print_exc()
      
if __name__ == "__main__":
  main(sys.argv[1:])