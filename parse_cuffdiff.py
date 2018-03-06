#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################
# Modules
#########################

import __main__ as main
from os import path
import pandas as pd
import sys

#########################
# Global Variables
#########################

input_genes_fpkm='genes.fpkm_tracking'
input_genes_fpkm_delim='\t'
output_genes_fpkm='hidata_genes.tsv'

input_genes_rg='genes.read_group_tracking'
input_genes_rg_delim='\t'
output_genes_rg_exp='expressed_genes.tsv'

input_gene_exp='gene_exp.diff'
input_gene_exp_delim='\t'
output_gene_exp_sig='significant_genes.tsv'

script_purpose = 'Parse Cuffdiff output into several files.'
script_name = path.basename(main.__file__)
script_path = path.dirname(path.realpath(__file__))

#########################
# Functions
#########################

def print_message(message, exit=False):
	"""Print a message to stdout. Optionally exit."""
	if exit:
		sys.exit(message)
	print(message)

def insert_row(df, data):
	"""Insert data as new row at top of df.

	Keyword arguments:
	df -- a Pandas dataframe.
	data -- a list of data to insert.
	"""
	# If data length is less then number columns, insert blanks.
	num_blanks = len(df.columns.values)-len(data)
	for i in range(num_blanks):
		data.extend([''])
	df.loc[-1] = data
	# Shift index
	df.index = df.index + 1
	# Sort by index
	df = df.sort_index()
	return df

def make_gct(pt):
	"""Make a GCT file from pt.

	Keyword arguments:
	pt -- a Pandas pivot table.
	"""

	# Remove mean column.
	pt.drop('mean', axis=1, inplace=True)
	# Convert pivot table to records. The index (tracking_id) is excluded.
	df = pd.DataFrame(pt.to_records(index=False))

	# Format column names: ('CAS', 0) becomes CAS_0, etc.
	# table contains the characters to remove from each column name.
	table = str.maketrans(dict.fromkeys("()' "))
	new_columns = []
	for c in df.columns:
		# Split by comma into list
		l = c.split(',')
		# Remove characters
		for i in range(len(l)):
			l[i] = l[i].translate(table)
		# Join using underscore
		new_columns.extend(['_'.join(l)])
	# Apply new column names.
	df.columns = new_columns

	# Rename gene to NAME.
	df.rename(columns={'gene_':'NAME'}, inplace=True)
	# Convert NAME column values to all caps.
	df['NAME'] = df['NAME'].str.upper()
	# Add Description column, equal to NAME column.
	df['Description'] = df['NAME']

	# Move NAME and Description to first two columns.
	cols = list(df.columns.values)
	# Remove NAME from list.
	cols.pop(cols.index('NAME'))
	# Remove Description from list.
	cols.pop(cols.index('Description'))
	# New dataframe with columns in the desired order.
	df = df[['NAME','Description']+cols]

	# If more than two groups, create pair-wise files.
	# ...

	# Insert column names.
	df = insert_row(df, list(df.columns.values))
	# Insert number rows, number of columns-2 (Name, Description).
	df = insert_row(df, [df.shape[0]-1, df.shape[1]-2])
	# Insert #1.2
	df = insert_row(df, ['#1.2'])

	return df

def write_output(file, data, delim, index=False, header=True):
	"""Write data to file.

	Keyword arguments:
	file -- output file to write data to
	data -- a Pandas dataframe or pivot table
	delim -- use this character to separate columns in file
	index -- set to True to include index column in file
	header -- write header to file, True or False.
	"""
	data.to_csv(path_or_buf=file, sep=delim, index=index, header=header)
	print_message('Output file: {0}'.format(file), False)

def df_info(df, message):
	"""Print info about df.

	Keyword arguments:
	df -- a Pandas dataframe
	message -- a string describing df
	"""
	print_message('\n{0}'.format(message), False)
	print_message('Rows: {0}'.format(df.shape[0]), False)
	print_message('Columns: {0}'.format(df.shape[1]), False)
	print_message('Empty cells: {0}'.format(df.isnull().sum(axis=0).sum()), False)

def get_gene(row, df_ed):
	"""Return gene if tracking_id exists in df_ed.

	Keyword arguments:
	row -- a dataframe row
	"""

	# row.name is the tracking_id
	result = df_ed.loc[(df_ed['test_id'] == row.name)]
	gene = 'NA'
	# If a result was returned, get first one then gene name.
	if result.shape[0] > 0:
		gene = result.iloc[0].gene
	return gene

def hidata_genes(df, file):
	"""Extract HIDATA genes from df. Write to file.

	Keyword arguments:
	df -- a Pandas dataframe
	file -- write results to this file
	"""
	print_message('\nHIDATA genes...', False)

	# Columns to remove from df.
	cols_rm = [
		'tracking_id',
		'class_code',
		'nearest_ref_id',
		'tss_id',
		'length',
		'coverage'
	]

	# Remove columns. axis 0 = rows, axis 1 = columns.
	df.drop(cols_rm, inplace=True, axis=1)
	print_message('Columns after removal: {0}'.format(df.shape[1]), False)

	# Column names containing '_status'
	status_cols = df.filter(like='_status').columns

	df_hidata_genes = pd.DataFrame()
	for col in status_cols:
		# Everything before _status
		col_name = col.split('_status')[0]
		# Get all rows containing HIDATA in col.
		# From that result, get the gene_short_name column.
		# Returns a series.
		sr_genes = df.loc[df[col] == 'HIDATA', 'gene_short_name']
		# Remove duplicates. Keep first occurance.
		sr_genes.drop_duplicates(keep='first', inplace=True)
		# Rename column from gene_short_name to col_name
		sr_genes.rename(col_name, inplace=True)
		# Sort alphabetically.
		sr_genes.sort_values(inplace=True)
		# Drop index (the row number). This minimizes NaN cells.
		sr_genes.reset_index(drop=True, inplace=True)
		# Add a new column to df_hidata_genes.
		df_hidata_genes = pd.concat([df_hidata_genes, sr_genes], axis=1)

	#print_message(df_hidata_genes, False)

	num_output_rows = df_hidata_genes.shape[0]
	print_message('Rows with HIDATA: {0}'.format(num_output_rows), False)
	if num_output_rows > 0:
		write_output(file, df_hidata_genes, '\t', False)

def read_input(file, delim):
	"""Read input file into a dataframe."""
	df = pd.read_csv(file, delimiter=delim)
	df_info(df, 'Input file: {0}'.format(file))
	return df

#########################
# Start script
#########################

def main():

	print_message('\n***\n* {0}\n* {1}\n***'.format(script_name,script_purpose), False)

	# Read input files into dataframes.
	df_fpkm = read_input(input_genes_fpkm, input_genes_fpkm_delim)
	df_rgt = read_input(input_genes_rg, input_genes_rg_delim)
	df_ed = read_input(input_gene_exp, input_gene_exp_delim)

	# From df_ed get rows where significant = yes
	df_ed = df_ed[df_ed["significant"] == 'yes']
	df_info(df_ed, 'gene_exp.diff significant=yes...')

	# Get the HIDATA genes.
	hidata_genes(df_fpkm, output_genes_fpkm)

	# Create pivot table.
	# Filter by status column, keeping rows not containing FAIL or HIDATA.
	pt = pd.pivot_table(df_rgt.loc[~df_rgt['status'].isin(['FAIL','HIDATA'])], index='tracking_id', columns=['condition','replicate'], values='FPKM')
	# Add a mean column.
	pt['mean'] = pt.mean(axis=1)
	df_info(pt, 'Pivot table & status=OK & mean column...')

	# Expressed genes: rows where the mean is greater than zero.
	pt = pt[pt['mean'] > 0]
	df_info(pt, 'Expressed genes (mean > 0)...')
	write_output(output_genes_rg_exp, pt, input_genes_rg_delim, index=True)

	# Significant genes.
	# Create column 'gene' set to NA or gene name. Apply function to each row.
	pt['gene'] = pt.apply(get_gene, axis=1, args=(df_ed,))
	# Remove rows where gene = NA
	pt = pt[pt.gene != 'NA']
	df_info(pt, 'Significant genes...')
	write_output(output_gene_exp_sig, pt, input_gene_exp_delim, index=True)

	# Make a GCT file from pt.
	gct_file = make_gct(pt)
	df_info(gct_file, 'Making GCT file...')
	write_output('gct.tsv', gct_file, '\t', header=False)

	print_message('\nDone!\n', False)

#########################
# Script entry point.
#########################

if __name__ == "__main__":
	main()



# Testing:

# import pandas as pd
# df_ed=pd.read_csv('/home/nick/Downloads/cuff_diff_summary/gene_exp.diff', delimiter='\t')
# df_ed = df_ed[df_ed["significant"] == 'yes']

# df = pd.read_csv('/home/nick/Downloads/cuff_diff_summary/genes.read_group_tracking', delimiter='\t')
# pt = pd.pivot_table(df.loc[~df['status'].isin(['FAIL','HIDATA'])], index='tracking_id', columns=['condition','replicate'], values='FPKM')
# pt['mean'] = pt.mean(axis=1)
# pt = pt[pt['mean'] > 0]

# def get_gene(row, df_ed):
# 	result = df_ed.loc[(df_ed['test_id'] == row.name)]
# 	gene = 'NA'
# 	if result.shape[0] > 0:
# 		gene = result.iloc[0].gene
# 	return gene

# pt['gene'] = pt.apply(get_gene, axis=1, args=(df_ed,))
# pt = pt[pt.gene != 'NA']

# pt.drop('mean', axis=1, inplace=True)

# pd.DataFrame(pt.to_records())
# result:
# tracking_id    ('CAS', 0)   ('CAS', 1)   ('CAS', 2)  ...
# 0    ENSMUSG00000000031.15     16.989900    17.526400    18.659900
