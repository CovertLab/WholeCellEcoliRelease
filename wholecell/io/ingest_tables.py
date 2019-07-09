import os
import sys
import numpy as np

from pymongo import MongoClient
from wholecell.io.tablereader import TableReader

INDEX_COLUMNS = [
	'key',
	'column',
	'subcolumn',
	'subcolumn_index']

def tolist(a):
	'''If value is np.ndarray, convert it to a list'''

	value = a
	if isinstance(a, np.ndarray):
		value = a.tolist()
	return value

def create_indexes(table):
	'''Create all of the necessary indexes for the given table name.'''

	for column in INDEX_COLUMNS:
		table.create_index(column)

def ingest_tables(key, sim_out, url='localhost:27017', database='wcm'):
	'''
	Ingest all simulation output into a mongo database.

	Args:
	    key (str): Unique key for this simulation generation.
	    sim_out (str): Path to the simulation output directory we are ingesting.
	    url (str): URL to the mongodb instance.
	    database (str): Name of database to ingest simulation output into.
	'''

	client = MongoClient(url)
	db = getattr(client, database)
	db.attributes.create_index('table')

	table_names = [
		name
		for name in os.listdir(sim_out)
		if os.path.isdir(os.path.join(sim_out, name))]

	for table_name in table_names:
		print('reading table {}'.format(table_name))

		reader = TableReader(os.path.join(sim_out, table_name))
		table = getattr(db, table_name)
		create_indexes(table)

		attribute_names = reader.attributeNames()
		attributes = {
			attribute_name: reader.readAttribute(attribute_name)
			for attribute_name in attribute_names}

		db.attributes.insert_one({
			"table": table_name,
			"attributes": attributes})

		column_names = reader.columnNames()
		subcolumns = attributes.get('subcolumns', {})

		for column_name in column_names:
			print('reading column {}'.format(column_name))

			column = reader.readColumn(column_name)
			subcolumn_header = subcolumns.get(column_name)

			if subcolumn_header:
				subcolumn_names = attributes.get(subcolumn_header)

				if subcolumn_names is None:
					raise ValueError('table {} column {} has missing subcolumn names (probably a typo): {}', table_name, column_name, subcolumns)

				for index, subcolumn_name in enumerate(subcolumn_names):
					print('reading subcolumn {} - {}'.format(index, subcolumn_name))

					document = {
						'key': key,
						'column': column_name,
						'subcolumn': subcolumn_name,
						'subcolumn_index': index,
						'data': tolist(column[:, index])}

					table.insert_one(document)
			else:
				document = {
					'key': key,
					'column': column_name,
					'subcolumn': column_name,
					'subcolumn_index': 0,
					'data': tolist(column)}

				table.insert_one(document)

class DBTableReader(object):
	'''
	This class emulates wholecell.io.tablereader.TableReader but provides results
	from a mongo database instead.
	'''

	def __init__(self, key, path, url='localhost:27017', database='wcm'):
		'''
		Initialize the DBTableReader object with a given unique key and "path" to the
		table of interest. In TableReader, this path is a filesystem path, but here we
		use only the `os.path.basename(path)` to identify the mongo collection.

		Args:
	        key (str): Unique key for this simulation generation.
		    path (str): Identifer for table.
		'''

		self.key = key
		self.table_name = os.path.basename(path)
		self.client = MongoClient(url)
		self.db = getattr(self.client, database)
		self.table = getattr(self.db, self.table_name)
		self.attributes = self.db.attributes.find_one({'table': self.table_name})['attributes']

	def readColumn(self, column_name, indices=None):
		'''
		Return all data associated with `column_name` for this table.

		Args:
		    column_name (str): Which column to retrieve.
		    indices (list(int)): Indexes of subcolumns to retrieve (all if `None`)
		'''

		subcolumn_header = self.attributes.get('subcolumns', {}).get(column_name)
		column = []

		if subcolumn_header:
			if indices:
				results = self.table.find({
					'key': self.key,
					'column': column_name,
					'subcolumn_index': {'$in': indices}})
			else:
				results = self.table.find({
					'key': self.key,
					'column': column_name}).sort('subcolumn_index', 1)

			column = np.transpose(np.array([result['data'] for result in results]))

		else:
			result = self.table.find_one({
				'key': self.key,
				'column': column_name})

			if result:
				column = np.array(result['data'])

		return column

if __name__ == '__main__':
	key = sys.argv[1]
	sim_out = sys.argv[2]
	ingest_tables(key, sim_out)
