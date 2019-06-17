import os
import sys
import numpy as np

from pymongo import MongoClient
from wholecell.io.tablereader import TableReader

def tolist(a):
	value = a
	if isinstance(a, np.ndarray):
		value = a.tolist()
	return value

def ingest_tables(key, sim_out, url='localhost:27017', database='wcm'):
	client = MongoClient(url)
	db = getattr(client, database)

	table_names = [
		name
		for name in os.listdir(sim_out)
		if os.path.isdir(os.path.join(sim_out, name))]

	for table_name in table_names:
		print('reading table {}'.format(table_name))

		reader = TableReader(os.path.join(sim_out, table_name))
		table = getattr(db, table_name)

		attribute_names = reader.attributeNames()
		attributes = {
			attribute_name: reader.readAttribute(attribute_name)
			for attribute_name in attribute_names}

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

if __name__ == '__main__':
	key = sys.argv[1]
	sim_out = sys.argv[2]
	ingest_tables(key, sim_out)
