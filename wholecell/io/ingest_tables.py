import os
from pymongo import MongoClient
from wholecell.io.tablereader import TableReader

def ingest_tables(key, sim_out, url='localhost:27017', database='wcm'):
	client = MongoClient(url)
	table_names = os.listdir(sim_out)

	for table_name in table_names:
		reader = TableReader(os.path.join(sim_out, table_name))

		attribute_names = reader.attributeNames()
		attributes = {
			attribute_name: reader.readAttribute(attribute_name)
			for attribute_name in attribute_names}

		column_names = reader.columnNames()
		
