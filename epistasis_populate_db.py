import pyodbc
import csv
import os
from functools import reduce
from decimal import Decimal
import argparse

parser = argparse.ArgumentParser(description="specifying argument for populating .gwas file")
# Specify the name of the table
parser.add_argument('table_name', action='store', help="the table to be populate the data")
parser.add_argument('-c', '--create', action='store_true', help="weather to create a new database", default=False)
# Specify the path of file folder
parser.add_argument('-p', '--path', action='store', help="path of the directory that contains .gwas file", default=".")
# Specify the database
parser.add_argument('-db', '--database', action='store', help="the destination database", default="Epistasis")
args = parser.parse_args()
path = args.path
database = args.database
table_name = args.table_name
create = args.create
'''
create a table to populate: 
    the current schema is fixed as (SNP1 varchar(20), SNP2 varchar(20), Trait varchar(50), Pvalue float)
'''
def create_table(db, table_name):
    global cursor
    query = "create table {!s} ".format(table_name) + "(" + \
            " SNP1 varchar(20)," \
            "SNP2 varchar(20)," \
            "Trait varchar(50)," \
            "Pvalue float, " \
			"CONSTRAINT PK PRIMARY KEY (SNP1, SNP2, Trait));" 
    #print(query)
    cursor.execute(query)
    cursor.commit()
    print("successfully create the table %s in database: %s" % (table_name, db))

def get_schema(db,table_name):
    global cursor
    query = " SELECT COLUMN_NAME" + \
            " FROM %s.INFORMATION_SCHEMA.COLUMNS \
            WHERE TABLE_NAME = " "\'%s\'" %(db, table_name)
    cursor.execute(query)
    result = cursor.fetchall()
    col_name = [tuple[0] for tuple in result]
    col_name = tuple(col_name)
    return col_name

def createConnect(server='PARKSLAB', database='Epistasis'):
    # username = 'myusername'
    # password = 'mypassword'
    cnxn = pyodbc.connect('DRIVER={SQL Server}' + \
                          ';SERVER=' + server + \
                          ';DATABASE=' + database + \
                          ';Trusted_Connection=Yes')
    return cnxn

# for splitting up gender attributes , commented out temporarily
# def readTrait_txt(traits):
#     os.chdir(path)
#     for trait in traits:
#         Gender = trait.split(".")[0]
#         with open (trait, 'r') as csvFile:
#             csvReader = csv.reader(csvFile, delimiter  = ',')
#             for row in csvReader:
#                 if Gender.strip() == 'female':
#                     query = "insert into table dbo.bxd_clinical_traits_linkage\
#                     (trait,RSID,female_lod)\
#                     values(trait,row[0], row[-1]);"
#                     cursor.execute(query)
#                 elif Gender.strip() == 'male':
#                     query = "insert into table dbo.bxd_clinical_traits_linkage\
#                     (trait,RSID,male_lod)\
#                     values(trait,row[0], row[-1]);"
#                     cursor.execute(query)
#                 else:
#                     print('I/O !')

if __name__ == '__main__':
    fileNames = []
    trait_Record = []

    for file in os.listdir(path):
        if file.endswith(".gwas"):
            fileNames.append(file)
    # connect to database
    cnxn = createConnect(database = database)
    print('connect the database: %s successfully!' %database)
    cursor = cnxn.cursor()
    # check the path
    if not os.path.isdir(path):
        print("path not exists")
        exit(1)

    os.chdir(path)
    # create table if specified
    if create:
        create_table(database, table_name)

    # search infomation from .gwas file through the directory specified by <path>
    for fileName in fileNames:
        # get the trait from file name
        trait = fileName.split(".gwas")[0]
        trait = trait.split("_")[:-1]
        trait = reduce((lambda x, y: x + "_" + y), trait)
        with open(fileName, 'r') as csvFile:
            csvReader = csv.reader(csvFile, delimiter='\t')
            # a dumb way to skip the header of a file. may re-write using f.readline () in the future
            count = 0;
            for row in csvReader:
                if count == 0:
                    count = 1
                    continue
                try:
                    query = "insert into dbo.{!s}".format(table_name) +\
                   "(SNP1, SNP2, Trait, Pvalue)" + \
                    " values ({!r}, {!r}, {!r}, {!s})".format(row[0], row[1], trait, str(Decimal(row[-1]))) 
                    cursor.execute(query)
                    # current setting is to commit every execution of the query
                    cursor.commit()

                # write errmsg if file I/O exception
                except ValueError as ex:
                    errmsg = "Warning: error in " + str(fileName) + "," + " row[-1] is: " + row[-1] + " type: " + str(
                        type(row[-1]))
                    print errmsg
                    print "current row:",
                    print row
                    f = open("Epistasis_{!r}_err.txt".format(table_name), "w")
                    f.write(errmsg + "\n")
                    #f.write(row) ## this was tripping up the code...not sure why.
                    f.close()

                except IndexError as iex:
                    print(row)

                except Exception as eex:
                    errmsg = "DUPLICATE KEY or WRONG VALUE in " + row[0] + " " + row[1]
                    print errmsg
                    print "current row:",
                    print row
                    f = open("Epistasis_{!r}_err.txt".format(table_name), "a+")
                    f.write(errmsg + "\n")
                    f.close()

    cursor.commit()
    print('Done!')
    cnxn.close()