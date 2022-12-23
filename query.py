### 0. Libraries
##  0.1. sqlalchemy: Sets connection to sql database.
##  0.2. pg8000: A driver for postgresql dialect. (Version 1.16.5 WORKS)
##  0.3. pandas: Easy translation from pandas dataframe to sql db.
import sqlalchemy as db
from sqlalchemy import create_engine
import pg8000  # IF ERROR OCCURS, TRY pip install pg8000==1.16.5
import pandas as pd
import pandas.io.sql as psql

### 1. Sets postgresql connection
##  1.1. postgresql server connection info
login = "u20210351"  # YOUR login ID
password = "kiyomiday1"  # YOUR login pw.
host = "biostar.kaist.ac.kr:5432"  # HOST for sql server. DO NOT MODIFY

##  1.2. Create connection with given information
engine = create_engine(
    "postgresql+pg8000://" + login + ":" + password + "@" + host + "/" + login
)
conn = engine.connect()

### 2. Interfaces
##  2.1. sendQuery(query): sends query to sql db
#        query(string): query in sql language to get the information
def sendQuery(query):
    try:
        df = psql.read_sql_query(query, conn)
        return df
    except db.exc.ResourceClosedError:
        return None
    except db.exc.ProgrammingError as progError:
        return "[BioJoin.ProgrammingError]: " + progError.orig.args[0]["M"]
    except db.exc.IntegrityError as intError:
        return "[BioJoin.IntegrityError]: " + intError.orig.args[0]["M"]


##  2.3. mkTableFromTsv(file, tableName, columnNames): makes table from given tab-seperated file and commits to db
#        file(string)
#        tablename(string)
#        columnNames(array)
def mkTableFromTsv(file, tableName, columnNames):
    data = pd.read_csv(file, sep="\t", header=0)
    data.columns = columnNames
    data.to_sql(name=tableName, con=engine, index=False)
