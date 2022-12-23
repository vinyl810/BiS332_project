#import requirements
from query import *
import pandas as pd
import sqlalchemy
from tqdm import tqdm
from sqlalchemy import create_engine

if __name__ == '__main__':
    column = ['tree_id', 'county', 'state', 'county_state', 'fips', 'snp_genotype']
    f = open("input_data.txt")
    f.readline()
    lines = f.readlines()
    to_insert = []
    query = "DELETE FROM snp_genotype_550 *;"
    sendQuery(query)
    for line in tqdm(lines):
        data = line.strip().split('\t')
        snp_genotype = ''.join(data[5:])
        to_insert = [data[0], data[1], data[2], data[3], data[4], snp_genotype]
        query = "INSERT INTO snp_genotype_550 VALUES" + str(tuple(to_insert)) + ";"
        sendQuery(query)
    f.close()