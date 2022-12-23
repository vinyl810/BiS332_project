from cmath import pi
from tqdm import tqdm
import math

nt = ["A", "G", "C", "T"]
INPUT_FILE = "input_data.txt"

from urllib.request import urlopen
import json

with urlopen(
    "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"
) as response:
    counties = json.load(response)

import pandas as pd
import plotly.express as px


def inputFile(file_name):
    f = open(file_name, "r")
    f.readline()
    fdata = f.readlines()
    f.close()
    npop = len(fdata)
    print("Parsing genetic distance input...")
    for i in range(npop):
        fdata[i] = fdata[i].strip().split("\t")
    nsnp = len(fdata[0]) - 5
    county = []
    county_idx = dict()
    county_pop = []
    county_fips = []
    for i in range(npop):
        if fdata[i][1] not in county:
            county_idx[fdata[i][1]] = len(county)
            county.append(fdata[i][1])
            county_pop.append(0)
            county_fips.append(fdata[i][4])
        county_pop[county_idx[fdata[i][1]]] += 1
    alFreq = [[dict() for j in range(nsnp)] for i in range(len(county))]
    for j in tqdm(range(nsnp)):
        for i in range(len(county)):
            for k in range(4):
                for l in range(4):
                    alFreq[i][j][nt[k] + nt[l]] = 0
        for i in range(npop):
            cou = county_idx[fdata[i][1]]
            allist = []
            if fdata[i][j + 5][0] == "N":
                for k in range(4):
                    alFreq[cou][j][nt[k] + fdata[i][j + 5][1]] += 1 / 4
            elif fdata[i][j + 5][1] == "N":
                for k in range(4):
                    alFreq[cou][j][fdata[i][j + 5][0] + nt[k]] += 1 / 4
            else:
                alFreq[cou][j][fdata[i][j + 5]] += 1

        for i in range(len(county)):
            for allele in alFreq[i][j].keys():
                alFreq[i][j][allele] /= county_pop[i]
    print("\033[2A\033[33C done.\033[1B")
    return county, county_fips, alFreq


def gen_dist(x_alf, y_alf):
    cavalli_sforza = 0
    assert len(x_alf) == len(y_alf)
    for i in range(len(x_alf)):
        original = cavalli_sforza
        for key in x_alf[i].keys():
            cavalli_sforza += math.sqrt(x_alf[i][key] * y_alf[i][key])
        if cavalli_sforza - original > 1:
            cavalli_sforza = original + 1
    cavalli_sforza = math.sqrt(2 * (1 - cavalli_sforza / len(x_alf))) / math.sqrt(2)
    return cavalli_sforza


def proc(county, alFreq):
    county_dist = [[0 for i in range(len(county))] for j in range(len(county))]
    for j in tqdm(range(len(county))):
        for i in range(len(county)):
            county_dist[i][j] = gen_dist(alFreq[i], alFreq[j])
    return county_dist


def indInput(ind_input):
    f = open(ind_input)
    ind_line = f.readline()
    f.close()
    return ind_line


def indInves(county, county_fips, alFreq, ind_line):
    ind_snp = ind_line.strip().split("\t")
    alF = [dict() for i in range(len(ind_snp))]
    for i in range(len(ind_snp)):
        for k in range(4):
            for l in range(4):
                alF[i][nt[k] + nt[l]] = 0
        if ind_snp[i][0] == "N":
            for k in range(4):
                alF[i][nt[k] + ind_snp[i][1]] = 1 / 4
        elif ind_snp[i][1] == "N":
            for k in range(4):
                alF[i][ind_snp[i][0] + nt[k]] = 1 / 4
        else:
            alF[i][ind_snp[i]] = 1
    ind_dist = []
    for cou in range(len(county)):
        ind_dist.append([gen_dist(alF, alFreq[cou]), county[cou], county_fips[cou]])
    ind_dist.sort(key=lambda lis: lis[0])

    return ind_dist


def create_county_dist_data(county, county_fips, county_data):
    cou_dist = []
    for cou in range(len(county)):
        cou_dist.append([county_data[cou], county[cou], county_fips[cou]])
    return cou_dist


def create_fig(data):
    df = pd.DataFrame(data, columns=["genetic_distance", "county", "fips"])
    fig = px.choropleth(
        df,
        geojson=counties,
        locations="fips",
        color="genetic_distance",
        color_continuous_scale="YlOrRd_r",
        range_color=(
            math.floor(df["genetic_distance"].min() * 100) / 100,
            math.ceil(df["genetic_distance"].max() * 100) / 100,
        ),
        scope="usa",
        labels={"genetic_distance": "genetic distance"},
    )
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
    return fig


def input_county_dist(cd_input):
    f = open(cd_input, "r")
    read_data = f.readlines()
    county_dist = []
    for i in range(len(read_data)):
        line = read_data[i].strip().split("\t")
        for j in range(len(line)):
            line[j] = float(line[j])
        county_dist.append(line)
    return county_dist


# county_dist.txt is generated when genetic_distance.py is executed, but not when imported.
if __name__ == "__main__":
    county, county_fips, alFreq = inputFile()
    county_dist = proc(county, alFreq)
    with open("county_dist.txt", "w") as external_file:
        for i in range(len(county_dist)):
            for j in range(len(county_dist[i])):
                print(f'{county_dist[i][j]:.15f}', end='\t', file = external_file)
            print('', file=external_file)
else:
    county, county_fips, alFreq = inputFile(INPUT_FILE)
    county_dist = input_county_dist("county_dist.txt")