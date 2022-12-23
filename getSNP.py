import requests, sys


def getSnpInfo(snp_id):
    server = "https://clinicaltables.nlm.nih.gov/api/snps/v3/search"
    ext = "?terms=rs" + str(snp_id)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    id = snp_id
    chromo = decoded[3][0][1]
    position = int(decoded[3][0][2]) + 1
    allele = decoded[3][0][3].split(",")[0].split("/")
    ancestor = allele[0]
    minor = allele[1]
    return [chromo, position, ancestor, minor, "rs" + str(id)]


# print(getSnpInfo(1921))
