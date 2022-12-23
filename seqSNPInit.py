from Bio import SeqIO
import json

fastaSeq = SeqIO.parse(open("GRCh38_latest_genomic.fna"), "fasta")
seqs = [
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
]

i = 0
for record in fastaSeq:
    if "NC_0000" in record.id:
        seqs[i] = seqs[i] + str(record.seq)
        i = i + 1
        print(record.id)

with open("seqs.json", "w") as json_file:
    json.dump(seqs, json_file)

print(seqs[22][101349768])
