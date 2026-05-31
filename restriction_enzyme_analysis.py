from Bio import SeqIO
import os
import re

sequence = ""
seq_name = ""
enzymes = {}
fasta_path = r"C:\Users\12262\PyCharmMiscProject\testseq.fasta"
enzymes_path = r"C:\Users\12262\PyCharmMiscProject\enzymes.txt"

for record in SeqIO.parse("testseq.fasta", "fasta"):
    sequence = str(record.seq)
    seq_name = record.id

with open("enzymes.txt", "r") as f:
    for line in f:
        enzyme = line.strip().split(";")
        enzymes[enzyme[0]] = enzyme[1]

fasta_filename = os.path.basename(fasta_path)
enzymes_filename = os.path.basename(enzymes_path)

msg = f'Restriction enzyme analysis of sequence from file {fasta_filename}.'
line_width = len(msg)

print(msg)
print(f'Cutting with enzymes found in file {enzymes_filename}.')

print("-" * line_width)

print(f'Sequence name: {seq_name}')
print(f'Sequence is {len(sequence)} bases long.')

for enzyme, cutting_site in enzymes.items():

    restriction_site = cutting_site.replace("^", "")
    matches = re.findall(restriction_site, sequence)
    site_count = len(matches)
    fragment_count = site_count + 1
    modified_sequence = sequence.replace(restriction_site, cutting_site)
    fragments = modified_sequence.split("^")

    if site_count > 0:
        print("-" * line_width)
        print(f'There are {site_count} cutting sites for {enzyme}, cutting at {cutting_site}')
        print(f'There are {fragment_count} fragments:')
        print()
        index = 1
        for fragment in fragments:
            frag_len = len(fragment)
            print(f'length: {frag_len}')
            print(f'{index:<8}', end="")
            for nt in range(0, frag_len, 10):
                print(fragment[nt:nt + 10], end=" ")
            print()
            index += len(fragment)
        print()
    else:
        print("-" * line_width)
        print(f'There are no sites for {enzyme}.')

