#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Selects best non-overlapping hits from a homology search hit table \
                                 (e.g. from tblastn), based on score or evalue')

parser.add_argument('--table',help = 'table of hits')
parser.add_argument('--sep', default = 'tab', help = 'table field delimiter, accepts "comma" and "tab" (default = tab')
parser.add_argument('--sortby',default = 'evalue', help = 'Criterion for "best" hit, accepts "evalue" \
                    (assumed to be field 11, minimizes) or "score" (assumed to be field 12, maximizes)')
parser.add_argument('--input_type', default = "blast_table", help = 'input file type- accepts "blast_table" (default) or \
                    "gtf"')
parser.add_argument('--output', default = "non-overlapping", help = 'specify "non-overlapping" (default) or "overlapping"')

args = parser.parse_args()
if args.sep == "comma":
    delimiter = ','
elif args.sep == "tab":
    delimiter = "\t"

sort_list = []
if args.input_type == "gtf":
    sortby_field = 3
elif args.sortby == 'evalue':
    sortby_field = 10
elif ags.sortby == 'score':
    sortby_field = 11
else:
    raise 'incorrect argument for "sortbby"'
    exit()

if args.input_type == "blast_table":
    coord1,coord2 = 8,9
elif args.input_type == "gtf":
    coord1,coord2 = 3,4

for line in open(args.table):
    fields = line.replace('\r','').replace('\n','').split(delimiter)
    sort_list.append([float(fields[sortby_field]),fields])

sort_list.sort()
if args.sortby == 'score':
    sort_list.reverse()

chunksize = 100000
lines_list = []
coords_dict = {}

for hit in sort_list:
    fields = hit[1]
    start, stop = min(int(fields[coord1]), int(fields[coord2])),max(int(fields[coord1]), int(fields[coord2]))
    if args.input_type == "blast_table":
        locus = fields[1]
    elif args.input_type == "gtf":
        locus = fields[0]
    coords = (start,stop)
    start_key = int(start / chunksize)
    stop_key = int(stop / chunksize)
    keys_list = [(locus,start_key)]
    if start_key != stop_key:
        for i in range(stop_key - start_key):
            keys_list.append((locus,keys_list[-1][-1] + 1))
    keep = True
    for key in keys_list:
        if key in coords_dict:
            for occupied_coord in coords_dict[key]:
                if occupied_coord[0] <= start <= occupied_coord[1] or \
                occupied_coord[0] <= stop <= occupied_coord[1] or \
                start < occupied_coord[0] < stop:
                    keep = False
                    if args.output == "overlapping":
                        print(delimiter.join(fields))
                    break
    if keep:
        for key in keys_list:       
            if not key in coords_dict:
                coords_dict[key] = [coords]
            else:
                coords_dict[key].append(coords)
        lines_list.append(delimiter.join(fields))

if args.output == "non-overlapping":
    print("\n".join(lines_list))


