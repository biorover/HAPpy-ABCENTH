#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas

#All of the random small functions used in other magot programs


import copy
import io


def tab2fasta(tab_file):
    tabs = ensure_file(tab_file)
    newlines = []
    for line in tabs:
        fields = line.replace('\n','').replace('\r','').split('\t')
        newlines.append('>'+fields[0]+'\n'+fields[1])
    return "\n".join(newlines)


def fasta2tab(fasta_file):
    fasta = ensure_file(fasta_file)
    newlines = []
    for line in fasta:
        if line[0] == '>':
            newlines.append(line[1:].replace('\n','').replace('\r','') + '\t')
        else:
            newlines[-1] = newlines[-1] + line.replace('\n','').replace('\r','')
    return '\n'.join(newlines)


def ensure_file(potential_file):
    """takes an input that can be either a file location, opened file, or string, and returns an opened file or file like object.
    Used to make downstream applications robust to different input types"""
    if potential_file == None:
        return None
    elif type(potential_file).__name__ == "file":
        return potential_file
    else:
        try:
            return open(potential_file)
        except IOError:
            return io.StringIO(potential_file)


def read_to_string(potential_file):
    """Tries to read "potential_file" into string- first as file location,
    then file, then string."""
    try:
        output_string = open(potential_file).read().replace('\r','')
    except IOError:
        try:
            output_string = potential_file.read().replace('\r','')
        except AttributeError:
            output_string = potential_file.replace('\r','')
    return output_string





