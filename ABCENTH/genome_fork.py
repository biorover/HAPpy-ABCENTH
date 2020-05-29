#!/usr/bin/python
#MAGOT: a functional and simple library for genomic analysis
#Sean McKenzie and Nelson Salinas
#See license on github (https://github.com/biorover/MAGOT/blob/master/LICENSE)
#Comments, feature requests, and friendly encouragement can be offered via github (https://github.com/biorover/MAGOT/issues)

#Known issues:
#Currently UTRs are not handled well. Most ways to import gff3 with UTRs encoded by exon feature regions which don't
#overlap with CDS feature regions will result in the UTRs being lost (because exon features are usually summarily ignored).
#Additionally, even if you manage to import UTRs, most output will likely ignore them. Bug us about it on a github ticket
#and we may well fix it.
#
#Importing GFFs still quite slow, likely having to do with ID assignment. Will try to fix soon.


import copy
import StringIO
from magot_smallfuncs import *

import numpy

verbose = True


def starjunc2gff(starjunc, output = 'string'):
    """parses junction file from STAR RNAseq aligner and returns an apollo-style gff3. Output can be "string", "list", or "print"."""
    if output != 'print':
        outlines = []
    junc_file = ensure_file(starjunc)
    strands = ['.','+','-']
    counter = 1
    for line in junc_file:
        ID = "StarAlignment_" + str(counter)
        counter = counter + 1
        atts = line.split('\t')
        stop1 = int(atts[1]) - 1
        start1 = stop1 - 20
        start2 = int(atts[2]) + 1
        stop2 = start2 + 20
        gffline1 = '\t'.join([atts[0],'star','match',str(start1),str(stop2),atts[6], strands[int(atts[3])],'.','ID=' + ID])
        gffline2 = '\t'.join([atts[0],'star','match_part',str(start1),str(stop1),atts[6], strands[int(atts[3])],'.','ID=' + ID + '-part1;Parent=' + ID])
        gffline3 = '\t'.join([atts[0],'star','match_part',str(start2),str(stop2),atts[6], strands[int(atts[3])],'.','ID=' + ID + '-part2;Parent=' + ID])
        if output == 'print':
            print gffline1
            print gffline2
            print gffline3
        else:
            outlines.append(gffline1)
            outlines.append(gffline2)
            outlines.append(gffline3)
    if output == "string":
        return "\n".join(outlines)
    elif output == "list":
        return outlines
    elif output == "print":
        pass
    else:
        print "invalid output format"


def apollo2genome(apollo_gff):
    apollo_list = read_to_string(apollo_gff).split('>')
    gff3 = apollo_list[0]
    fasta = '>' + apollo_list[1]
    return Genome(fasta,gff3,annotation_format='gff3')


def vulgar2gff(vulgarlist, feature_types=['match','match_part'],source='exonerate'):
    """takes vulgar alignment list (e.g. vulgarstring.split() ) and outputs gff lines. For eventual use with a read_exonerate function"""
    #sets variables to be added to gff lines
    qname = vulgarlist[0] + '-against-' + vulgarlist[4]
    qstart = vulgarlist[1]
    qend = vulgarlist[2]
    qstrand = vulgarlist[3]
    tname = vulgarlist[4]
    tstart = vulgarlist[5]
    tend = vulgarlist[6]
    tstrand = vulgarlist[7]
    score = vulgarlist[8]
    vulgartrips = vulgarlist[9:]
    addfeat = False
    IDnum = 1
    if tstrand == "+":
        tposition = int(tstart) + 1
    else:
        tposition = int(tstart)
        tend = str(int(tend)+1)
    #makes top level gff line for hit
    gfflines = ["\t".join([tname,source,feature_types[0],str(tposition),tend,score,tstrand,'.','ID='+qname])]
    #iterates over matches within hit to make bottom level feature gff lines
    for i in range(len(vulgartrips)):
        field = vulgartrips[i]
        if i % 3 == 0:
            if field == 'M' or field == 'S' or field == 'G' or field == 'F':
                if not addfeat:
                    addfeat = True
                    line_to_add = [tname,source,feature_types[1]]
                    coords = [str(tposition)]
                else:
                    pass
            elif addfeat:
                gfflines.append('\t'.join(line_to_add+[str(min(coords)),str(max(coords)),score,tstrand,'.',
                                                       'ID='+qname+'_'+feature_types[1]+str(IDnum)+';Parent='+qname]))
                IDnum = IDnum + 1
                addfeat = False
            else:
                pass
        if i % 3 == 2:
            if tstrand == "+":
                tposition = tposition + int(field)
            elif tstrand == "-":
                tposition = tposition - int(field)
            if addfeat == True:
                if tstrand == '+':
                    coords.append(str(tposition-1))
                elif tstrand == '-':
                    coords.append(str(tposition+1))
    if addfeat:
        gfflines.append('\t'.join(line_to_add+[str(min(coords)),str(max(coords)),score,tstrand,'.',
                                                       'ID='+qname+'_'+feature_types[1]+str(IDnum)+';Parent='+qname]))
    return '\n'.join(gfflines)


def read_exonerate(exonerate_output,annotation_set_to_modify = None,feature_types = ['match','match_part']):
    #Reads an exonerate file with a little help from vulgar2gff. 
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    exonerate_lines = ensure_file(exonerate_output)
    gfflines = []
    IDdic = {}
    qname = ""
    tname = ""
    for original_line in exonerate_lines:
        line = original_line.replace('\r','').replace('\n','')
        if line[:16] == "         Query: ":
            qname = line[16:]
        elif line[:16] == "        Target: ":
            tname = line[16:].replace(':[revcomp]','').replace('[revcomp]','')
            if tname[-1] == " ":
                tname = tname[:-1]
        elif line[:16] == "         Model: ":
            exonerate_model = line.split()[1].split(':')[0]
        elif line[:8] == "vulgar: ":
            vulgar_line_list = line[8:].split()
            #trying to makesure IDs are unique
            vulgar_line_list[0] = qname
            vulgar_line_list[4] = tname
            ID = vulgar_line_list[0] + '-against-' + vulgar_line_list[4]
            if ID in IDdic:
                vulgar_line_list[0] = vulgar_line_list[0] + "-" + str(IDdic[ID])
                IDdic[ID] = IDdic[ID] + 1
            else:
                IDdic[ID] = 1
            gfflines.append(vulgar2gff(vulgar_line_list,feature_types = feature_types,source = exonerate_model))
    read_gff("\n".join(gfflines), annotation_set_to_modify = annotation_set)
    if annotation_set_to_modify == None:
        return annotation_set
    
    
def write_longform_gff(annotation_set,keep_UTR_features = False):
    """returns gff string formated for compatability with Apollo genome annotation """
    gfflines = []
    fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
    #adds matches to gff
    if 'match' in annotation_set.__dict__:
        for match in annotation_set.match:
            match_obj = annotation_set.match[match]
            newline_list = []
            for field in fields:
                try:
                    newline_list.append(str(eval('match_obj.' + field)))
                except:
                    newline_list.append('.')
            attribute_list = ['ID=' + match_obj.ID]
            for attribute in match_obj.__dict__:
                if not attribute in fields+['annotation_set','parent','child_list','ID']:
                    attribute_list.append(attribute + '=' + eval('match_obj.' + attribute))
            newline_list.append(';'.join(attribute_list))
            gfflines.append('\t'.join(newline_list))
            for match_part in match_obj.child_list:
                match_part_obj = annotation_set[match_part]
                newline_list = []
                for field in fields:
                    try:
                        newline_list.append(str(eval('match_part_obj.' + field)))
                    except:
                        newline_list.append('.')
                attribute_list = ['ID=' + match_part_obj.ID,'Parent=' + match_part_obj.parent]
                for attribute in match_part_obj.__dict__:
                    if not attribute in fields+['annotation_set','parent','child_list','ID','coords']:
                        attribute_list.append(attribute + '=' + eval('match_part_obj.' + attribute))
                newline_list.append(';'.join(attribute_list))
                gfflines.append('\t'.join(newline_list))
        #adds genes to gff
    if 'gene' in annotation_set.__dict__:
        for gene in annotation_set.gene:
            gene_obj = annotation_set.gene[gene]
            newline_list = []
            for field in fields:
                try:
                    newline_list.append(str(eval('gene_obj.' + field)))
                except:
                    newline_list.append('.')
            attribute_list = ['ID='+gene_obj.ID]
            for attribute in gene_obj.__dict__:
                if not attribute in fields+['annotation_set','parent','child_list','ID']:
                    attribute_list.append(attribute + '=' + eval('gene_obj.' + attribute))
            newline_list.append(';'.join(attribute_list))
            gfflines.append('\t'.join(newline_list))
            for gene_child in gene_obj.child_list:
                gene_child_obj = annotation_set[gene_child]
                if gene_child_obj.feature_type == 'transcript':
                    transcript_obj = gene_child_obj
                    newline_list = []
                    for field in fields:
                        try:
                            newline_list.append(str(eval('transcript_obj.'+field)).replace('transcript','mRNA'))
                        except:
                            newline_list.append('.')
                    attribute_list = ['ID='+transcript_obj.ID,'Parent='+transcript_obj.parent]
                    for attribute in gene_obj.__dict__:
                        if not attribute in fields+['annotation_set','parent','child_list','ID']:
                            attribute_list.append(attribute + '=' + eval('transcript_obj.' + attribute))
                    newline_list.append(';'.join(attribute_list))
                    gfflines.append('\t'.join(newline_list))
                    exondict = {}
                    CDS_UTR_dict = {}
                    for transcript_child in transcript_obj.child_list:
                        transcript_child_obj = annotation_set[transcript_child]
                        line_base_list = []
                        for field in fields:
                            try:
                                line_base_list.append(str(eval('transcript_child_obj.'+ field)))
                            except:
                                line_base_list.append('.')
                        exon_attributes = 'ID=' + transcript_child_obj.ID + '-exon;Parent=' + transcript_child_obj.parent
                        transcript_child_attribute_list = ['ID=' + transcript_child_obj.ID, 'Parent=' + transcript_child_obj.parent]
                        for attribute in transcript_child_obj.__dict__:
                            if not attribute in fields+['annotation_set','parent','child_list','ID','coords']:
                                transcript_child_attribute_list.append(attribute + '=' + eval('transcript_child_obj.' + attribute))
                        transcript_child_attributes = ';'.join(transcript_child_attribute_list)
                        exondict[transcript_child_obj.coords] = '\t'.join(line_base_list).replace('CDS','exon').replace('UTR','exon') + '\t' + exon_attributes
                        CDS_UTR_dict[transcript_child_obj.coords] = '\t'.join(line_base_list) +'\t' + transcript_child_attributes
                    exondict_list = list(exondict)
                    exondict_list.sort()
                    CDS_UTR_dict_list = list(CDS_UTR_dict)
                    CDS_UTR_dict_list.sort()
                    #merges adjacent exons from abbutting UTRs and CDSs and writes exon gff lines
                    for i in range(len(exondict_list) - 1):
                        if exondict_list[i][1] + 1 == exondict_list[i+1][0]:
                            del exondict[exondict_list[i]]
                            new_exon_list = exondict[exondict_list[i+1]].split('\t')
                            exondict[exondict_list[i+1]] = '\t'.join(new_exon_list[:3]+[str(exondict_list[i][0]),str(exondict_list[i+1][1])]+new_exon_list[5:])
                        else:
                            gfflines.append(exondict[exondict_list[i]])
                    gfflines.append(exondict[exondict_list[-1]])
                    for i in CDS_UTR_dict_list:
                        if CDS_UTR_dict[i].split('\t')[2] == 'CDS' or keep_UTR_features:
                            gfflines.append(CDS_UTR_dict[i])
                else:
                    print 'ERROR: currently only accepts AnnotationSets with gene format CDS/UTR -> transcript -> gene'
                    break
    return '\n'.join(gfflines)

def write_gff(annotation_set, gff_format = "simple gff3"):
    gff_lines = []
    for attribute in annotation_set.__dict__:
        if type(annotation_set.__dict__[attribute]) == dict:
            if len(annotation_set.__dict__[attribute]) > 0:
                if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "ParentAnnotation" or annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "BaseAnnotation":
                    if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].parent == None:
                        for annotationID in annotation_set.__dict__[attribute]:
                            annotation = annotation_set.__dict__[attribute][annotationID]
                            gff_lines.append(annotation.get_gff(gff_format))
    return '\n'.join(gff_lines)

def write_bed(annotation_set, columns = 12):
    bed_lines = []
    for attribute in annotation_set.__dict__:
        if type(annotation_set.__dict__[attribute]) == dict:
            if len(annotation_set.__dict__[attribute]) > 0:
                if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "ParentAnnotation" or annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].__class__.__name__ == "BaseAnnotation":
                    if annotation_set.__dict__[attribute][list(annotation_set.__dict__[attribute])[0]].parent == None:
                        for annotationID in annotation_set.__dict__[attribute]:
                            annotation = annotation_set.__dict__[attribute][annotationID]
                            bed_lines.append(annotation.get_bed(columns=columns))
    return '\n'.join(bed_lines)
    


def read_bed(bed, annotation_set_to_modify = None, parents_hierarchy = ['region_set'], base_features = 'region'):
    bed_file = ensure_file(bed)
    ID_set = set()
    new_ID_dict = {}
    if annotation_set_to_modify:
        annotation_set = annotation_set_to_modify
        for attribute in annotation_set.__dict__:
            if type(annotation_set.__dict__[attribute]) == dict:
                for feature_id in annotation_set.__dict__[attribute]:
                    ID_set.add(feature_id)
                    new_ID_dict[feature_id] = 2
    else:
        annotation_set = AnnotationSet()
    possible_fields = ['seqid','start','stop','name','score','strand','thickStart','thickEnd','itemRgb']
    for line in bed_file:
        if '\t' in line:
            fields = line.split('\t')
        else:
            fields = line.split()
        if len(fields) < 3:
            print "Error: line with fewer than three columns. First three BED columns are manditory."
            print line
            continue
        #Figures out ID
        if len(fields) > 3:
            ID = fields[3]
        else:
            ID = fields[0] + ':' + fields[1] + '-' + fields[2]
        bol = False
        while ID in ID_set:
            bol = True
            ID = ID + '-' + str(new_ID_dict[ID])
            new_ID_dict[ID] = new_ID_dict[ID] + 1
        if not bol:
            new_ID_dict[ID] = 2
        ID_set.add(ID)
        #create parents
        for parent_index in range(len(parents_hierarchy)):
            parent_type = parents_hierarchy[parent_index]
            if parent_index == 0:
                child_list = []
            else:
                child_list = [parent_ID]
            if parent_index == len(parents_hierarchy) - 1:
                parent = None
                parent_ID = ID
            else:
                if parent_index == len(parents_hierarchy) - 2:
                    parent = ID
                else:
                    parent = ID + '-' + parents_hierarchy[parent_index + 1]
                parent_ID = ID + '-' + parent_type
            annotation_set[ID] = ParentAnnotation(parent_ID, fields[0], parent_type, child_list = child_list, parent = None, strand = ".", annotation_set = None, other_attributes = {})
    
    



def read_gff(gff,annotation_set_to_modify = None, base_features = ['CDS','match_part','similarity','region'],features_to_ignore = ['exon'],
    gff_version = "auto", parents_hierarchy = [], features_to_replace = [], IDfield = "ID", parent_field = "Parent", presets = None):
    """Assumptions of file:    
    Relatively consistent, nothing crazy like mixed gff3 and gff2
    
    If a parent's children are base features, they should NOT overlap. E.g. if your file has exons and CDS, only
    give CDS as base_features and give exons as features_to_ignore.
    
    Although parent feature types may be used multiply in hierarchies, any given base feature
    will only be present in one hierarchy
    
    If parent features are not represented by lines (e.g. minimal gff2 format or cegma), they
    can be inferred from defline attributes given in parents_hierarchy. It is expected that the order they
    are given is consistent with the hierarchy (e.g. ['mRNA','gene']). If an underscore is used in
    this defline, the word before the underscore is used as the parent feature type. The "Parent"
    defline attribute is reserved for the imediate parent of the feature- if not given, this will be assumed to
    be the lowest-level parent specified in defline.
    """
    replace_dict = {"\n":"","\r":""}
    presets_dict = {}
    presets_dict["augustus"] = "features_to_ignore = ['gene','transcript','stop_codon','terminal','internal','initial','intron',\
        'start_codon','single']\nparent_field = None\nparents_hierarchy = ['transcript_id','gene_id']\nIDfield = None\ngff_version=2"
    presets_dict["RepeatMasker"] = "IDfield = None\nparent_field = None\nparents_hierarchy = ['match']\nfeatures_to_replace = [['similarity','match_part']]"
    presets_dict["CEGMA"] = "features_to_replace = [['First','CDS'],['Internal','CDS'],['Terminal','CDS'],['Single','CDS']]\n"
    if presets in presets_dict:
        exec(presets_dict[presets])
    version = gff_version
    gff_file = ensure_file(gff)
    for feature in features_to_replace:
        replace_dict["\t" + feature[0] + "\t"] = "\t" + feature[1] + "\t"
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    #this dictionary and this list helps generate names if features passed with identical ID fields
    generate_new_ID_dict = {}
    allready_present_IDs = set()
    if annotation_set_to_modify:
        for attribute in annotation_set.__dict__:
            if type(annotation_set.__dict__[attribute]) == dict:
                for feature_id in annotation_set.__dict__[attribute]:
                    allready_present_IDs.add(feature_id)
    #figures out which feature types are base annotations
    base_dict = {}
    #here we go folks
    for original_line in gff_file:
        if original_line[0] != "#" and original_line.count('\t') == 8:
            line = original_line
            for string in replace_dict:
                line = line.replace(string,replace_dict[string])
            fields = line.split('\t')
            #Tries to detect GFF type if no presets or if presets are vague
            #
            if version == "auto":
                if "=" in fields[8]:
                    version = 3
                elif " " in fields[8]:
                    version = 2
                    if not IDfield == None:
                        if not " " + IDfield + " " in " " + fields[8].replace(';',' ') and parents_hierarchy == []:
                            IDfield = None
                            parent_field = None
                            if "gene_id" in fields[8] and "transcript_id" in fields[8]:
                                parents_hierarchy = ['transcript_id','gene_id']
                            elif "gene_id" in fields[8]:
                                parents_hierarchy = ['transcript_id','gene_id']
                else:
                    version = None
                    parents_hierarchy = ['transcript_id','gene_id']
                    IDfield = None
                    parent_field = None
            #Sets up basic variables- same for all gff types
            #
            ID = None
            parent = None
            other_attributes = {}
            seqid = fields[0]
            other_attributes['source'] = fields[1]
            feature_type = fields[2]
            if feature_type in features_to_ignore:
                continue
            coords = [int(fields[3]),int(fields[4])]
            coords.sort()
            coords = tuple(coords)
            try:
                other_attributes['score'] = float(fields[5])
            except ValueError:
                pass
            strand = fields[6]
            if fields[7] in ['0','1','2']:
                other_attributes['phase'] = int(fields[7])
            defline_dict = {}
            #Parsing the defline
            #
            for defline_field in fields[8].split(';'):
                if defline_field != "":
                    if version == None:
                        defline_dict['gene_id'] = defline_field
                    elif version == 2:
                        if '"' in defline_field:
                            defline_dict[defline_field.split()[0]] = defline_field.split('"')[1]
                        else:
                            try: defline_dict[defline_field.split()[0]] = defline_field.split()[1]
                            except:
                                print defline_field
                                return None
                    elif version == 3:
                        defline_dict[defline_field.split('=')[0]] = defline_field.split('=')[1]
            #Tries to figure out parent
            if parent_field != None:
                if parent_field in defline_dict:
                    parent = defline_dict[parent_field]
            elif parents_hierarchy != []:
                for parent_type in parents_hierarchy:
                    if parent_type in defline_dict:
                        parent = defline_dict[parent_type]
                        break
            #Assigns ID
            if IDfield != None:
                ID = defline_dict[IDfield]
            elif parent != None:
                ID = parent + '-' + feature_type
            else:
                ID = seqid + '-' + feature_type + fields[3]            
            #checks if ID already in annotation_set and adjusts ID if necessary
            if ID in allready_present_IDs:
                try:
                    generate_new_ID_dict[ID] = generate_new_ID_dict[ID] + 1
                    ID = ID + "-" + str(generate_new_ID_dict[ID])
                except KeyError:
                    generate_new_ID_dict[ID] = 2
                    ID = ID + '2'
            allready_present_IDs.add(ID)
            #If no parent assigned the first time, but one should be created, sets parent name
            if parent == None and parents_hierarchy != []:
                parent = ID + '-' + parents_hierarchy[0]
            #Creates parents if necessary, adds to parent if exists
            if parent != None:
                child_to_assign = ID
                for parent_feature_index in range(len(parents_hierarchy)):
                    parent_feature = parents_hierarchy[parent_feature_index]
                    if parent_feature in defline_dict:
                        parent_feature_ID = defline_dict[parent_feature]
                        parent_feature_type = parent_feature.split('_')[0]
                    #guesses a parent id and creates when not in defline dict
                    else:
                        parent_feature_type = parent_feature.split('_')[0]
                        if parents_hierarchy[-1] in defline_dict:
                            parent_feature_ID = defline_dict[parents_hierarchy[-1]] + "-" + parent_feature_type
                        else:
                            parent_feature_ID = ID + "-" + parent_feature_type
                        #resets parent if this should be parent
                        if parent_feature_index == 0:
                            parent = parent_feature_ID
                    parents_parent = None
                    if parent_feature_index != len(parents_hierarchy) - 1:
                        for parents_parent_feature in parents_hierarchy[parent_feature_index + 1:]:
                            if parents_parent_feature in defline_dict:
                                parents_parent = defline_dict[parents_parent_feature]                        
                    if not parent_feature_type in annotation_set.__dict__:
                        annotation_set.__dict__[parent_feature_type] = {}
                    if parent_feature_ID in annotation_set.__dict__[parent_feature_type]:
                        if not child_to_assign in annotation_set.__dict__[parent_feature_type][parent_feature_ID].child_list:
                            annotation_set.__dict__[parent_feature_type][parent_feature_ID].child_list.append(child_to_assign)
                    else:
                        annotation_set.__dict__[parent_feature_type][parent_feature_ID] = ParentAnnotation(parent_feature_ID, seqid,
                                                                                        parent_feature_type, child_list = [child_to_assign],
                                                                                        parent = parents_parent, strand = strand,
                                                                                        annotation_set = annotation_set, other_attributes = {'source':fields[1]})
                    child_to_assign = parent_feature_ID
                #In case of no parent in parent hierarchy:
                if len(parents_hierarchy) == 0:
                    try:
                        annotation_set[parent].child_list.append(ID)
                    except KeyError:
                        print """It seems that this line has a parent attribute but that that parent doesn't have a line itself nor
                        does this line have a defline attribute that specifies a parent type. I'm afraid this function can't currently
                        deal with that."""
                        print ID
                        print parent
                        return None
            #fills other_attributes from defline
            for defline_attribute in defline_dict:
                if not defline_attribute in [IDfield,parent_field]:
                    other_attributes[defline_attribute] = defline_dict[defline_attribute]
            #And now to create the feature!
            if not feature_type in annotation_set.__dict__:
                annotation_set.__dict__[feature_type] = {}
            if feature_type in base_features:
                annotation_set.__dict__[feature_type][ID] = BaseAnnotation( ID, seqid, coords, feature_type, parent,
                                                                           strand , other_attributes, annotation_set)
            else:
                child_list = []
                annotation_set.__dict__[feature_type][ID] = ParentAnnotation(ID, seqid, feature_type, child_list,
                                                                             parent, strand, annotation_set, other_attributes)
    if annotation_set_to_modify == None:
        return copy.deepcopy(annotation_set)


def read_cegma_gff(cegma_gff,annotation_set_to_modify = None):
    """reads gff produced by CEGMA and returns AnnotationSet populated by CEGMA predictions"""
    annotation_set = read_gff(cegma_gff,annotation_set_to_modify = annotation_set_to_modify, presets = "CEGMA")
    if annotation_set_to_modify == None:
        return annotation_set


def read_blast_csv(blast_csv,annotation_set_to_modify = None,hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
    """Reads csv output from blast (-outfmt 10) into an AnnotationSet object. Currently does not string hits together because I'm
    biased towards working on genes in tandem arrays where stringing hits together is annoying. May add option in future."""
    #reads blast_csv from file location, file, or string
    blast_file = ensure_file(blast_csv)
    #checks if annotation_set is given and creates annotation_set if not
    if annotation_set_to_modify == None:
        annotation_set = AnnotationSet()
    else:
        annotation_set = annotation_set_to_modify
    id_generator_dict = {}
    feature_type = hierarchy[-1]
    create_parents_chain = hierarchy[:-1]
    create_parents_chain.reverse()
    if not feature_type in annotation_set.__dict__:
        setattr(annotation_set, feature_type, {})
    if find_truncated_locname:
        if annotation_set.genome == None:
            print '"warning: find_truncated_locname" was set to true, but annotation set has no associated genome object so this cannot be done'
            find_truncated_locname = False
        else:
            genome_seqids = annotation_set.genome.get_seqids()
    for whole_line in blast_file:
        line = whole_line.replace('\r','').replace('\n','')
        fields = line.split(',')
        if len(fields) > 8:
            seqid = fields[1]
            if find_truncated_locname:
                if not seqid in genome_seqids:
                    for genome_seqid in genome_seqids:
                        if seqid == genome_seqid.split()[0]:
                            seqid = genome_seqid
                            break
            tstart = int(fields[8])
            tend = int(fields[9])
            if tstart < tend:
                coords = (tstart,tend)
                strand = '+'
            else:
                coords = (tend,tstart)
                strand = '-'
            score = fields[11]
            IDbase = fields[0]
            if IDbase in eval('annotation_set.' + feature_type):
                ID = IDbase + '-' + str(id_generator_dict[IDbase])
                id_generator_dict[IDbase] = id_generator_dict[IDbase] + 1
                while ID in eval('annotation_set.' + feature_type):
                    ID = IDbase + '-' + str(id_generator_dict[IDbase])
                    id_generator_dict[IDbase] = id_generator_dict[IDbase] + 1
            else:
                ID = IDbase
                id_generator_dict[IDbase] = 1
            other_attributes = {}
            other_attributes['evalue'] = fields[10]
            other_attributes['score'] = score
            parent = ID + '-' + create_parents_chain[0]
            child_to_set = ID
            for parent_index in range(len(create_parents_chain)):
                parent_feature = create_parents_chain[parent_index]
                if not parent_feature in annotation_set.__dict__:
                    setattr(annotation_set, parent_feature, {})
                if parent_index != len(create_parents_chain) - 1:
                    parent_to_set = ID + '-' + create_parents_chain[parent_index + 1]
                else:
                    parent_to_set = None
                annotation_set.__dict__[parent_feature][ID + '-' + parent_feature] = ParentAnnotation(ID + '-' + parent_feature, seqid, parent_feature,
                                                                                                    [child_to_set], parent_to_set, strand,
                                                                                                    annotation_set, other_attributes = {})
                child_to_set = ID + '-' + parent_feature
            eval('annotation_set.' + feature_type)[ID] = BaseAnnotation(ID, seqid, coords, feature_type, parent, strand,
                                                                        other_attributes, annotation_set)

            
    if annotation_set_to_modify == None:
        return annotation_set


class Genotype(tuple):
    """individual genotype to populate a GenotypeDict"""
    def __init__(self, phased = False, phase_group = None, qual = None):
        self.phased = phased
        self.phase_group = phase_group
        self.qual = qual
    

class Variant():
    """individual variants from reference genome sequence. Can be SNPs or polymorphism"""
    def __init__(self, ref_allele, alt_alleles, genotypes = None, variant_set = None, vcf_header_index = None):
        self.variant_set = variant_set
        self.ref_allele = ref_allele
        self.alt_alleles = alt_alleles
        self.genotypes = genotypes
    

class VariantSet(dict):
    """set of variants from reference genome sequence. Can be SNPs or polymorphism"""
    def __init__(self, genome = None, vcf_headers = None):
        self.genome = genome
        self.vcf_headers = vcf_headers
    

class AnnotationSet():
    """A set of annotations of a single genome. Each feature type (e.g. gene, transcript, exon, etc.)
    is stored in it's own dictionary as Annotations with their ID as their key (see "Annotation" class).
    The AnnotationSet itself also functions losely as a dictionary, in that any feature can be returned
    by indexing the AnnotationSet with the ID as a key (e.g. my_annotation_set["my_feature_ID"])"""
    def __init__(self, genome = None):
        self.gene = {}
        self.transcript = {}
        self.CDS = {}
        self.UTR = {}
        self.genome = genome
    
    def __getitem__(self,item):
        all_dicts = {}
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                try:
                    all_dicts[item] = eval("self." + attribute)[item]
                except:
                    pass
        return all_dicts[item]
    
    def read_gff(self, gff, *args, **kwargs):
        kwargs["annotation_set_to_modify"] = self
        read_gff(gff, *args, **kwargs)
    
    def get_seqid(self, seqid):
        seqid_annotation_set = AnnotationSet()        
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                setattr(seqid_annotation_set,attribute,{})
                for feature in eval('self.'+ attribute):
                    feature_obj = eval('self.' + attribute)[feature]
                    if feature_obj.seqid == seqid:
                        eval('seqid_annotation_set.' + feature_obj.feature_type)[feature] = feature_obj
        return seqid_annotation_set
    
    def get_all_seqids(self):
        seqid_list = []
        for attribute in dir(self):
            if type(eval("self." + attribute)) == dict:
                for feature in eval('self.'+ attribute):
                    seqid_list.append(eval('self.' + attribute)[feature].seqid)
        return list(set(seqid_list))
    
    def read_exonerate(self, exonerate_output):
        read_exonerate(exonerate_output,annotation_set_to_modify = self)
    
    def read_blast_csv(self, blast_csv, hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
        read_blast_csv(blast_csv, annotation_set_to_modify = self, hierarchy = hierarchy, source = source, find_truncated_locname = find_truncated_locname)
    
    def read_cegma_gff(self, cegma_gff):
        read_cegma_gff(cegma_gff, annotation_set_to_modify = self)
    
    def get_fasta(self,feature,seq_type = "nucleotide",longest=False,genomic = False):
        fasta_list = []
        for annotation in eval('self.' + feature):
                fasta_list.append(eval('self.' + feature)[annotation].get_fasta(seq_type = seq_type,longest = longest, genomic = genomic))
        return "\n".join(fasta_list)
        
    

class BaseAnnotation():
    """Bottom-most level annotation on a genome, for example CDS, UTR, Match_part, etc. Anything that should have no children"""
    def __init__(self, ID, seqid, coords, feature_type, parent = None, strand = ".", other_attributes = {}, annotation_set = None):
        #Sets up most attributes
        self.ID = ID
        self.seqid = seqid
        self.coords = coords
        self.feature_type = feature_type
        self.annotation_set = annotation_set
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
        self.parent = parent
        self.strand = strand

    def get_coords(self):
        return self.coords
    
    def get_seq(self):
        try:
            if self.strand == '+' or self.strand == '.':
                return Sequence(self.annotation_set.genome.genome_sequence[self.seqid][self.coords[0]-1:self.coords[1]])
            elif self.strand == '-':
                return Sequence(self.annotation_set.genome.genome_sequence[self.seqid][self.coords[0]-1:self.coords[1]]).reverse_compliment()
            else:
                print self.ID + ' has invalid strand value "' + self.strand + '"'
        except:
            print "either base_annotation has not annotation_set, or annotation_set has no genome, or genome has no\
            genome sequence, or genome sequence has no matching seqid, or coords are out of range on that seqid"
            print self.seqid
    
    def get_gff(self, gff_format = "simple gff3"):
        if self.annotation_set != None:
            fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
            fields_list = []
            for field in fields:
                try:
                    fields_list.append(str(eval('self.' + field)))
                except AttributeError:
                    fields_list.append('.')
            if gff_format in ["simple gff3","extended gff3","apollo gff3"]:
                defline = 'ID=' + self.ID
                if self.parent != None:
                    defline = defline + ';Parent=' + self.parent                
            if gff_format == "extended gff3":
                for attribute in self.__dict__:
                    if type(self.__dict__[attribute]).__name__ == 'str' and not attribute in ['ID','Parent','score','strand','seqid','feature_type','phase','source']:
                        defline = defline + ';' + attribute + '=' + self.__dict__[attribute]
            elif gff_format[:13] == "augustus hint":
                gff_format_fields = gff_format.split()
                fields_list[2] = gff_format_fields[2]
                defline = "src=" + gff_format_fields[3]
                if len(gff_format_fields) > 4:
                    defline = defline + ";pri=" + gff_format_fields[4]
            elif gff_format == "gtf":
                defline = 'transcript_id ' + self.parent + ';gene_id ' + self.annotation_set[self.parent].parent
            fields_list.append(defline)
            if gff_format == "apollo gff3" and fields_list[2] == "CDS":
                return '\t'.join(fields_list).replace('\tCDS\t','\texon\t').replace('ID=','ID=ExonOf') + "\n" + '\t'.join(fields_list)
            else:
                return '\t'.join(fields_list)
    


class ParentAnnotation():
    """Parent of any BaseAnnotation. Examples include genes and transcripts. Suggested hierarchy for genes is
    CDS (as BaseAnnotation) -> transcript -> gene."""
    def __init__(self, ID, seqid, feature_type, child_list = [], parent = None, strand = ".", annotation_set = None, other_attributes = {}):
        self.ID = ID
        self.seqid = seqid
        self.feature_type = feature_type
        self.child_list = copy.copy(child_list)
        self.parent = parent
        self.annotation_set = annotation_set
        self.strand = strand
        for attribute in other_attributes:
            setattr(self, attribute, other_attributes[attribute])
    
    def get_coords(self):
        if len(self.child_list) > 0 and self.annotation_set != None:
            coords_list = []
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.__class__.__name__ == 'ParentAnnotation':
                    coords_list = coords_list + list(child_object.get_coords())
                elif child_object.__class__.__name__ == 'BaseAnnotation':
                    coords_list = coords_list + list(child_object.coords)
                else:
                    print "for some reason you have children in ParentAnnotation " + self.ID + " which are neither \
                    ParentAnnotation objects nor BaseAnnotation object. Get your act together"
            return (min(coords_list),max(coords_list))
    
    def get_fasta(self, seq_type = "nucleotide", longest=False, genomic = False, name_from = 'ID'):
        """Returns fasta of this annotation's sequence. If this feature has multiple subfeatures (e.g. this is a gene
        and it has multiple transcripts), the sequence of each subfeature will be an entry in the fasta string."""
        if genomic == True:
            if self.annotation_set.genome != None:
                return ">" + self.ID + '\n' + self.annotation_set.genome.genome_sequence[self.seqid][self.get_coords()[0] - 1:self.get_coords()[1]] + '\n'
        elif len(self.child_list) > 0 and self.annotation_set != None:
            if self.annotation_set.genome != None:
                fasta_list = []
                child_type = self.annotation_set[self.child_list[0]].__class__.__name__
                if child_type == 'BaseAnnotation':
                    seq_list = []
                    child_dict = {}
                    for child in self.child_list:
                        child_obj = self.annotation_set[child]
                        try:
                            child_dict[child_obj.coords] = child_obj.get_seq()
                        except AttributeError:
                            print "ParentAnnotation has both ParentAnnotation and BaseAnnotation children!"
                            print self.ID
                        strand = child_obj.strand
                    children_in_correct_order = list(child_dict)
                    children_in_correct_order.sort()
                    if strand == '-':
                        children_in_correct_order.reverse()
                    for child in children_in_correct_order:
                        seq_list.append(child_dict[child])
                    if seq_type == "nucleotide":
                        new_seq = Sequence("".join(seq_list))
                    elif seq_type == "protein":
                        new_seq = Sequence("".join(seq_list)).translate()
                    else:
                        print seq_type + ' is not valid seq_type. Please specify "protein" or "nucleotide".'
                    fasta_list.append('>' + self.__dict__[name_from] + '\n' + new_seq)
                else:
                    for child in self.child_list:
                        try:
                            child_fasta = self.annotation_set[child].get_fasta(seq_type=seq_type, name_from = name_from)
                            if child_fasta != "":
                                fasta_list.append(child_fasta)
                        except AttributeError:
                            print "ParentAnnotation has both ParentAnnotation and BaseAnnotation children!"
                            print self.ID
                if longest == True:
                    seqlens = {}
                    for seq in fasta_list:
                        seqlens[len("".join(seq.split('\n')[1:]))] = seq
                    return seqlens[max(list(seqlens))]
                else:
                    try:
                        return '\n'.join(fasta_list)
                    except:
                        return ""
            else: return ""
        else: return ""
    
    def get_gff(self, gff_format = "simple gff3"):
        """presets currently include: "simple gff3", "extended gff3", "gtf", "apollo gff3", "maker hint", and "augustus hint". Will eventually include more by request"""
        if self.annotation_set != None:
            fields = ['seqid','source','feature_type','get_coords()[0]','get_coords()[1]','score','strand','phase']
            fields_list = []
            parent_line = True
            for field in fields:
                try:
                    fields_list.append(str(eval('self.' + field)))
                except AttributeError:
                    fields_list.append('.')
            if gff_format in ["simple gff3","extended gff3","apollo gff3"]:
                defline = 'ID=' + self.ID
                if self.parent != None:
                    defline = defline + ';Parent=' + self.parent
            if gff_format == "extended gff3":
                for attribute in self.__dict__:
                    if type(self.__dict__[attribute]).__name__ == 'str' and not attribute in ['ID','Parent','score','strand','seqid','feature_type','phase','source']:
                        defline = defline + ';' + attribute + '=' + self.__dict__[attribute]
            elif gff_format == "apollo gff3":
                if 'Name' in self.__dict__:
                    defline = defline + ';Name=' + self.Name
                elif 'name' in self.__dict__:
                    defline = defline + ';Name=' + self.name
                else:
                    defline = defline + ';Name=' + self.ID
            elif gff_format[:13] == "augustus hint":
                parent_line = False
            elif gff_format == 'gtf':
                parent_line = False
            if parent_line:
                fields_list.append(defline)
                lines_list = ['\t'.join(fields_list)]
            else:
                lines_list = []
            child_dict = {}
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.get_coords() not in child_dict:
                    child_dict[child_object.get_coords()] = child_object.get_gff(gff_format)
                else:
                    child_dict[(child_object.get_coords()[0],child_object.get_coords()[1] + len(child_dict))] = child_object.get_gff(gff_format)
            child_coords = list(child_dict)
            child_coords.sort()
            for child_index in child_coords:
                lines_list.append(child_dict[child_index])
            if gff_format == "apollo gff3":
                lines_list = '\n'.join(lines_list).split('\n')
                for line in lines_list[:]:
                    if line.split('\t')[2] == "CDS":
                        lines_list.remove(line)
                        lines_list.append(line)
            return '\n'.join(lines_list)
    
    def get_bed(self, columns = 12):
        """If children are base annotations, returns a bed-format line with the number of columns specified. Otherwise,
        returns bed-format lines for each child"""
        if self.annotation_set != None:
            lines_list = []
            child_blocks = []
            for child in self.child_list:
                child_object = self.annotation_set[child]
                if child_object.__class__.__name__ == "ParentAnnotation":
                    lines_list.append(child_object.get_bed(columns = columns))
                elif child_object.__class__.__name__ == "BaseAnnotation":
                    child_blocks.append((child_object.get_coords()[0] - self.get_coords()[0],
                                         child_object.get_coords()[1] - child_object.get_coords()[0] + 1))
            if child_blocks != []:
                if lines_list != []:
                    print "Weirdo feature " + self.ID + ", it has both ParentAnnotation and BaseAnnotation children.\
                    I'll output the line for this feature and it's ParentAnnotation children, but\
                    something is probably wrong"
                child_blocks.sort()
                child_block_sizes = []
                child_block_starts = []
                for child_block in child_blocks:
                    child_block_starts.append(str(child_block[0]))
                    child_block_sizes.append(str(child_block[1]))
                potential_fields = ['seqid','get_coords()[0] - 1','get_coords()[1]','ID','score','strand','thickStart','thickEnd','itemRgb']
                fields = []
                for field_number in range(columns):
                    if field_number < 9:
                        try:
                            fields.append(str(eval("self." + potential_fields[field_number])))
                        except AttributeError:
                            fields.append('.')
                    elif field_number > 11:
                        print "Only 12 fields in bed format!"
                    elif field_number == 9:
                        fields.append(str(len(child_block_sizes)))
                    elif field_number == 10:
                        fields.append(','.join(child_block_sizes))
                    elif field_number == 11:
                        fields.append(','.join(child_block_starts))
                lines_list.append('\t'.join(fields))
            return "\n".join(lines_list)
    


class Sequence(str):
    """DNA sequence. Has methods allowing reverse complimenting,
        translating, etc."""
    def reverse_compliment(self):
        """returns reverse compliment of self"""
        new_sequence_list = []
        compliment_dict = {'a':'t','t':'a','g':'c','c':'g','A':'T','T':'A','G':'C','C':'G','n':'n','N':'N','-':'-'}
        for residue in self[::-1]:
            try:
                new_sequence_list.append(compliment_dict[residue])
            except KeyError:
                new_sequence_list.append('n')
        return Sequence(''.join(new_sequence_list))
    
    def translate(self,library = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                                  'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                                  'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                                  'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                                  'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                                  'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                                  'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                                  'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'},
                  frame = 0, strand = '+',trimX = True):
        triplet = ""
        newseq = ""
        if strand == '+':
            seq = self
        elif strand == '-':
            seq = self.reverse_compliment()
        if len(seq) > (2 + frame):
            for residue_position in range(frame, len(self)):
                triplet = triplet + seq[residue_position].upper()
                if (residue_position + frame) % 3 == 2:
                    try:
                        newseq = newseq + library[triplet]
                    except KeyError:
                        newseq = newseq + 'X'
                    triplet = ""
            if trimX:
                if newseq[0] == 'X':
                    newseq = newseq[1:]
            return newseq
    
    def get_orfs(self, longest = False, strand = 'both', from_atg = False):
        if len(self) > 3:
            orflist = []
            if longest:
                candidate_list = []
                longest_orf_len = 0
            for frame in [0,1,2]:
                for strand in ['-','+']:
                    translated_seq = self.translate(frame=frame,strand=strand)
                    if translated_seq:
                        translated_seq_list = translated_seq.split('*')
                        for orf in translated_seq_list:
                            if from_atg:
                                try:
                                    output_orf = 'M' + ''.join(orf.split('M')[1:])
                                except IndexError:
                                    continue
                            else:
                                output_orf = orf
                            if longest:
                                if len(output_orf) > longest_orf_len:
                                    candidate_list.append(output_orf)
                                    longest_orf_len = len(output_orf)
                            else:
                                orflist.append(output_orf)
            if longest:
                return candidate_list[-1]
            else:
                return orflist


class GenomeSequence(dict):
    """genome sequence class, currently takes input in multi-fasta format."""
    def __init__(self,genome_sequence = None, truncate_names = False):
        #reads input file location, file, or string
        sequence_file = ensure_file(genome_sequence)
        #breaks sequence into sequence blocks (contigs, scaffolds, or chromosomes), adds sequence 
        #   from each block as dictionary entry into self with block name as key.
        if sequence_file != None:
            seq = ""
            seqname = ""
            for line in sequence_file:
                if line[0] == ">":
                    if truncate_names == True:
                        seqid = line[1:].replace('\r','').replace('\n','').split()[0]
                    else:
                        seqid = line[1:].replace('\r','').replace('\n','')
                    if seq != "":
                        self[seqname] = seq
                        seq = ""
                    seqname = seqid
                else:
                    seq = seq + line.replace('\r','').replace('\n','')
            if seq != "":
                self[seqname] = seq


class Genome():
    """genome class, which contains sequence and annotations. Annotations can be given as annotation_set object, gff3, cegma_gff,
    blast_csv, or exonerate_output (just set annotation_format)."""
    def __init__(self,genome_sequence = None, annotations = None, varients = None, annotation_format = 'annotation_set', truncate_names = False):
        if genome_sequence.__class__.__name__ == 'GenomeSequence' or genome_sequence == None:
            self.genome_sequence = genome_sequence
        else:
            self.genome_sequence = GenomeSequence(genome_sequence, truncate_names = truncate_names)
        if annotations != None:
            if annotations.__class__.__name__ == "AnotationSet" and annotation_format == 'annotation_set':
                self.annotations = annotations
                self.annotations.genome = self
            elif annotation_format == 'gff3':
                self.annotations = read_gff(annotations)
                self.annotations.genome = self
            elif annotation_format == 'cegma_gff':
                self.annotations = read_cegma_gff(annotations)
                self.annotations.genome = self
            elif annotation_format == 'blast_csv':
                self.annotations = read_blast_csv(annotations)
                self.annotations.genome = self
            elif annotation_format == 'exonerate_output':
                self.annotations = read_exonerate(annotations)
                self.annotations.genome = self
        else:
            self.annotations = annotations
    
    def get_scaffold_fasta(self, seqid):
        return '>' + seqid + '\n' + self.genome_sequence[seqid]
    
    def get_genome_fasta(self, remove_spaces = False):
        fasta_list = []
        for seqid in self.genome_sequence:
            if remove_spaces:
                fasta_header = seqid.split()[0]
            else:
                fasta_header = seqid
            fasta_list.append('>' + fasta_header + '\n' + self.genome_sequence[seqid])
        return "\n".join(fasta_list)
    
    def write_apollo_gff(self, seqid, suppress_fasta = False):
        if self.genome_sequence != None and self.annotations != None:
            try:
                apollo_gff = write_gff(self.annotations.get_seqid(seqid),gff_format = 'apollo gff3')
                if not suppress_fasta:
                    apollo_gff = apollo_gff + '\n' + self.get_scaffold_fasta(seqid)
                return apollo_gff
            except:
                if not suppress_fasta:
                    return self.get_scaffold_fasta(seqid)
                else:
                    return ""
        else:
            print "genome object is either missing genome_sequence or annotations"
    
    def get_seqids(self, from_annotations = False):
        seqid_list = []
        warning = False
        if self.genome_sequence != None:
            for seqid in self.genome_sequence:
                seqid_list.append(seqid)
        if self.annotations != None and from_annotations:
            for seqid in self.annotations.get_all_seqids():
                if seqid not in seqid_list:
                    seqid_list.append(seqid)
                    warning = True
        if warning:
            print "warning, some annotations possessed seqids not found in sequence dictionary"
        return seqid_list
    
    def read_exonerate(self, exonerate_output):
        if self.annotations != None:
            self.annotations.read_exonerate(exonerate_output)
        else:
            self.annotations = read_exonerate(exonerate_output)
            self.annotations.genome = self
    
    def read_blast_csv(self, blast_csv, hierarchy = ['match','match_part'], source = 'blast', find_truncated_locname = False):
        if self.annotations == None:
            self.annotations = AnnotationSet()
            self.annotations.genome = self
        self.annotations.read_blast_csv(blast_csv, hierarchy = hierarchy, source = source, find_truncated_locname = find_truncated_locname)
    
    def read_cegma_gff(self, cegma_gff):
        if self.annotations != None:
            self.annotations.read_cegma_gff(cegma_gff)
        else:
            self.annotations = read_cegma_gff(cegma_gff)
            self.annotations.genome = self

    def read_gff(self, gff, *args, **kwargs):
        if self.annotations != None:
            self.annotations.read_gff(gff, *args, **kwargs)
        else:
            self.annotations = read_gff(gff, *args, **kwargs)
            self.annotations.genome = self
    
    def read_vcf(self, vcf):
        self.variants = read_vcf(vcf)
    

class position_dic(dict):
    def __init__(self, genome_sequence, dtype=bool):
        for seqid in genome_sequence:
            self[seqid] = numpy.zeros(len(genome_sequence[seqid]),dtype=dtype)

    def fill_from_annotations(self, annotation_set, feature, fill_type = "coords",fill_with = "1"):
        """accepts "start" and "coords" for fill_type"""
        for annotation in eval("annotation_set." + feature):
            annotation_obj = eval("annotation_set." + feature)[annotation]
            seqid = annotation_obj.seqid
            coords = annotation_obj.get_coords()
            try:
                parent = annotation_obj.parent
            except:
                pass
            ID = annotation_obj.ID
            if fill_type == "coords":
                for position in range(coords[0] - 1,coords[1]):
                    self[seqid][position] = eval(fill_with)
            elif fill_type == "start":
                self[seqid][coords[0] - 1] = eval(fill_with)
    
    def count_from_annotations(self, annotation_set, feature):
        count_list = []
        for annotation in eval("annotation_set." + feature):
            annotation_obj = eval("annotation_set." + feature)[annotation]
            seqid = annotation_obj.seqid
            coords = annotation_obj.get_coords()
            featureID = annotation_obj.ID
            try:
                parent = annotation_obj.parent
            except:
                pass
            ID = annotation_obj.ID
            count0 = 0
            count1 = 0
            for position in range(coords[0] - 1,coords[1]):
                if self[seqid][position] == 0:
                    count0 = count0 + 1
                elif self[seqid][position] == 1:
                    count1 = count1 + 1
            count_list.append([featureID,count0,count1])
        return count_list
    

    def at_content(self, genome_sequence):
        for seqid in self:
            for position in range(len(self[seqid])):
                if genome_sequence[seqid][position] in "ATat":
                    self[seqid][position] = 1

    def sliding_window_calculate(self, window_size, window_jump = 1, operation = "sum", output = "dict",
                                 threshold = 1, seqs_to_exclude = []):
        """operation may be "sum", "average", or "set". Output may be "dict","coords", or "annotation_set"."""
        if output == "dict":
            new_dic = {}
        elif output == "annotation_set":
            annotation_set = AnnotationSet()
            annotation_set.region = {}
        elif output == "coords":
            coords_list = []
        for seqid in self:
            if len(self[seqid]) > window_size and not seqid in seqs_to_exclude:
                if output == "dict":
                    new_dic[seqid] = []
                window_in = False
                coords_list = []
                for position in range(len(self[seqid][:-1 * window_size]) / window_jump):
                    window_start = position * window_jump
                    if operation == "set":
                        new_dic[seqid].append(set(self[seqid][window_start:window_start + window_size]))
                        continue
                    window_sum = numpy.sum(self[seqid][window_start:window_start + window_size])
                    if operation == "sum":
                        value = window_sum
                    elif operation == "average":
                        value = window_sum * 1.0 / window_size
                    if output == "dict":
                        new_dic[seqid].append(value)
                    elif output == 'annotation_set':
                        if value >= threshold:
                            if window_in == False:
                                window_in = True
                                if len(coords_list) == 0:
                                    coords_list.append([1 + window_start])
                                elif position > coords_list[-1][1]:
                                    coords_list.append([1 + window_start])
                        else:
                            if window_in == True:
                                window_in = False
                                if len(coords_list[-1]) == 1:
                                    coords_list[-1].append(window_start + window_size)
                                else:
                                    coords_list[-1][1] = window_start + window_size
                    elif output == "coords":
                        if threshold[0] <= value <= thredshold[1]:
                            coords_list.append([seqid, 1 + window_start, window_start + window_size])
                    elif type(output) == file:
                        output.write(seqid + '\t' + str(value) + '\n')
                    if window_start % 10000000 == 0 and verbose:
                        print "processed " + seqid + " to position " + str(position)
                if output == "annotation_set":
                    if len(coords_list) > 0:
                        if len(coords_list[-1]) == 1:
                            coords_list[-1].append(len(self[seqid]))
                        for coords in coords_list:
                            ID = seqid + "-window" + str(coords[0])
                            annotation_set.region[ID] = BaseAnnotation(ID, seqid, tuple(coords), "region", annotation_set = annotation_set)
                if verbose:
                    print "processed " + seqid
                    if output == 'annotation_set':
                        for region in annotation_set.region:
                            if annotation_set.region[region].seqid == seqid:
                                print str(annotation_set.region[region].coords)                    
        if output == "dict":
            return new_dic
        elif output == "annotation_set":
            return copy.deepcopy(annotation_set)
        elif output == "coords":
            return coords_list





