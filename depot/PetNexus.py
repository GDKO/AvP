from Bio import SeqIO
from ete3 import Tree
from collections import Counter

result_nohgt = "no_hgt"
result_unknown = "unknown_topology"
result_complex = "complex_topology"
result_hgt = "hgt"

def add_color(str_id, colors):

    # Color scheme in colors.txt
    if ("TOI" in str_id):
        return str_id + "[&!color=#" + colors['TOI'] + "]"
    elif("StudiedOrganism" in str_id):
        return str_id + "[&!color=#" + colors['StudiedOrganism'] + "]"
    elif("EGP" in str_id):
        return str_id + "[&!color=#" + colors['EGP'] + "]"

    for name in colors:
        if name in str_id:
            return str_id + "[&!color=#" + colors[name] + "]"

    return str_id + "[&!color=#0048B3]"

def make_nexus_file(gene, group, lineage, gene_nexus_path, group_file, phylogeny_file, colors):

    f_nexus = open(gene_nexus_path,'w')

    #BLOC TAXA
    ntax = 0
    nchar = 0
    tax = []

    for seq_record in SeqIO.parse(group_file, "fasta"):
        geneid = seq_record.id.split('@')[0]
        geneid.replace(">","")

        if geneid in lineage:
            tax.append(add_color(seq_record.id+str(lineage[geneid]),colors))
        elif geneid == gene:
            tax.append(seq_record.id+"[&!color=#000000]")
        else:
            tax.append(add_color(seq_record.id,colors))

        nchar = len(seq_record.seq)
        ntax+= 1

    f_nexus.write("#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS ntax="+str(ntax)+";")
    f_nexus.write("\n\tTAXLABELS\n"+ '\n'.join(tax)+";" +"\nEND;\n\n")

    #BLOC CHARACTERS
    f_nexus.write("BEGIN CHARACTERS;\n\tDIMENSIONS nchar="+str(nchar)+';\n\tFORMAT datatype=Protein gap=-;\n\tMATRIX\n')
    for seq_record in SeqIO.parse(group_file, "fasta"):
        f_nexus.write("\t"+seq_record.id+"\n"+str(seq_record.seq)+"\n")
    f_nexus.write(";\nEND;\n\n")

    #BLOC TREES

    t = Tree(phylogeny_file)

    R = t.get_midpoint_outgroup()
    if(R != None):
        t.set_outgroup(R)

    f_nexus.write("BEGIN TREES;\n\tTREE tree_1="+t.write()+"\n")

    f_nexus.write("END;")
    f_nexus.close()


def analyze_tree(tree_filename, full_name_studied_gene, node_support, complex_per):
    global result_nohgt
    global result_hgt
    global result_complex
    global result_unknown

    # Load a tree structure from a newick file
    gene_tree = Tree(tree_filename,format=0)


    if node_support != 0:
        node_supports = []
        for node in gene_tree.traverse("preorder"):
            node_supports.append(node.support)
        if all(i <= 1 for i in node_supports):
            node_support = node_support/100
        for node in gene_tree.traverse("preorder"):
            if "@" not in node.name:
                if node.support < node_support:
                    node.delete()

    no_TOI = True
    only_TOI = True

    # Check if no_TOI or only_TOI to speed up calculations
    for node in gene_tree:
        if "@TOI" in str(node):
            no_TOI = False
        elif "EGP" not in str(node) or "StudiedOrganism" not in str(node):
            only_TOI = False


    if only_TOI:
        return "only_TOI"

    if no_TOI:
        return "no_TOI"

    # Root the tree using the midpoint
    R = gene_tree.get_midpoint_outgroup()
    if(R != None):
        gene_tree.set_outgroup(R)

    return analysis(gene_tree, full_name_studied_gene, complex_per)

def analysis(gene_tree, full_name_studied_gene, complex_per):

    global result_nohgt
    global result_hgt
    global result_complex
    global result_unknown

    classify_nodes(gene_tree,full_name_studied_gene)
    node2labels = gene_tree.get_cached_content(store_attr="name")
    def collapsed_leaf(node):
        if len(node2labels[node]) == 1:
            return True
        else:
            return False

    list_of_ranks = ["H", "T"]

    sister_nodes_list = get_sister_classification(gene_tree,full_name_studied_gene,list_of_ranks)

    sister_clade_full = sister_nodes_list[0]
    ancestral_clade_full = sister_nodes_list[1]

    sister_clade = sister_clade_full.keys()
    sister_per = (sister_clade_full["H"]/(sister_clade_full["H"]+sister_clade_full["T"]))*100
    if sister_per < complex_per[0]:
        sister_clade = ["T"]
    elif sister_per > (complex_per[1]):
        sister_clade = ["H"]

    ancestral_clade = ancestral_clade_full.keys()
    if ancestral_clade:
        ancest_per = (ancestral_clade_full["H"]/(ancestral_clade_full["H"]+ancestral_clade_full["T"]))*100
        if ancest_per < complex_per[0]:
            ancestral_clade = ["T"]
        elif ancest_per > (complex_per[1]):
            ancestral_clade = ["H"]

    if "H" in sister_clade and "T" in sister_clade:
        if "T" in ancestral_clade and "H" not in ancestral_clade:
            return result_nohgt
        else:
            return result_complex
    elif "T" in sister_clade:
        return result_nohgt
    elif "H" in sister_clade:
        if "T" in ancestral_clade:
            return result_complex
        else:
            return result_hgt

    return result_unknown



def classify_nodes(gene_tree,full_name_studied_gene):

    # Label leaves
    for node in gene_tree.traverse("preorder"):
        if("@" in node.name):
            if "@TOI" in node.name:
                node.name = "@T"
            elif "@EGP" in node.name:
                node.name = "@E"
            elif("Unknown" in node.name or "Other" in node.name or "Unclassified" in node.name):
                node.name = "@U"
            elif full_name_studied_gene in node.name:
                pass
            elif "@StudiedOrganism" in node.name:
                node.name = "@E"
            else:
                node.name = "@H"

    # Label internal nodes
    for node in gene_tree.traverse("preorder"):
        if("@" not in node.name):
            if(len(node.children)>0):
                children = node.get_leaf_names()
                if len(set(children)) == 1:
                    node.name = set(children).pop()
                else:
                    node.name = "@L"

def classify_tree(tree_filename,full_name_studied_gene,ranks,complex_per_node):

    gene_tree = Tree(tree_filename,format=0)
    R = gene_tree.get_midpoint_outgroup()
    if(R != None):
        gene_tree.set_outgroup(R)

    sister_nodes_list = get_sister_classification(gene_tree,full_name_studied_gene,ranks)

    final_nodes = sister_nodes_list[0] + sister_nodes_list[1]

    per_of_most_common_rank = (final_nodes.most_common(1)[0][1]/sum(final_nodes.values()))*100

    if per_of_most_common_rank > complex_per_node:
        final_nodes = {final_nodes.most_common(1)[0][0]:final_nodes.most_common(1)[0][1]}

    return list(final_nodes.keys())


def get_sister_classification(gene_tree,full_name_studied_gene,ranks):
    sister_nodes = {}
    sister_nodes["sister"]= []
    sister_nodes["ancestral"]= []
    gene_node = gene_tree.get_leaves_by_name(full_name_studied_gene)[0]
    ancestral_nodes = gene_node.get_ancestors()
    sisters_node_flag = True

    sister_children = []
    for ancestral_node in ancestral_nodes:
        if sisters_node_flag:
            children_old = ancestral_node.get_leaf_names()
            children_new = []
            for child in children_old:
                tc = child.split("@")
                if tc[-1] in ranks:
                    children_new.append(tc[-1])
                    sister_children = children_new
                    sisters_node_flag = False
            sister_nodes["sister"] = sister_children
        else:
            children_old = ancestral_node.get_leaf_names()
            children_new = []
            for child in children_old:
                tc = child.split("@")
                if tc[-1] in ranks:
                    children_new.append(tc[-1])
            for item in sister_children:
                children_new.remove(item)
            sister_nodes["ancestral"]= children_new
            if len(children_new) > 0:
                break
    sister_nodes_list = []
    sister_nodes["sister"]=Counter(sister_nodes["sister"])
    sister_nodes["ancestral"]=Counter(sister_nodes["ancestral"])
    sister_nodes_list = [sister_nodes["sister"], sister_nodes["ancestral"]]

    return sister_nodes_list
