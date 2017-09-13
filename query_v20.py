import os, time, glob
from data_v11 import geneidchoice   #, genenamechoice

from model_v12 import InputForm_gui_v12
form = InputForm_gui_v12()

gene_id_t = {}
gene_chr_t = {}
for item in geneidchoice:
    gene_id_t[item[1]] = str(item[0]);
    gene_chr_t[item[1]] = item[2]

def gene2chr(id):
    try:
        return gene_chr_t[id];
    except KeyError:
        return 'MT'; # return MT for unfound id.

def variant2chr(id):
    if id>=1 and id<=147569: return str(10);
    elif id>=147570 and id<=399060: return str(11);
    elif id>=399061 and id<=504981: return str(12);
    elif id>=504982 and id<=623452: return str(13);
    elif id>=623453 and id<=725584: return str(14);
    elif id>=725585 and id<=841302: return str(15);
    elif id>=841303 and id<=942575: return str(16);
    elif id>=942576 and id<=1083978: return str(17);
    elif id>=1083979 and id<=1165513: return str(18);
    elif id>=1165514 and id<=1259101: return str(19);
    elif id>=1259102 and id<=1467268: return str(1);
    elif id>=1467269 and id<=1749663: return str(2);
    elif id>=1749664 and id<=1904274: return str(3);
    elif id>=1904275 and id<=2136473: return str(4);
    elif id>=2136474 and id<=2346838: return str(5);
    elif id>=2346839 and id<=2523882: return str(6);
    elif id>=2523883 and id<=2753944: return str(7);
    elif id>=2753945 and id<=2903006: return str(8);
    elif id>=2903007 and id<=3072372: return str(9);
    elif id>=3072373 and id<=3175693: return str(X);
    elif id>=3175694 and id<=3176668: return str(Y);
    else: return str(MT); # also return MT for unfound id.

def variant_list2chr(variant_id_list):
    variant_list = variant_id_list.replace(" ", "").split(',')
    variant_list = [ int(s) for s in variant_list]
    variant_list.sort()
    # Find chr file need for variant
    chr_list = [ variant2chr(s) for s in variant_list ]
    return zip(chr_list,variant_list)

def gene_list2chr(gene_id_list):
    gene_list = gene_id_list.replace(" ", "").split(',')
    chr_list=[]
    for s in gene_list:
        try:
            chr_list.append(gene_chr_t[s])
        except:
            chr_list.append('MT')
#chr_list = [ gene_chr_t[s] for s in gene_list ]
    return zip(chr_list,gene_list)

def generatequery_gui_sub(form):
    
    # Input: a form object from corresponding model
    variantinputtype = form.variantinputtype.data
    view = form.view.data
    variant_id =form.variant_id_text.data
    chr = form.chr.data
    strnum = int(form.strnum.data)
    endnum = int(form.endnum.data)
    genesearch=form.genesearch.data
    selectCollum1=form.selectCollum1.data
    selectCollum2=form.selectCollum2.data
    selectCollum3=form.selectCollum3.data
    selectCollum4= form.selectCollum4.data
    prob_cutoff= form.prob_cutoff.data
    ismax= form.ismax.data
    homo = form.homo.data
    consequence_list=form.consequence.data
    strainlistnew = form.strainlistnew.data
    strainp1new=form.strainp1new.data
    strainp2new=form.strainp2new.data
    limnum=form.limnum.data
    
    FROM = ""
    SELECT = []
    
    consequence_list_b = []
    for x in consequence_list:
        for y in x.split(','):
            consequence_list_b.append(y)
    consequence_list = list(set(consequence_list_b))
    
    if variantinputtype=='1':
        #if variantfile != 1:
        variant_list = variant_id.replace(" ", "").split(',')
    else:
        pass
    if variantinputtype=='2':
        #if bedfile == None:
        [chr,strnum,endnum]=[chr,strnum,endnum]
    if variantinputtype == '3':
        #if genefile != 1:
        gene_list=[]
        for s in genesearch.replace(" ", "").split(','):
            try:
                gene_list.append(gene_id_t[s])
            except KeyError:
                gene_list.append(0)

    if view == '1':
        WHERE_list=[]
        selectCollum = selectCollum1
        FROM = 'exon_1410E_' + str(chr) + '.genotype_transcript_view'

        if 'consequence' in selectCollum:
            index = selectCollum.index('consequence')
            selectCollum = selectCollum[0:index]+ ["consequence_1","consequence_2"]+ selectCollum[index+1:]
        if 'allele' in selectCollum:
            index = selectCollum.index('allele')
            selectCollum = selectCollum[0:index]+ ["allele_1","allele_2"]+ selectCollum[index+1:]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["allele_id_1 <=> allele_id_2"]
        if homo == 'Hetero':
            WHERE_list += ["NOT allele_id_1 <=> allele_id_2"]
        if len(consequence_list) > 0 :
        # WHERE_list += ["( " + " OR ".join(["consequence_id_1=" + str(s) + " or consequence_id_2=" + str(s) for s in consequence_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["consequence_1 LIKE '%" + str(s) + "%' or consequence_2 LIKE '%" + str(s) + "%'" for s in consequence_list]) + " )"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["prob>" + str(prob_cutoff)]
        if ismax:
            WHERE_list += ["is_max=1"]
    
        if len(strainlistnew)>0:
            WHERE_list += ["( " + " OR ".join(["strain_id = " + str(s) for s in strainlistnew]) + " )"]

        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT distinct "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT distinct "+ SELECT + " FROM " + FROM

    ###################################################
    if view == '2':
        WHERE_list=[]
        selectCollum = selectCollum2
        FROM = 'exon_1410E_' + str(chr) + '.diplotype_view'

        if 'founder_name' in selectCollum:
            index = selectCollum.index('founder_name')
            selectCollum = selectCollum[0:index]+ ["founder_name_1","founder_name_2"]+ selectCollum[index+1:]
        
        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["founder_name_1 <=> founder_name_2"]
        if homo == 'Hetero':
            WHERE_list += ["NOT founder_name_1 <=> founder_name_2"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["prob>" + str(prob_cutoff)]

        if len(strainlistnew)>0:
            WHERE_list += ["( " + " OR ".join(["strain_id = " + str(s) for s in strainlistnew]) + " )"]

        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT distinct "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT distinct "+ SELECT + " FROM " + FROM


    ###################################################
    if view == '3':
        WHERE_list=[]
        selectCollum = selectCollum3
        FROM = 'exon_1410E_' + str(chr) + ".genotype_sampling_view as s1 inner join " + 'exon_1410E_' + str(chr) + ".genotype_sampling_view as s2 on s1.variant_id=s2.variant_id and s1.transcript_id=s2.transcript_id";
        
        strain_pair = [(strainp1new,strainp2new)]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["s1.allele_index <=> s2.allele_index"]
        if homo == 'Hetero':
            WHERE_list += ["NOT s1.allele_index <=> s2.allele_index"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["s1.variant_id = " + str(s) for s in variant_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["s1.chrom='"+str(chr)+"'"]
                WHERE_list += ["s2.chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["s1.pos>="+ str(strnum)]
                WHERE_list += ["s2.pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["s1.pos<="+ str(endnum)]
                WHERE_list += ["s2.pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["s1.gene_id = " + str(s) for s in gene_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["s1.prob*s2.prob >" + str(prob_cutoff)]
        
        WHERE_list_strain = []
        for s2 in strain_pair:
            WHERE_list_strain +=  ["(s1.strain_id=" + str(s2[0]) + " AND s2.strain_id=" + str(s2[1]) + " )"]
        if len(consequence_list) > 0 :
            WHERE_list += ["( " + " OR ".join(["s1.consequence_name like '%" + str(s) + "%' OR s2.consequence_name like '%" + str(s) + "%'" for s in consequence_list]) + " )"]

        if ismax:
            WHERE_list += [" s1.is_max=1 AND s2.is_max=1 "]

        WHERE_list += ["( " + " OR ".join([ s for s in WHERE_list_strain]) + " )"]


        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT "+ SELECT + " FROM " + FROM

    ###################################################
    if view == '4':
        WHERE_list=[]
        selectCollum = selectCollum4
        FROM = 'exon_1410E_' + str(chr) + ".diplotype_sampling_view as s1 inner join "+ 'exon_1410E_' + str(chr) + ".diplotype_sampling_view as s2 on s1.variant_id=s2.variant_id";
        
        strain_pair = [(strainp1new,strainp2new)]

        SELECT = ",".join(selectCollum)

        if homo == 'Homo':
            WHERE_list += ["s1.founder_id <=> s2.founder_id"]
        if homo == 'Hetero':
            WHERE_list += ["NOT s1.founder_id <=> s2.founder_id"]

        if variantinputtype=='1':
            WHERE_list += ["( " + " OR ".join(["s1.variant_id = " + str(s) for s in variant_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.variant_id = " + str(s) for s in variant_list]) + " )"]
        if variantinputtype=='2':
            if chr != "ALL":
                WHERE_list += ["s1.chrom='"+str(chr)+"'"]
                WHERE_list += ["s2.chrom='"+str(chr)+"'"]
            if strnum != None:
                WHERE_list += ["s1.pos>="+ str(strnum)]
                WHERE_list += ["s2.pos>="+ str(strnum)]
            if endnum != None:
                WHERE_list += ["s1.pos<="+ str(endnum)]
                WHERE_list += ["s2.pos<="+ str(endnum)]
        if variantinputtype == '3':
            WHERE_list += ["( " + " OR ".join(["s1.gene_id = " + str(s) for s in gene_list]) + " )"]
            WHERE_list += ["( " + " OR ".join(["s2.gene_id = " + str(s) for s in gene_list]) + " )"]

        if prob_cutoff>0:
            WHERE_list += ["s1.prob*s2.prob >" + str(prob_cutoff)]
        
        WHERE_list_strain = []
        for s2 in strain_pair:
            WHERE_list_strain +=  ["(s1.strain_id=" + str(s2[0]) + " AND s2.strain_id=" + str(s2[1]) + " )"]

        WHERE_list += ["( " + " OR ".join([ s for s in WHERE_list_strain]) + " )"]


        if len(WHERE_list)>0:
            WHERE = " and ".join(WHERE_list)
            query = "SELECT "+ SELECT + " FROM " + FROM + " WHERE " + WHERE
        else:
            query = "SELECT "+ SELECT + " FROM " + FROM

    if (int(limnum)==0):
        return str(query)
    else:
        return str(query) + " limit " + str(int(limnum))


def generatequery_gui_v20(form):
    # Record the original value of form
    variantinputtype = form.variantinputtype.data
    view = form.view.data
    variant_id =form.variant_id_text.data
    chr = form.chr.data
    genesearch=form.genesearch.data
    limnum=form.limnum.data
    
    form.limnum.data='0' # important, the limnum should be added at the end of Union
    
    if variantinputtype=='2' and chr=='ALL':
        form.chr.data=1
        OUT = generatequery_gui_sub(form);
        OUT = str(" (") + OUT + str(") ");
        for i in range(1,20)+['X','Y','MT']:
            form.chr.data=i
            OUT = OUT + " UNION " + " (" + generatequery_gui_sub(form)+ ") ";
    elif variantinputtype=='1':
        input_in= variant_list2chr(variant_id)[0]
        form.variant_id_text.data= str(input_in[1])
        form.chr.data= str(input_in[0])
        
        OUT = generatequery_gui_sub(form);
        
        if len(variant_list2chr(variant_id))>1:
            OUT = str(" (") + OUT + str(") ");
            for input_in in variant_list2chr(variant_id)[1:]:
                form.variant_id_text.data= str(input_in[1])
                form.chr.data= input_in[0]
                OUT = OUT + " UNION " + " (" + generatequery_gui_sub(form)+ ") ";
    elif variantinputtype=='3':
        # Search for the proper chr to use
        input_in= gene_list2chr(genesearch)[0]
        form.genesearch.data= str(input_in[1])
        form.chr.data= str(input_in[0])
        
        OUT = generatequery_gui_sub(form);
        if len(gene_list2chr(genesearch))>1:
            OUT = str(" (") + OUT + str(") ");
            for input_in in gene_list2chr(genesearch)[1:]:
                form.genesearch.data= str(input_in[1])
                form.chr.data= str(input_in[0])
            
                OUT = OUT + " UNION " + " (" + generatequery_gui_sub(form) + ") ";
    else:
        OUT = generatequery_gui_sub(form);

    if (int(limnum)==0): # add limit at the end;
        return str(OUT)
    else:
        return str(OUT) + " limit " + str(int(limnum))

    return OUT

if __name__ == '__main__':
    print generatequery_gui_sub(form)

