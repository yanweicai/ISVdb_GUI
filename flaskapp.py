import copy
import csv
import string
import random
import time

from flask import Flask, render_template, request, Response, stream_with_context
from urllib import quote,unquote


from model_v10 import InputForm_gui_v10
from query_v10 import generatequery_gui_v10

from model_v11 import InputForm_gui_v11
from query_v11 import generatequery_gui_v11

from model_v12 import InputForm_gui_v12
from query_v12 import generatequery_gui_v12

from model_v20 import InputForm_gui_v20
from query_v20 import generatequery_gui_v20

app = Flask(__name__)


def connection(): # Establish connection with ilvdb database/
    # These two modules exist in the server but can't get it installed locally
    import MySQLdb
    import MySQLdb.cursors # important
    
    conn = MySQLdb.connect(host="valdardb.its.unc.edu",user = "valdar_user",db="exon_1410E_11",passwd="u6GENwMGMPkDE6",cursorclass = MySQLdb.cursors.SSCursor)
    c = conn.cursor()

    return c, conn

def tablename(s):
    # change names in the result_header
    if 'gene_name' in s:
        s[s.index('gene_name')]='gene_id'
    if 'transcript_name' in s:
        s[s.index('transcript_name')]='transcript_id'
    if 'founder_name_1' in s:
        s[s.index('founder_name_1')]='haplotype_1'
    if 'founder_name_2' in s:
        s[s.index('founder_name_2')]='haplotype_2'
    if 'strain_name' in s:
        s[s.index('strain_name')]='strain'
    if 'strain_name_1' in s:
        s[s.index('strain_name_1')]='strain_1'
    if 'strain_name_2' in s:
        s[s.index('strain_name_2')]='strain_2'
    return(s)

def generate_d(c): # The function that take a connection variable and streamming MyQTL to download file
    # input c. connection that already executed QTL query
    result_header= [x[0] for x in c.description]
    result_header=tablename(result_header)
    yield ','.join(str(e) for e in result_header) + '\n' # yield header
    
    row = c.fetchone()
    while row is not None:
        yield ','.join(str(e) for e in row) + '\n'  # yield data in a memoryless way
        row = c.fetchone()

def strain_sep(result_table, strain_id):
    for i in range(0,len(result_table)):
        a = list(result_table[i])
        strnld=str(a[strain_id])
        if strnld=='A_J' or strnld =='CAST_EiJ' or strnld =='129S1_SvlmJ' or strnld =='NOD_ShiLtJ' or strnld =='NZO_HlLtJ' or strnld =='PWK_PhJ' or strnld =='WSB_EiJ' or strnld =='C57BL6J':
            continue
        else:
            a[strain_id] = strnld.split('_')[0]
            result_table[i] = a
    return result_table

########################### The Main GUI System #########################
@app.route('/getfile/') # Download route with <name> as the quoted version of query name
# No. ssh is not an option!
def getfile():
    import paramiko
    import select
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    #client.connect('valdardb.its.unc.edu',username='isvdbguest', password='sn$P2^Et5XpqZg')

    client.connect('valdardb.its.unc.edu',username='ywcai', password='Yan2015+1s')

    sftp_client = client.open_sftp()
    remote_file = sftp_client.open('/data/ywcai/dumpc/strain_list.txt')
    try:
        for line in remote_file:
            print line
    finally:
        remote_file.close()
    
    return "File Sucess"

def getfileconnet():
    import paramiko
    import select
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('valdardb.its.unc.edu',username='isvdbguest', password='sn$P2^Et5XpqZg')
    
    return client

def makefig(client,fname,form):
    #######################################################
    # make config file for R readin
    # client: a client object from getfileconnet()
    # fname: file name of config file
    # form: a form object from flask structure
    #######################################################
    client.exec_command("echo -e '"+ 'variantinputtype'+'\t'+form.variantinputtype.data+"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'view'+'\t'+form.view.data +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'variant_id'+'\t'+form.variant_id_text.data +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'chr'+'\t'+form.chr.data +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'strnum'+'\t'+str(form.strnum.data) + "' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'endnum'+'\t'+str(form.endnum.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'genesearch'+'\t'+form.genesearch.data +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'prob_cutoff'+'\t'+str(form.prob_cutoff.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'ismax'+'\t'+str(form.ismax.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'homo'+'\t'+str(form.homo.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'consequence_list'+'\t'+','.join(form.consequence.data) +"' >>"+fname+".cfig")
    
    client.exec_command("echo -e '"+ 'strainp1new'+'\t'+str(form.strainp1new.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'strainp2new'+'\t'+str(form.strainp2new.data) +"' >>"+fname+".cfig")

    client.exec_command("echo -e '"+ 'strainlistnew'+'\t'+','.join(form.strainlistnew.data) +"' >>"+fname+".cfig")
    #client.exec_command("echo -e '"+ 'selectCollum1'+'\t'+','.join(form.selectCollum1.data) +"' >>"+fname+".cfig")
    #client.exec_command("echo -e '"+ 'selectCollum2'+'\t'+','.join(form.selectCollum2.data) +"' >>"+fname+".cfig")
    #client.exec_command("echo -e '"+ 'selectCollum3'+'\t'+','.join(form.selectCollum3.data) +"' >>"+fname+".cfig")
    #client.exec_command("echo -e '"+ 'selectCollum4'+'\t'+','.join(form.selectCollum4.data) +"' >>"+fname+".cfig")
    client.exec_command("echo -e '"+ 'limnum'+'\t'+str(form.limnum.data) +"' >>"+fname+".cfig")

@app.route('/v2.0',methods=['GET','POST'])
def gui_page_v20():
    
    client = getfileconnet()
    form = InputForm_gui_v20(request.form)
    result='SUCCESS RUN'
    if request.method == 'POST':
        # randomly generate a filename that only exists in memory.
        fname="DB"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        # creat a configure file with unique name in valderdb server
        makefig(client,fname,form)
        time.sleep(1) # wait for configure file to make
    
        # in the valdardb server, make the result table by data.table filtering and joining
        stdin,stdout,stderr=client.exec_command('/home/doreper/bin/Rscript makefile.r '+fname)
        time.sleep(20)
        se='FILE not exist'
        i=0;
        while (se!='' and i<=100): # Check if the file exists
            time.sleep(10)
            stdin,stdout,stderr=client.exec_command('ls -l '+fname+'.res') # if there are result file
            se=stderr.read();
            i=i+1
    else:
        fname="DFLT" # use default result table

    result_table=list()

    remote_file = client.open_sftp().open(fname+'.res')
    i=1
    for line in remote_file:
        if i==1:
            result_header=line.split(',');
        else:
            result_table.append(line.split(','));
        i=i+1

    #result="FILE haven't been finished"

    return render_template("index_v20.html", form=form, result=result, result0q=fname, result_table=result_table, result_header=result_header)

@app.route('/downloadf/<name>') # Download route with <name> as the quoted version of query name
def downloadf(name):
    import os
    from tempfile import SpooledTemporaryFile
    # find if there exits configure file
    client = getfileconnet()
    transfer=client.open_sftp()
    path=name+'.res'

    # First, run the unlimited version make result file
    #stdin,stdout,stderr=client.exec_command('/home/doreper/bin/Rscript makefile.r '+fname)

    # Check if the file exists
    stdin1,stdout1,stderr1=client.exec_command('ls -l '+name+'.res') # if there are result file
    stdin2,stdout2,stderr2=client.exec_command('ls -l '+name+'.cfig') # if there are cfig file
    se1=stderr1.read();
    se2=stderr2.read();
    if (se1!='' and se2==''):
        return "Task still running"
    if (se1!='' and se2!=''):
        return "unknown task ID"
    if (se1=='' and se2==''):
        with SpooledTemporaryFile(1024000) as f:  # example max size before moving to file = 1MB
            transfer.getfo(path, f)
            f.seek(0)
            r = app.response_class(f.read(), mimetype='application/octet-stream')
    
        r.headers.set('Content-Disposition', 'attachment', filename=os.path.basename(path))
        return r

@app.route('/download/<name>') # Download route with <name> as the quoted version of query name
def download(name):
    c, conn = connection()
    try:
        c.execute(unquote(name))
        return Response(stream_with_context(generate_d(c)), mimetype='text/csv')
    except:
        return "Unknown query"

@app.route('/v1.0',methods=['GET','POST'])
def gui_page():
    c, conn = connection()
    form = InputForm_gui_v10(request.form)
    if request.method == 'POST':
        result=generatequery_gui_v10(variantinputtype = form.variantinputtype.data, view = form.view.data,
                                 variant_id =form.variant_id_text.data,
                                 chr = form.chr.data, strnum = form.strnum.data, endnum = form.endnum.data,
                                 genesearch=form.genesearch.data,
                                 selectCollum1=form.selectCollum1.data, selectCollum2=form.selectCollum2.data, selectCollum3=form.selectCollum3.data, selectCollum4= form.selectCollum4.data,
                                 prob_cutoff= form.prob_cutoff.data, ismax= form.ismax.data, homo = form.homo.data, consequence_list=form.consequence.data,
                                 strainlistnew = form.strainlistnew.data, strainp1new=form.strainp1new.data, strainp2new=form.strainp2new.data,
                                 limnum=form.limnum.data)
    else:
        result = "SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM exon_1410b.genotype_transcript_view WHERE chrom='19' and pos>=30000000 and pos<=35000000 and ( strain_id = 1 OR strain_id = 2 OR strain_id = 3 OR strain_id = 4 OR strain_id = 5 OR strain_id = 6 OR strain_id = 7 OR strain_id = 8 OR strain_id = 9 ) limit 100";
    c.execute(str(result))
    result_table = c.fetchall()
    result_header= [x[0] for x in c.description]
    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query
    result_header=tablename(result_header)

    return render_template("index_v10.html", form=form, result=result, result0q=result0q,result_table=result_table, result_header=result_header)

@app.route('/help/',methods=['GET','POST'])
def help_my():
    return render_template("help_ilvdb.html")

@app.route('/data/',methods=['GET','POST'])
def data_my():
    return render_template("data_ilvdb.html")

@app.route('/archive/',methods=['GET','POST'])
def archive_my():
    return render_template("archive_ilvdb.html")

@app.route('/',methods=['GET','POST'])
def gui_page_test():
    
    form = InputForm_gui_v11(request.form)
    c, conn = connection()
    if request.method == 'POST':
        result=generatequery_gui_v11(variantinputtype = form.variantinputtype.data, view = form.view.data,
                                 variant_id =form.variant_id_text.data,
                                 chr = form.chr.data, strnum = form.strnum.data, endnum = form.endnum.data,
                                 genesearch=form.genesearch.data,
                                 selectCollum1=form.selectCollum1.data, selectCollum2=form.selectCollum2.data, selectCollum3=form.selectCollum3.data, selectCollum4= form.selectCollum4.data,
                                 prob_cutoff= form.prob_cutoff.data, ismax= form.ismax.data, homo = form.homo.data, consequence_list=form.consequence.data,
                                 strainlistnew = form.strainlistnew.data, strainp1new=form.strainp1new.data, strainp2new=form.strainp2new.data,
                                 limnum=form.limnum.data)
    else:
        result = "SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM exon_1410E_11.genotype_transcript_view WHERE chrom='11' and pos>=30000000 and pos<=35000000 and ( strain_id = 1 OR strain_id = 2 OR strain_id = 3 OR strain_id = 4 OR strain_id = 5 OR strain_id = 6 OR strain_id = 7 OR strain_id = 8 OR strain_id = 9 ) limit 100";

    c.execute(str(result))
    result_table = list(c)
    result_header= [x[0] for x in c.description]

    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query
    result_header=tablename(result_header)

    result_table_b = copy.deepcopy(result_table) # Important! Must copy to list tuple to change it!
    if 'prob' in result_header: # chop prob to 0.000 total 5 charactors
        prob_id = result_header.index('prob')
        for i in range(0,len(result_table)):
            a = list(result_table_b[i])
            a[prob_id] = str(a[prob_id])[:5]
            result_table_b[i]=a
    if 'strain' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain'))
    if 'strain_1' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_1'))
    if 'strain_2' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_2'))


    return render_template("index_v11.html", form=form, result=result, result0q=result0q, result_table=result_table_b, result_header=result_header)


@app.route('/v1.2',methods=['GET','POST'])
def gui_page_v12():
    
    form = InputForm_gui_v12(request.form)
    c, conn = connection()
    if request.method == 'POST':
        result=generatequery_gui_v12(form)
    else:
        result = "SELECT distinct variant_id,chrom,pos,strain_name,allele_1,allele_2,prob,is_max,consequence_1,consequence_2,gene_name,transcript_name FROM exon_1410D_11.genotype_transcript_view WHERE chrom='11' and pos>=30000000 and pos<=35000000 and ( strain_id = 1 OR strain_id = 2 OR strain_id = 3 OR strain_id = 4 OR strain_id = 5 OR strain_id = 6 OR strain_id = 7 OR strain_id = 8 OR strain_id = 9 ) limit 100";

    result0 = result.split(" limit ",1)[0]  # The result0 query doesnot set limit
    result0q= quote(result0) # URLable version of query

    c.execute(str(result))
    result_table = list(c)
    result_header= [x[0] for x in c.description]
    result_header=tablename(result_header)

    result_table_b = copy.deepcopy(result_table) # Important! Must copy to list tuple to change it!
    if 'prob' in result_header: # chop prob to 0.000 total 5 charactors
        prob_id = result_header.index('prob')
        for i in range(0,len(result_table)):
            a = list(result_table_b[i])
            a[prob_id] = str(a[prob_id])[:5]
            result_table_b[i]=a
    if 'strain' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain'))
    if 'strain_1' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_1'))
    if 'strain_2' in result_header:
        result_table = strain_sep(result_table_b, result_header.index('strain_2'))

    return render_template("index_v12.html", form=form, result=result, result0q=result0q, result_table=result_table_b, result_header=result_header)

if __name__ == "__main__":
	app.run(debug=True)

