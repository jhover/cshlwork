#
# Sections from fastcafa removed for clarity.
#
#


def do_evaluate_pr(config, predictdf, goaspect,  threshold ):
    """
    i    cid           goterm       pest    cgid
    0    G960600000001 GO:0086041   53.0   Q9Y3Q4_HUMAN
    1    G960600000001 GO:0086089   49.0   Q9Y3Q4_HUMAN
    2    G960600000001 GO:0086090   49.0   Q9Y3Q4_HUMAN
    
    For each CID of predictions:
        # create correct numbers
        correctv = zeros()
        Get corresponding 'gid' from uniprot.  Q9Y3Q4_HUMAN
        For each goterm for that gid:
            Get goterm boolean vector gbv for that goterm
                correctv = correctv + gbv
        
        # process predictions (at this threshold?)
        predictv = zeros()
        For each goterm in prediction:
            Get goterm boolean vector gbvp for that goterm
                predictv = predictv + gbvp
        
        num-predicted = predictv.sum()
        num-correct = (correctv AND predictv).sum()
        num-annotated = correctv.sum() 
   
    cid    numpredict     numcorrect   numannotated    
    0001       27             18             37                   
    0002       30             30             36     
    
    #-----------------------------------------------
    #         totalpred     totalcorr    totalannot
    #
    #   % correct prediction     numcorrect / numpredict
    #   % correct answers        numcorrect / totalannot         
    #
    
    
    """
    logging.debug(f"goaspect={goaspect} threshold={threshold} predictdf=\n{predictdf}\n")

    outlist = []
    
    logging.debug("getting experimentally validated uniprot_byterm_df..")
    udf = get_uniprot_byterm_exponly_df(config, usecache=True)
    
    logging.debug("getting gomatrix...")
    gomatrix = get_ontology_matrix(config, usecache=True)
    logging.debug("getting goterm index for lookup...")
    gotermidx = get_gotermidx(config, usecache=True )
    
    cidlist = list(predictdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")

    for cid in cidlist:
        # build correct govector for all annotated terms in cid
        logging.debug(f"build correct govector for cid {cid}")
        cdf = predictdf[predictdf.cid == cid]
        # get gene ide. 
        cgid = cdf.cgid.unique()[0]
        #logging.debug(f"geneid for this target is is {cgid}")
        cgv = np.zeros(len(gotermidx), dtype=bool)
        # get true goterms for geneid
        cudf = udf[udf.pid == cgid]
        for (i, ser) in cudf.iterrows():
            # ser is the row. 
            #logging.debug(f"goterm is {ser.goterm}")  
            try:
                cgv = cgv + gomatrix[gotermidx[ser.goterm]].astype(np.int64)
            except KeyError:
                logging.debug("got keyerror for goterm index lookup. trying alt...")
                realgt = altidx[upser.goterm]
                logging.debug(f"real goterm is {realgt}")
                realidx = gotermidx[realgt]
                cgv2 = gomatrix[realidx].astype(np.int64)
                cgv = cgv + cgv2
        # cgv is now full vector of correct terms. 

        # build govector for prediction
        logging.debug(f"build prediction govector for cid {cid}")
        pgv = np.zeros(len(gotermidx), dtype=bool)
        pdf = predictdf[predictdf.cgid == cgid]
        for (i, ser) in pdf.iterrows():
            # ser is the row. 
            #logging.debug(f"goterm is {ser.goterm}")  
            try:
                pgv = pgv + gomatrix[gotermidx[ser.goterm]].astype(np.int64)
            except KeyError:
                logging.debug("got keyerror for goterm index lookup. trying alt...")
                realgt = altidx[upser.goterm]
                logging.debug(f"real goterm is {realgt}")
                realidx = gotermidx[realgt]
                pgv2 = gomatrix[realidx].astype(np.int64)
                pgv = pgv + pgv2
        # pgv is now full vector of predicted terms.         
        #num-predicted = predictv.sum()
        #num-correct = (correctv AND predictv).sum()
        #num-annotated = correctv.sum() 
        #cid    numpredict     numcorrect   numannotated    
        #0001       27             18             37                   
        
        numpredict = pgv.sum()
        numcorrect = np.logical_and(pgv,cgv).sum()
        numannotated = cgv.sum()
        outlist.append([cid, numpredict, numcorrect, numannotated])
    
    logging.debug(f"outlist={outlist}")
    return outlist


def do_evaluate_map(config, predictdf, goaspect):
    """
    i    cid           goterm       pest    cgid
    0    G960600000001 GO:0086041   53.0   Q9Y3Q4_HUMAN
    1    G960600000001 GO:0086089   49.0   Q9Y3Q4_HUMAN
    2    G960600000001 GO:0086090   49.0   Q9Y3Q4_HUMAN

    predictions:    p1    p2    p3     p4     p5    p6
    correct?        Y     N     Y      Y      N     Y 
    
    I               1     2     3      4      5     6 
    P              1/1    0    2/3    3/4     0     4/6                   
    AP              1     0    .666   .75     0     .666
    APS             1     1    1.666  2.416 2.416  3.08                   
    MAP             1 +  0 + .666 + .75 + 0 + .666  = .5138
    
    Return:
    0    cid           cgid             MAP
    1 G960600000001    Q9Y3Q4_HUMAN     0.467
    2 G960600000022    TF7L2_HUMAN      0.76 
    3 G960600000003
    
    
    """
    logging.debug(f"got predictdf {predictdf}evaluating...")
    ubgo = get_uniprot_bygene_object(config, usecache=True)
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug(f"got known uniprot and ontology object.")
    
    cidlist = list(predictdf.cid.unique())
    logging.debug(f"cid list: {cidlist}")
    
    listoflists = []
        
    for cid in cidlist:
        cdf = predictdf[predictdf.cid == cid]
        # get gene id. 
        cgid = cdf.cgid.unique()[0]
        #logging.debug(f"geneid for this target is is {cgid}")
        
        numcorrect = 0

        numseen = 0
        apsum = 0
        
        for (i, row) in cdf.iterrows():
            # ser is the row.
            #logging.debug(f"goterm is {row.goterm}")
            iscorrect = ubgo.contains(cgid, row.goterm) 
            numseen += 1 
            if iscorrect:
                numcorrect += 1
                prec = numcorrect / numseen 
                apsum = apsum + prec
        map = apsum / numseen  
        newrow = [ cid,  cgid, map, numseen ]
        listoflists.append(newrow)    
    
    df = pd.DataFrame(listoflists, columns=['cid','cgid','map', 'numgt'])
    return df






def run_evaluate_map(config, predictfile, outfile, goaspect=None):
    """
    Consume a prediction.csv file, and score based on accuracy. 
    X.prediction.csv
   
    """

    df = pd.read_csv(os.path.expanduser(predictfile), index_col=0)
    edf = do_evaluate_map(config, df, goaspect)
    logging.debug(f"got evaluation df:\n{edf}")
    edf.to_csv(outfile)
    mapsum = edf.map.sum()
    totalgt = edf.numgt.sum()
    logging.info(f"overall average MAP of prediction: {mapsum / totalgt}")
    
    max_goterms = config.get('global','max_goterms')
    eval_threshold = config.get('phmmer','eval_threshold')
    topx_threshold = config.get('phmmer','topx_threshold')
    score_method = config.get('phmmer','score_method')
    
    
    logging.info(f"hyperparams:\nmax_goterms={max_goterms}\neval_threshold={eval_threshold}\ntopx_threshold={topx_threshold}\nscore_method={score_method}  ")
    print(edf)



def build_uniprot_test(config, usecache):
    """
   
    [ {'proteinid': '001R_FRG3G', 
       'protein': '001R', 
       'species': 'FRG3G', 
       'proteinacc': 'Q6GZX4', 
       'taxonid': '654924', 
       'goterms': {'GO:0047043': 'IEA', 'GO:0006694': 'IEA'}, 
       'seqlength': 256, 
       'sequence': 'MAFSAEDVL......LYDDSFRKIYTDLGWKFTPL'},
       'gene' : 'LRRK2'
       .
       .
    ]
   
    Create redundant dataframe for later slimming. Include sequence. Cache. 
   
    """    
    cachedir = os.path.expanduser(config.get('uniprot','cachedir'))
    cachefile = f"{cachedir}/uniprottest.pickle"    
    lodt = None
    
    if os.path.exists(cachefile) and usecache:
        logging.debug("Cache hit. Using existing info...")    
        try:
            cf = open(cachefile, 'rb')    
            lodt = pickle.load(cf)
        except Exception:
            logging.error(f'unable to load via pickle from {cachefile}')
            traceback.print_exc(file=sys.stdout)    
        finally:
            cf.close()       
    else:    
        lod = build_uniprot(config, usecache=True)

        lodt = []
        for p in lod:
            newgts = {}
            for gt in p['goterms'].keys():
                evcode = p['goterms'][gt]
                item = [ p['proteinacc'],
                         p['protein'],
                         p['species'],
                         gt, 
                         evcode,
                         p['seqlength'],
                         p['sequence'] 
                      ]
                lodt.append(item)
                
        logging.debug(f"saving listofdicts: to {cachefile}")
        try:
            cf = open(cachefile, 'wb')    
            pickle.dump(lodt, cf )
            logging.debug(f"saved listofdicts: to {cachefile}")
        except Exception as e:
            logging.error(f'unable to dump via pickle to {cachefile}')
            traceback.print_exc(file=sys.stdout)     
        finally:
            cf.close() 
        logging.debug(f"created uniprot test source with {len(lodt)} entries.")
    return lodt


def build_uniprot_bypid(config):
    """
    Since this is only used by do_evaluate we will only use experimentally
    validated annotations. 
    
    builds dictionary geneid -> propagated govector. 
    
    Returns dict for later use..
    """
    ontobj = get_ontology_object(config, usecache=True)
    ubt = get_uniprot_byterm_df(config, usecache=True, exponly=True)
    pids = list(ubt.pid.unique())
    logging.debug(f"pid list, length={len(pids)} e.g. : {pids[:50]}")
    ubtd = ubt.to_dict(orient = 'index')
    logging.debug(f"converted DF to dict: e.g. \n{ [ubtd[i] for i in range(0,3)] } ")
   
    bypiddict = {}
    sumreport = 1
    suminterval = 10000
    repthresh = sumreport * suminterval
    gtlength = len(ontobj.gotermidx)
    logging.debug(f"gv length is {gtlength}")
    
    i = 0
    currentpid = None
    currentgv = np.zeros(gtlength, dtype=bool)
    while i < len(ubtd):
        row = ubtd[i]
        pid = row['pid']
        #logging.debug(f"row {i} is pid {pid}")
        if currentpid is None:
            currentpid = pid
            currentgv = currentgv + ontobj[row['goterm']]

        elif currentpid == pid:
            # same pid, continue...
            currentgv = currentgv + ontobj[row['goterm']] 
        else:
            # new pid
            bypiddict[currentpid] = currentgv 
            currentpid = pid
            currentgv = np.zeros(gtlength, dtype=bool)

        if len(bypiddict) >= repthresh:
            logging.info(f"Processed {len(bypiddict)} entries... ")
            sumreport +=1
            repthresh = sumreport * suminterval    
        
        i += 1
        
    samplekeys = list(bypiddict.keys())[:3]
    logging.debug(f"Made dict by proteinid: {[ bypiddict[k] for k in samplekeys]} ")
    return bypiddict


def build_uniprot_bygene(config):
    """
    Since this will be used for inference, we will use all evidence codes. 

        pacc species      goterm goev          pid   gene
0     P0DJZ0   PAVHV  GO:0030430  IDA    11K_PAVHV    11K
1     P32234   DROME  GO:0005525  IDA  128UP_DROME  128UP
2     P83011   SCYCA  GO:0043231  IDA  13KDA_SCYCA    NaN
7143  P24224   ECOLI  GO:0008897  IDA   ACPS_ECOLI   ACPS
7144  P24224   ECOLI  GO:0018070  IDA   ACPS_ECOLI   ACPS    

NOTES:  NaN expected for gene, omit...
        gene will not be unique. handle all...
    
    builds dictionary gene -> propagated govector. 
    
    Returns dict for later use..
    """
    ontobj = get_ontology_object(config, usecache=True)
    ubtdf = get_uniprot_byterm_df(config, usecache=True, exponly=False)
    ubtdf = ubtdf[ubtdf.gene.notna()]
    ubtdf.sort_values(by='gene', inplace=True)
    ubtdf.reset_index(drop=True, inplace=True)
        
    ubtd = ubtdf.to_dict(orient = 'index')
    logging.debug(f"converted DF to dict: e.g. \n{ [ubtd[i] for i in range(0,3)] } ")
   
    bygenedict = {}
    sumreport = 1
    suminterval = 10000
    repthresh = sumreport * suminterval
    gtlength = len(ontobj.gotermidx)
    logging.debug(f"gv length is {gtlength}")
    
    i = 0
    currentgid = None
    currentgv = np.zeros(gtlength, dtype=bool)
    while i < len(ubtd):
        row = ubtd[i]
        gid = row['gene']
        #logging.debug(f"row {i} is pid {pid}")
        if currentgid is None:
            currentgid = gid
            currentgv = currentgv + ontobj[row['goterm']]
            
        elif currentgid == gid:
            # same gid, continue...
            currentgv = currentgv + ontobj[row['goterm']] 
        else:
            # new gid
            bygenedict[currentgid] = currentgv 
            currentgid = gid
            currentgv = np.zeros(gtlength, dtype=bool)

        if len(bygenedict) >= repthresh:
            logging.info(f"Processed {len(bygenedict)} entries... ")
            sumreport +=1
            repthresh = sumreport * suminterval    
        i += 1
        
    samplekeys = list(bygenedict.keys())[:3]
    logging.debug(f"Made dict by geneid: {[ (k, bygenedict[k]) for k in samplekeys]} ")
    return bygenedict

