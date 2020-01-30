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
