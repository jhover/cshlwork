#!/usr/bin/env python
#
#  Quick n dirty. 
#  Takes .CSV made from NCBI all taxon info...
# 
# ,species,kingdom,taxonid,lineanname,commonname
# 0,AADNV,V,648330,Aedes albopictus densovirus (isolate Boublik/1994),AalDNV
# 1,AAV2 ,V,10804,Adeno-associated virus 2,AAV-2
# 2,AAV2S,V,648242,Adeno-associated virus 2 (isolate Srivastava/1982),AAV-2
# 3,ABAMA,E,118452,Abacion magnum,Millipede
# 4,ABANI,E,72259,Abaeis nicippe,Sleepy orange butterfly
#
# Outputs only proper linnean names, no dupes.  
#
import argparse
import logging
import traceback
import sys

def parse_latinnames(infile):
    try:
        filehandle = open(infile, 'r')
        for line in filehandle:
            logging.debug(f"handling line {line}")
            flist = line.split(",")
            idx = flist[0]
            tcode = flist[1]                
            kingdom = flist[2]
            taxid = flist[3]
            rawname = flist[4]
            commonname = flist[5]    
            
            sflist = rawname.split()
            if len(sflist) < 2:
                pass
            elif len(sflist) > 2:
                pass
            else:
                fixed = rawname
                logging.debug(f"raw: {rawname}")
                logging.debug(f"fixed: {fixed} ")
                print(fixed)
            
            
            
               
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
        
    finally:
        if filehandle is not None:
            filehandle.close()
    logging.debug("done")
    


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='Pandas .CSV')
    
                    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
        
    logging.info("dendrotool...")
    
    parse_latinnames(args.infile)