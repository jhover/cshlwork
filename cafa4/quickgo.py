#!/usr/bin/env python
#
#  Look up GO annoatations from UniProtKB:entrys. 
#  E.g.  P35213  (1422B_RAT  gene: Ywhab
# 

import requests, sys

annotationURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/about"

def print_body(req):
    
    if not req.ok:
      req.raise_for_status()
      sys.exit()

    responseBody = r.text
    print(responseBody)


r = requests.get(annotationURL, headers={ "Accept" : "application/json"})

print_body(r)

requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?targetSet=Ywhab";

headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
r = requests.post(requestURL, headers=headers, data='{  "and": {"goTerms": ["string"], "goUsage": "string", "goUsageRelationships": "string"},  "not": { "goTerms": ["string"], "goUsage": "string", "goUsageRelationships": "string"  }}')

print_body(r)





