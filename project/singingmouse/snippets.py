#
#
#  Singing mouse analysis

# split query gene name. 
stdf[ ['probe','code','q_gene','end']] = stdf['name_read'].str.split(':',expand=True)
mmdf[ ['probe','code','q_gene','end']] = mmdf['name_read'].str.split(':',expand=True)

# split reference gene name
mmdf[ ['pre', 'r_gene']] = mmdf['name_align'].str.split('-',expand=True)




