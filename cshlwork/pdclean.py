
import logging
import pandas as pd

def convert_columns(dataframe, col_list, lookup_dlists, dest_type, inplace=False ):
    '''
    Set values to those looked up from dict and set to dtype. 
    
    @param Dataframe to clean. 
    @param List of one or more columns. 
    @param List of value lookup dictionaries. 
    @param Destination type as *string* to ensure all are. 
    
    @return altered dataframe
        
    '''
    czip = zip (col_list, lookup_dlists)
    for (col, lookup)  in czip:
        def fixfunc(row):
            try:
                v = row[col]
                if isinstance(v, str):
                    v = v.strip()
                return lookup[row[col]]
            except:
                # if value is not handled, just return existing value. 
                # makes function idempotent. 
                return row[col]
            
        dataframe[col] = dataframe.apply(lambda row: fixfunc(row), axis=1)
        dataframe[col] = dataframe[col].astype(dest_type)
    return dataframe


def cast_columns(dataframe, col_list, dest_type, inplace=False ):
    '''
    Cast columns in col_list to selected type, setting
    any values un-convertable to NaNs.     
    @param Dataframe
    @param list of columns to be cast.
    @param destination type as string ( 'int32' '
    
    '''
    for col in col_list:
        def fixfunc(row):
            try:
                out = pd.to_numeric(row[col])
                return row[col]
            except:
                # if value is not handled, just return existing value. 
                # makes function idempotent. 
                return pd.NA
            
        dataframe[col] = dataframe.apply(lambda row: fixfunc(row), axis=1)
        dataframe[col] = dataframe[col].astype(dest_type)
    return dataframe
    
    
            
    
    
