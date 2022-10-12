import json
import os
import pandas as pd

#name = 'PTF12dam'


temp = os.path.abspath("second_cut.csv")
names = pd.read_csv(temp, sep=",")
names = pd.DataFrame(names)


temp = os.path.abspath("redshifts.csv")
used_r = pd.read_csv(temp, sep=",")
used_r = pd.DataFrame(used_r)

temp = os.path.abspath("table_of_sample_copy.csv")
table = pd.read_csv(temp, sep=",")
table = pd.DataFrame(table)

source = []

for name in names['Name']:
    r = float(used_r.loc[ used_r['Name'] == name]['z'])
    bands = table.loc[ table['Name'] == name]['bandpasses'].tolist()[0].split(', ')
    
    with open('sne/' + name + '.json') as f:
        sn = json.load(f)
    
    ############################################
    #поиск по z
    alias_list = []
    for cont in sn[name]['redshift']:
        if round(float(cont['value']), 4) == round(r, 4):
            for a in cont['source'].split(','):
                if a not in alias_list:
                    alias_list.append(a)
    
    source_str = ''
    for alias in alias_list:
        # for cont in sn[name]['sources']:
        #     if cont['alias'] == alias:
        #         source_list.append(cont['name'])
        #         continue
        try:
            source_str +=  sn[name]['sources'][int(alias) - 1]['bibcode'] + ', '
        except:
            continue
    
    ############################################
    #поиск по фотометрии
    alias_list = []
    for cont in sn[name]['photometry']:
        try:
            if cont['band'] in bands:
                for a in cont['source'].split(','):
                    if a not in alias_list:
                        alias_list.append(a)
        except:
            continue
        
        
    for alias in alias_list:
        try:
            source_str += sn[name]['sources'][int(alias) - 1]['bibcode'] + ', '
        except:
            continue
    
    source.append(source_str[:-2])

data = pd.DataFrame(columns=['name', 'sources'])
data['name'], data['sources'] = names, source

data.to_csv('references', index=None, sep='\t')