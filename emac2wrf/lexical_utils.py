import re


def parse_mapping_rules(spc_map):
    pipe_to_process = []

    for a in spc_map:
        m = re.split('->|;', a)
        # print(m)
        ar = re.findall(r'(-?\ *\.?[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)\*\[?(\w+)\]?', m[1])
        # m=re.findall(r'(\w+) (\-\>)((-?\ *\.?[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)\*\[?(\w+)\]?', a[0])
        # http://stackoverflow.com/questions/18152597/extract-scientific-number-from-string
        m[0] = m[0].strip()
        m[2] = float(m[2])
        for r in ar:
            rule_vo = {}
            rule_vo['merra_key'] = r[1]
            rule_vo['wrf_key'] = m[0]
            rule_vo['wrf_multiplier'] = float(r[0])
            rule_vo['wrf_exponent'] = m[2]  # TODO: this is reduntant field, should be moved to the coefficient
            pipe_to_process.append(rule_vo)

    return pipe_to_process


def parse_mapping_rule(rule):
    pipe_to_process = []

    m = re.split('->|;', rule)
    # print(m)
    ar = re.findall(r'(-?\ *\.?[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)\*\[?(\w+)\]?', m[1])
    # m=re.findall(r'(\w+) (\-\>)((-?\ *\.?[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)\*\[?(\w+)\]?', a[0])
    # http://stackoverflow.com/questions/18152597/extract-scientific-number-from-string
    m[0] = m[0].strip()
    m[2] = float(m[2])
    for r in ar:
        rule_vo = {}
        rule_vo['merra_key'] = r[1]
        rule_vo['wrf_key'] = m[0]
        rule_vo['wrf_multiplier'] = float(r[0])
        rule_vo['wrf_exponent'] = m[2]  # TODO: this is redundant field, should be moved to the coefficient
        pipe_to_process.append(rule_vo)

    return pipe_to_process


def get_unique_wrf_keys_from_mappings(mappings):
    unique_wrf_keys = []
    for mapping in mappings:
        pipe_to_process = parse_mapping_rule(mapping.mapping_rule_str)
        for vo in pipe_to_process:
            if vo['wrf_key'] not in unique_wrf_keys:
                unique_wrf_keys.append(vo['wrf_key'])
    return unique_wrf_keys


def get_emac_poa_keys():
    # build the list of POA (primary organic aerosols in EMAC)
    types = 'f,bb'.split(',')
    subtypes = np.arange(1, 6)
    modes = 'ks,as,cs'.split(',')

    keys = []
    nco_string = ''
    for type in types:
        for subtype in subtypes:
            for mode in modes:
                key = '{}POA0{}_{}_ave'.format(type, subtype, mode)
                keys.append(key)
                nco_string += key + '+'

    return keys


