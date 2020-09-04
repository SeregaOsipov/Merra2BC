import re
import config

# mapping between MERRA2 species and WRF species
pipe_to_process = []  # set of mapping rules from merra to wrf


def initialise():
    for a in config.spc_map:
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
            # TODO: this is reduntant field, should be moved to the coefficient
            rule_vo['wrf_exponent'] = m[2]
            pipe_to_process.append(rule_vo)

    print("\nConversion MAP:")
    for item in pipe_to_process:
        print(str(item) + ":\t" )


def get_list_of_wrf_spec_by_merra_var(name):
    return pipe_to_process.get(name)

