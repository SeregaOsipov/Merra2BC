def get_cams2wrf_mapping():
    '''
    chem 100/106 case (MADE SORGAM)

    To prescribe aerosol in the MADE scheme we need to specify two parameters:
    mass concentration (linked to 3rd moment) and number concentration (0th moment).
    3rd moment (mass) can be derived from MERRA2. For 0th moment (number concentration) I assumed fixed size.

    See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient


    Unit Conversion: CAMS provides trace gases and aerosols in mass mixing ratio ($kg/kg$). WRF-Chem usually expects volume mixing ratio ($ppmv$) for trace gases and mass mixing ratio ($ukg/kg$) for aerosols.
    '''


    ml_stream = 'ml'
    mappings = [
        # --- Inorganics (Gases) ---
        # Ozone: 28.97 / 48.0
        IcBcMappingItem('o3 -> 0.604*[go3];1.e6', ml_stream),
        # Carbon Monoxide: 28.97 / 28.01
        IcBcMappingItem('co -> 1.034*[co];1.e6', ml_stream),
        # Nitric Oxide: 28.97 / 30.01
        IcBcMappingItem('no -> 0.965*[no];1.e6', ml_stream),
        # Nitrogen Dioxide: 28.97 / 46.01
        IcBcMappingItem('no2 -> 0.630*[no2];1.e6', ml_stream),
        # Sulfur Dioxide: 28.97 / 64.06
        IcBcMappingItem('so2 -> 0.452*[so2];1.e6', ml_stream),
        # Nitric Acid: 28.97 / 63.01
        IcBcMappingItem('hno3 -> 0.460*[hno3];1.e6', ml_stream),
        # Hydrogen Peroxide: 28.97 / 34.01
        IcBcMappingItem('h2o2 -> 0.852*[h2o2];1.e6', ml_stream),
        # Methane: 28.97 / 16.04
        IcBcMappingItem('ch4 -> 1.806*[ch4_c];1.e6', ml_stream),
        # Ammonia: 28.97 / 17.03
        IcBcMappingItem('nh3 -> 1.701*[nh3];1.e6', ml_stream),
        # Organic Nitrates (ONIT): 28.97 / 119.0 (assuming average weight for RADM2 lumping)
        IcBcMappingItem('onit -> 0.243*[onit];1.e6', ml_stream),
        # --- Nighttime / Reservoir Species ---
        # Dinitrogen Pentoxide: 28.97 / 108.01
        IcBcMappingItem('n2o5 -> 0.268*[n2o5];1.e6', ml_stream),
        # Pernitric Acid: 28.97 / 79.01
        IcBcMappingItem('hno4 -> 0.367*[ho2no2];1.e6', ml_stream),
        # DMS (Dimethyl Sulfide): 28.97 / 62.13
        # IcBcMappingItem('dms -> 0.466*[dms];1.e6', ml_stream),  # DMS is not in the chem_opt=100
        # PAN (Peroxyacetyl nitrate): 28.97 / 121.05
        IcBcMappingItem('pan -> 0.239*[pan];1.e6', ml_stream),
        # NO3 Radical (Gas phase): 28.97 / 62.00
        IcBcMappingItem('no3 -> 0.467*[no3];1.e6', ml_stream),
        # Hydroxyl Radical (OH): 28.97 / 17.01
        IcBcMappingItem('ho -> 1.703*[oh];1.e6', ml_stream),
        # Hydroperoxy Radical (HO2): 28.97 / 33.01
        IcBcMappingItem('ho2 -> 0.878*[ho2];1.e6', ml_stream),
        # Methyl Peroxy Radical (CH3O2): 28.97 / 47.03
        # IcBcMappingItem('mo2 -> 0.616*[ch3o2];1.e2', ml_stream),  # MO2 is not in the chem_opt=100

        # --- VOCs (Lumped RADM2) ---
        # Formaldehyde: 28.97 / 30.03
        IcBcMappingItem('hcho -> 0.965*[hcho];1.e6', ml_stream),
        # Acetaldehyde: 28.97 / 44.05
        IcBcMappingItem('ald -> 0.658*[ald2];1.e6', ml_stream),
        # Isoprene: 28.97 / 68.12
        IcBcMappingItem('iso -> 0.425*[c5h8];1.e6', ml_stream),
        # Ethane: 28.97 / 30.07
        IcBcMappingItem('eth -> 0.963*[c2h6];1.e6', ml_stream),
        # Propane (HC3): 28.97 / 44.10
        IcBcMappingItem('hc3 -> 0.657*[c3h8];1.e6', ml_stream),
        # Ethene (OL2): 28.97 / 28.05
        # IcBcMappingItem('ol2 -> 1.033*[c2h4];1.e6', ml_stream),  # ol2 is not in the chem_opt=100. DO not confuse ol2 and OL2, different variables
        # Propene (OLT): 28.97 / 42.08 + # Alkenes (OLE): 28.97 / 30.0 (Average lumped Alkenes)
        IcBcMappingItem('olt -> 0.688*[c3h6]+0.965*[ole];1.e6', ml_stream),
        # Acetone (KET): 28.97 / 58.08  +  # Ethanol (Mapped to KET or ALD depending on lumping preference): 28.97 / 46.07
        IcBcMappingItem('ket -> 0.499*[ch3coch3]+0.628*[c2h5oh];1.e6', ml_stream),
        # Formic Acid (ORA2): 28.97 / 46.03 + # Methanol: 28.97 / 32.04 (Often mapped to ORA2 in RADM2) + # Methylcarboxylic acid (MCOOH): 28.97 / 60.05
        IcBcMappingItem('ora2 -> 0.629*[hcooh]+0.904*[ch3oh]+0.482*[mcooh];1.e6', ml_stream),
        # Paraffins (Lumped Alkanes): 28.97 / 1.0 (PAR is often treated as unitless/molar based in CAMS, but for RADM2 lumping we use 1.0 or the HC5/HC8 equivalent)
        IcBcMappingItem('hc5 -> 1.0*[par];1.e6', ml_stream),
        # Methylglyoxal (MGLY): 28.97 / 72.06
        IcBcMappingItem('mgly -> 0.402*[ch3cocho];1.e6', ml_stream),
        # Methyl Hydroperoxide (MOP): 28.97 / 48.04
        # IcBcMappingItem('mop -> 0.603*[ch3ooh];1.e6', ml_stream),  # mop is not in the chem_opt=100
        # Monoterpenes (C10H16): 28.97 / 136.23
        IcBcMappingItem('api -> 0.213*[c10h16];1.e6', ml_stream),
        # Methanesulfonic acid (MSA): 28.97 / 96.11
        # IcBcMappingItem('msa -> 0.301*[msa];1.e6', ml_stream),  # MSA is not in the chem_opt=100

        # # --- Aerosols (GOCART Bins) ---
        # Sulfate (Bulk)
        IcBcMappingItem('ac0 -> 3.56*[aermr11];1.e11', ml_stream),
        IcBcMappingItem('so4aj -> 1.0*[aermr11];1.e9', ml_stream),
        # Dust Bins [0.03-0.55 um; 0.55-0.9 um; 0.9-20 um], put all dust into the coarse mode
        IcBcMappingItem('corn -> 1.68*[aermr04]+1.68*[aermr05]+1.68*[aermr06];1.e7', ml_stream),
        IcBcMappingItem('soila -> 1.0*[aermr04]+1.0*[aermr05]+1.0*[aermr06];1.e9', ml_stream),
        # Sea Salt Bins [0.03-0.5 um; 0.5-5 um; 5-20 um], put all sea salt into the coarse mode.
        IcBcMappingItem('corn -> 1.98*[aermr01]+1.98*[aermr02]+1.98*[aermr03];1.e7', ml_stream),
        IcBcMappingItem('seas -> 1.0*[aermr01]+1.0*[aermr02]+1.0*[aermr03];1.e9', ml_stream),
        # Black Carbon (Hydrophilic and Hydrophobic)  # put all BC into acc mode  # put BC into elemental carbon, consistent with emissions
        IcBcMappingItem('ac0 -> 6.41*[aermr09]+6.41*[aermr10];1.e11', ml_stream),
        IcBcMappingItem('ecj -> 1.0*[aermr09]+1.0*[aermr10];1.e9', ml_stream),
        # Organic Carbon (Hydrophilic and Hydrophobic organic matter), assign into primary orgpaj in WRF, similar to the emissions
        IcBcMappingItem('ac0 -> 6.41*[aermr07]+6.41*[aermr08];1.e11', ml_stream),
        IcBcMappingItem('orgpaj -> 1.0*[aermr07]+1.0*[aermr08];1.e9', ml_stream),
        # NH4, assign into acc mode in WRF. Strictly speaking one have to split between Aitken and Acc
        IcBcMappingItem('ac0 -> 3.7*[nh4];1.e11', ml_stream),
        IcBcMappingItem('nh4aj -> 1.0*[nh4];1.e9', ml_stream),
        # NO3, assign into acc mode in WRF.
        IcBcMappingItem('ac0 -> 4.24*[no3_a];1.e11', ml_stream),
        IcBcMappingItem('no3aj -> 1.0*[no3_a];1.e9', ml_stream),
    ]

    # mappings = [  # debug
    #     IcBcMappingItem('so2 -> 0.452*[so2];1.e6', ml_stream),
    #     # IcBcMappingItem('soila -> 1.0*[aermr04]+1.0*[aermr05]+1.0*[aermr06];1.e9', ml_stream),
    # ]
    return mappings


def get_emac2wrf_mapping():  # this is the new mapping of streams, where Andrea put everything specifically into WRF_bc_met & WRF_bc_chem
    '''
    see channel.nml
    ECHAM5 met (tm1_ave, pressi,geopoti_ave,
    gboxarea_ave, grvol_ave [m^3], rho_air_dry_ave[kg m^-3]
    u & v could be cos(phi) weighted
    aer_oracle.nc
    tr_terp.nc (terpines)
    '''

    '''
    chem_100 case (MADE SORGAM)
    Prescribe the organic aerosols from EMAC
    See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient

    Below is some communication with Andrea Pozzer regarding the VOCs/SOA.
    "SOAv0201_* is the aerosol phase equivalent of the first entry in SOGv02 in oracle.nml, i.e. LaSOGv02. This is a product of anthropogenic gases (see final part of ~/MESSy/Nussbaumer/messy_2.54.0/messy/mbm/caaba/mecca/mecca.eqn).
    SOAv0202_* in instead the equivalent of the second entry of SOGv02, i.e. LbSOGv02, which are the *BIOGENIC* emissions."
    '''

    mappings = [  # see channel.nml
        IcBcMappingItem('so2 -> 1*[SO2];1.e6', 'WRF_bc_chem'),  # SO2_ave
        IcBcMappingItem('ch4 -> 1*[CH4];1.e6', 'WRF_bc_chem'),  # _ave
        IcBcMappingItem('co -> 1*[CO];1.e6', 'WRF_bc_chem'),
        IcBcMappingItem('o3 -> 1*[O3];1.e6', 'WRF_bc_chem'),
        # TODO: I can map H2SO4 too, but the life time is short compared to SO4, I can ignore it
        # VOCs, but currently only longer lived ones
        IcBcMappingItem('hono -> 1*[HONO];1.e6', 'WRF_bc_chem'),
        IcBcMappingItem('hcho -> 1*[HCHO];1.e6', 'WRF_bc_chem'),
        IcBcMappingItem('eth -> 1*[C2H6];1.e6', 'WRF_bc_chem'),  # ethane
        IcBcMappingItem('hc3 -> 1*[C3H8];1.e6', 'WRF_bc_chem'),  # propane in Volatility Basis Set

        # aerosols, see aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
        # aerosols in EMAC are in ppv, in WRF in mmr.
        # Compute mmr = vmr * Ma/Mb, where Mb is dry air molar mass = 28.97 g/mol
        # Check EMAC nc for molar masses & densities

        # Sulfate, split into nu + acc
        # TODO: what is mm in SO4mm_?
        IcBcMappingItem('nu0 -> 2.99*[SO4mm_ns]+2.99*[SO4mm_ks];1.e14', 'WRF_bc_chem'),
        IcBcMappingItem('so4ai -> 3.315*[SO4mm_ns]+3.315*[SO4mm_ks];1.e9', 'WRF_bc_chem'),  # 96.06/28.97 = 3.315
        IcBcMappingItem('ac0 -> 3.56*[SO4mm_as]+3.56*[SO4mm_cs];1.e11', 'WRF_bc_chem'),
        IcBcMappingItem('so4aj -> 3.315*[SO4mm_as]+3.315*[SO4mm_cs];1.e9', 'WRF_bc_chem'),
        # Dust, only coarse
        IcBcMappingItem('corn -> 1.65*[DU_as]+1.65*[DU_ai]+1.65*[DU_cs]+1.65*[DU_ci];1.e7', 'WRF_bc_chem'),
        IcBcMappingItem('soila -> 1.38*[DU_as]+1.38*[DU_ai]+1.38*[DU_cs]+1.38*[DU_ci];1.e9', 'WRF_bc_chem'),  # 40.08/28.97=1.38
        # Sea salt
        IcBcMappingItem('corn -> 2.01*[SS_ks]+2.01*[SS_as]+2.01[SS_cs];1.e7', 'WRF_bc_chem'),  # coarse
        IcBcMappingItem('seas -> 2.017*[SS_ks]+2.017*[SS_as]+2.017*[SS_cs];1.e9', 'WRF_bc_chem'),  # SS molar mass is 58.44; 58.44/28.97=2.017
        # TODO: Cations: Mg, Ca and such (from Dust and Sea Salt)
        # Organics, BC. BC probably should be in the nucleation mode, not acc
        IcBcMappingItem('ac0 -> 3.2*[BC_ks]+3.2*[BC_as]+3.2*[BC_cs]+3.2*[BC_ki];1.e11', 'WRF_bc_chem'),  # accumulation
        IcBcMappingItem('ecj -> 0.414*[BC_ks]+0.414*[BC_as]+0.414*[BC_cs]+0.414*[BC_ki];1.e9', 'WRF_bc_chem'),  # molar mass 12.01 =>

        # Organics, POA. Check EMAC in nml/oracle.nml.
        # TODO: this rule is slow due to large number of variables. USE nco to preprocess and sum up all POA. Same for SOA by bins
        IcBcMappingItem('ac0 -> 6.41*[fPOA01_ks]+6.41*[fPOA02_ks]+6.41*[fPOA03_ks]+6.41*[fPOA04_ks]+6.41*[fPOA05_ks]+'
                        '6.41*[fPOA01_as]+6.41*[fPOA02_as]+6.41*[fPOA03_as]+6.41*[fPOA04_as]+6.41*[fPOA05_as]+'
                        '6.41*[bbPOA01_ks]+6.41*[bbPOA02_ks]+6.41*[bbPOA03_ks]+6.41*[bbPOA04_ks]+'
                        '6.41*[bbPOA01_as]+6.41*[bbPOA02_as]+6.41*[bbPOA03_as]+6.41*[bbPOA04_as];1.e11', 'WRF_bc_chem'),
        IcBcMappingItem('orgpaj -> 8.63*[fPOA01_ks]+8.63*[fPOA02_ks]+8.63*[fPOA03_ks]+8.63*[fPOA04_ks]+8.63*[fPOA05_ks]+'
                        '8.63*[fPOA01_as]+8.63*[fPOA02_as]+8.63*[fPOA03_as]+8.63*[fPOA04_as]+8.63*[fPOA05_as]+'
                        '8.63*[bbPOA01_ks]+8.63*[bbPOA02_ks]+8.63*[bbPOA03_ks]+8.63*[bbPOA04_ks]+'
                        '8.63*[bbPOA01_as]+8.63*[bbPOA02_as]+8.63*[bbPOA03_as]+8.63*[bbPOA04_as];1.e9', 'WRF_bc_chem'),  # molar mass 250 =>
        # Organics, SOA, have to assign according to VBS bins (saturation vapor pressure).
        # EMAC vars: *sv semivolatile, *iv intermidiate, *v volatile
        # WRF has 4 bins from 1-10**3 ug/m^3 @ 300K. EMAC bins are in oracle.nml
        # SOA, VBS bin 1, asoa1j. An example only
        # IcBcMappingItem('ac0 -> 6.41*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e11', 'WRF_bc_chem'),
        # IcBcMappingItem('orgpaj -> 8.63*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e9', 'WRF_bc_chem'),  # molar mass 250 =>

        # Nitrate
        IcBcMappingItem('nu0 -> 3.56*[NO3m_ns]+3.56*[NO3m_ks];1.e14', 'WRF_bc_chem'),
        IcBcMappingItem('no3ai -> 2.14*[NO3m_ns]+2.14*[NO3m_ks];1.e9', 'WRF_bc_chem'),  # 62.01/28.97
        IcBcMappingItem('ac0 -> 4.24*[NO3m_as]+4.24*[NO3m_cs];1.e11', 'WRF_bc_chem'),
        IcBcMappingItem('no3aj -> 2.14*[NO3m_as]+2.14*[NO3m_cs];1.e9', 'WRF_bc_chem'),
        # Ammonium
        IcBcMappingItem('nu0 -> 3.11*[NH4p_ns]+3.11*[NH4p_ks];1.e14', 'WRF_bc_chem'),
        IcBcMappingItem('nh4ai -> 0.62*[NH4p_ns]+0.62*[NH4p_ks];1.e9', 'WRF_bc_chem'),  # 18.05/28.97
        IcBcMappingItem('ac0 -> 3.70*[NH4p_as]+3.70*[NH4p_cs];1.e11', 'WRF_bc_chem'),
        IcBcMappingItem('nh4aj -> 0.62*[NH4p_as]+0.62*[NH4p_cs];1.e9', 'WRF_bc_chem'),
    ]
    return mappings


# def get_emac2wrf_mapping_OUTDATED():  # original mapping of EMAC streams
#
#     '''
#     see channel.nml
#     ECHAM5 met (tm1_ave, pressi,geopoti_ave,
#     gboxarea_ave, grvol_ave [m^3], rho_air_dry_ave[kg m^-3]
#     u & v could be cos(phi) weighted
#     aer_oracle.nc
#     tr_terp.nc (terpines)
#     '''
#
#     '''
#     chem_100 case (MADE SORGAM)
#     Prescribe the organic aerosols from EMAC
#     See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
#
#     Below is some communication with Andrea Pozzer regarding the VOCs/SOA.
#     "SOAv0201_* is the aerosol phase equivalent of the first entry in SOGv02 in oracle.nml, i.e. LaSOGv02. This is a product of anthropogenic gases (see final part of ~/MESSy/Nussbaumer/messy_2.54.0/messy/mbm/caaba/mecca/mecca.eqn).
#     SOAv0202_* in instead the equivalent of the second entry of SOGv02, i.e. LbSOGv02, which are the *BIOGENIC* emissions."
#     '''
#
#
#     mappings = [  # see channel.nml
#         IcBcMappingItem('so2 -> 1*[SO2];1.e6', 'tr_sulfur'),  # SO2_ave
#         IcBcMappingItem('ch4 -> 1*[CH4];1.e6', 'tr_hycarbs'),  # _ave
#         IcBcMappingItem('co -> 1*[CO];1.e6', 'tr_alks'),
#         IcBcMappingItem('o3 -> 1*[O3];1.e6', 'tr_Ox_HOx'),
#         # TODO: I can map H2SO4 too, but the life time is short compared to SO4, I can ignore it
#         # VOCs, but currently only longer lived ones
#         IcBcMappingItem('hono -> 1*[HONO];1.e6', 'tr_NOx_NOy'),
#         IcBcMappingItem('hcho -> 1*[HCHO];1.e6', 'tr_alks'),
#         IcBcMappingItem('eth -> 1*[C2H6];1.e6', 'tr_hycarbs'),  # ethane
#         IcBcMappingItem('hc3 -> 1*[C3H8];1.e6', 'tr_hycarbs'),  # propane in Volatility Basis Set
#
#         # aerosols, see aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
#         # aerosols in EMAC are in ppv, in WRF in mmr.
#         # Compute mmr = vmr * Ma/Mb, where Mb is dry air molar mass = 28.97 g/mol
#         # Check EMAC nc for molar masses & densities
#
#         # Sulfate, split into nu + acc
#         # TODO: what is mm in SO4mm_?
#         IcBcMappingItem('nu0 -> 2.99*[SO4mm_ns]+2.99*[SO4mm_ks];1.e14', 'aer_sulf'),
#         IcBcMappingItem('so4ai -> 3.315*[SO4mm_ns]+3.315*[SO4mm_ks];1.e9', 'aer_sulf'),    # 96.06/28.97 = 3.315
#         IcBcMappingItem('ac0 -> 3.56*[SO4mm_as]+3.56*[SO4mm_cs];1.e11', 'aer_sulf'),
#         IcBcMappingItem('so4aj -> 3.315*[SO4mm_as]+3.315*[SO4mm_cs];1.e9', 'aer_sulf'),
#         # Dust, only coarse
#         IcBcMappingItem('corn -> 1.65*[DU_as]+1.65*[DU_ai]+1.65*[DU_cs]+1.65*[DU_ci];1.e7', 'aer_du_ss'),
#         IcBcMappingItem('soila -> 1.38*[DU_as]+1.38*[DU_ai]+1.38*[DU_cs]+1.38*[DU_ci];1.e9', 'aer_du_ss'),  # 40.08/28.97=1.38
#         # Sea salt
#         IcBcMappingItem('corn -> 2.01*[SS_ks]+2.01*[SS_as]+2.01[SS_cs];1.e7', 'aer_du_ss'),  # coarse
#         IcBcMappingItem('seas -> 2.017*[SS_ks]+2.017*[SS_as]+2.017*[SS_cs];1.e9', 'aer_du_ss'),  # SS molar mass is 58.44; 58.44/28.97=2.017
#         # TODO: Cations: Mg, Ca and such (from Dust and Sea Salt)
#         # Organics, BC. BC probably should be in the nucleation mode, not acc
#         IcBcMappingItem('ac0 -> 3.2*[BC_ks]+3.2*[BC_as]+3.2*[BC_cs]+3.2*[BC_ki];1.e11', 'aer_bc_oc'),  # accumulation
#         IcBcMappingItem('ecj -> 0.414*[BC_ks]+0.414*[BC_as]+0.414*[BC_cs]+0.414*[BC_ki];1.e9', 'aer_bc_oc'),  # molar mass 12.01 =>
#
#         # Organics, POA. Check EMAC in nml/oracle.nml.
#         # TODO: this rule is slow due to large number of variables. USE nco to preprocess and sum up all POA. Same for SOA by bins
#         IcBcMappingItem('ac0 -> 6.41*[fPOA01_ks]+6.41*[fPOA02_ks]+6.41*[fPOA03_ks]+6.41*[fPOA04_ks]+6.41*[fPOA05_ks]+'
#                         '6.41*[fPOA01_as]+6.41*[fPOA02_as]+6.41*[fPOA03_as]+6.41*[fPOA04_as]+6.41*[fPOA05_as]+'
#                         '6.41*[bbPOA01_ks]+6.41*[bbPOA02_ks]+6.41*[bbPOA03_ks]+6.41*[bbPOA04_ks]+'
#                         '6.41*[bbPOA01_as]+6.41*[bbPOA02_as]+6.41*[bbPOA03_as]+6.41*[bbPOA04_as];1.e11', 'aer_OA'),
#         IcBcMappingItem('orgpaj -> 8.63*[fPOA01_ks]+8.63*[fPOA02_ks]+8.63*[fPOA03_ks]+8.63*[fPOA04_ks]+8.63*[fPOA05_ks]+'
#                         '8.63*[fPOA01_as]+8.63*[fPOA02_as]+8.63*[fPOA03_as]+8.63*[fPOA04_as]+8.63*[fPOA05_as]+'
#                         '8.63*[bbPOA01_ks]+8.63*[bbPOA02_ks]+8.63*[bbPOA03_ks]+8.63*[bbPOA04_ks]+'
#                         '8.63*[bbPOA01_as]+8.63*[bbPOA02_as]+8.63*[bbPOA03_as]+8.63*[bbPOA04_as];1.e9', 'aer_OA'),  # molar mass 250 =>
#         # Organics, SOA, have to assign according to VBS bins (saturation vapor pressure).
#         # EMAC vars: *sv semivolatile, *iv intermidiate, *v volatile
#         # WRF has 4 bins from 1-10**3 ug/m^3 @ 300K. EMAC bins are in oracle.nml
#         # SOA, VBS bin 1, asoa1j. An example only
#         # IcBcMappingItem('ac0 -> 6.41*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e11', 'aer_OA'),
#         # IcBcMappingItem('orgpaj -> 8.63*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e9', 'aer_OA'),  # molar mass 250 =>
#
#         # Nitrate
#         IcBcMappingItem('nu0 -> 3.56*[NO3m_ns]+3.56*[NO3m_ks];1.e14', 'aer_nitr'),
#         IcBcMappingItem('no3ai -> 2.14*[NO3m_ns]+2.14*[NO3m_ks];1.e9', 'aer_nitr'),  # 62.01/28.97
#         IcBcMappingItem('ac0 -> 4.24*[NO3m_as]+4.24*[NO3m_cs];1.e11', 'aer_nitr'),
#         IcBcMappingItem('no3aj -> 2.14*[NO3m_as]+2.14*[NO3m_cs];1.e9', 'aer_nitr'),
#         # Ammonium
#         IcBcMappingItem('nu0 -> 3.11*[NH4p_ns]+3.11*[NH4p_ks];1.e14', 'aer_nitr'),
#         IcBcMappingItem('nh4ai -> 0.62*[NH4p_ns]+0.62*[NH4p_ks];1.e9', 'aer_nitr'),  # 18.05/28.97
#         IcBcMappingItem('ac0 -> 3.70*[NH4p_as]+3.70*[NH4p_cs];1.e11', 'aer_nitr'),
#         IcBcMappingItem('nh4aj -> 0.62*[NH4p_as]+0.62*[NH4p_cs];1.e9', 'aer_nitr'),
#         ]
#     return mappings


def get_merra2wrf_mapping():
    '''
    chem 100/106 case (MADE SORGAM)

    To prescribe aerosol in the MADE scheme we need to specify two parameters:
    mass concentration (linked to 3rd moment) and number concentration (0th moment).
    3rd moment (mass) can be derived from MERRA2. For 0th moment (number concentration) I assumed fixed size.

    See aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient

    Assume that MERRA2 provides
    '''


    gas_stream = 'inst3_3d_chm_Nv'
    aerosol_stream = 'inst3_3d_aer_Nv'
    mappings = [
        # chem 100: Gases
        IcBcMappingItem('so2 -> 0.453*[SO2];1.e6', aerosol_stream),
        IcBcMappingItem('o3 -> 0.604*[O3];1.e6', gas_stream),  # O3:units = "kg kg-1" ;
        IcBcMappingItem('co -> 0.966*[CO];1.e6', gas_stream),  # TODO: CO:units = "mol mol-1" ;
        # MERRA2 does not have CH4

        # chem 100: Aerosols
        # TODO: BC & SO4 should likely be assigned to the i-mode. Check emissions (speciation) break down.
        # TODO: Here currently everything goes to j mode.
        # put all SO4 into accumulation mode
        IcBcMappingItem('ac0 -> 3.56*[SO4];1.e11', aerosol_stream),  # MERRA2_400.inst3_3d_aer_Nv.20*
        IcBcMappingItem('so4aj -> 1.0*[SO4];1.e9', aerosol_stream),
        # put all dust into the coarse mode
        IcBcMappingItem('corn -> 1.68*[DU001]+1.68*[DU002]+1.68*[DU003]+1.68*[DU004]+1.68*[DU005];1.e7', aerosol_stream),
        IcBcMappingItem('soila -> 1.0*[DU001]+1.0*[DU002]+1.0*[DU003]+1.0*[DU004]+1.0*[DU005];1.e9', aerosol_stream),
        # put all sea salt into the coarse mode.
        # TODO: Check what should be the SD, for now assume the same
        IcBcMappingItem('corn -> 1.98*[SS001]+1.98*[SS002]+1.98*[SS003]+1.98*[SS004]+1.98*[SS005];1.e7', aerosol_stream),
        IcBcMappingItem('seas -> 1.0*[SS001]+1.0*[SS002]+1.0*[SS003]+1.0*[SS004]+1.0*[SS005];1.e9', aerosol_stream),
        # put all organic aerosols into acc mode
        # put BC into elemental carbon, consistent with emissions
        IcBcMappingItem('ac0 -> 6.41*[BCPHOBIC]+6.41*[BCPHILIC];1.e11', aerosol_stream),
        IcBcMappingItem('ecj -> 1.0*[BCPHOBIC]+1.0*[BCPHILIC];1.e9', aerosol_stream),
        # assign OC into primary orgpaj in WRF, similar to the emissions
        IcBcMappingItem('ac0 -> 6.41*[OCPHOBIC]+6.41*[OCPHILIC];1.e11', aerosol_stream),
        IcBcMappingItem('orgpaj -> 1.0*[OCPHOBIC]+1.0*[OCPHILIC];1.e9', aerosol_stream),
    ]
    return mappings


# old implementation for varios chem opts preserved for reference purposes
# def get_merra2wrf_config():
#     config = Config()
#
#     config.root_path = '/project/k1090/osipovs'  # SHAHEEN
#     config.root_path = '/work/mm0062/b302074'  # MISTRAL
#
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/"
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2017-09-01"
#     config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2017-10-31"
#     # config.wrf_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/2018-03-08"
#     config.wrf_input_file = "wrfinput_d01"
#     config.wrf_bdy_file = "wrfbdy_d01"
#
#     config.wrf_met_dir = config.root_path + "/Data/AirQuality/AQABA/IC_BC/met_em/"
#     config.wrf_met_files = "met_em.d01.2017-0*"
#
#     config.mera_dir = config.root_path + "/Data/NASA/MERRA2/inst3_3d_aer_Nv/"
#     config.mera_files = "MERRA2_400.inst3_3d_aer_Nv.201*"  # "MERRA2_400.inst3_3d_aer_Nv.20170*"
#
#     config.do_IC = False  # True
#     config.do_BC = True
#
#     # GOCART DUST ONLY
#     config.spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
#                'DUST_2 -> 1.0*[DU002];1.e9',
#                'DUST_3 -> 1.0*[DU003];1.e9',
#                'DUST_4 -> 1.0*[DU004];1.e9',
#                'DUST_5 -> 1.0*[DU005];1.e9']
#
#     # GOCART FULL
#     config.spc_map = ['DUST_1 -> 1.0*[DU001];1.e9',
#                'DUST_2 -> 1.0*[DU002];1.e9',
#                'DUST_3 -> 1.0*[DU003];1.e9',
#                'DUST_4 -> 1.0*[DU004];1.e9',
#                'DUST_5 -> 1.0*[DU005];1.e9',
#                'SEAS_1 -> 1.0*[SS002];1.e9',
#                'SEAS_2 -> 1.0*[SS003];1.e9',
#                'SEAS_3 -> 1.0*[SS004];1.e9',
#                'SEAS_4 -> 1.0*[SS005];1.e9',
#                'so2 -> 0.453*[SO2];1.e6',
#                'sulf -> 0.302*[SO4];1.e6',
#                'BC1 -> 0.4143*[BCPHOBIC];1.e9', 'BC2 -> 0.4143*[BCPHILIC];1.e9',
#                'OC1 -> 0.4143*[OCPHOBIC];1.e9', 'OC2 -> 0.4143*[OCPHILIC];1.e9',
#                'dms -> 0.467*[DMS];1.e6']
#     # ,'msa -> 0.302*[MSA];1.e6'
#
#     # config.spc_map = ['so2 -> 0.453*[SO2];1.e6', 'sulf -> 0.302*[SO4];1.e6']
#     # config.spc_map = ['o3 -> 0.604*[O3];1.e6', 'co -> 1*[CO];1.e6']
#
#     '''
#     #CBMZ-MOSAIC_8bins  DUST only
#     config.spc_map =['oin_a02->4.92e-4*[DU001];1.e9',
#               'oin_a03->8.94e-3*[DU001];1.e9',
#               'oin_a04->0.14300*[DU001];1.e9',
#               'oin_a05->0.84740*[DU001]+0.1520*[DU002];1.e9',
#               'oin_a06->0.84800*[DU002]+0.4055*[DU003];1.e9',
#               'oin_a07->0.59450*[DU003]+0.4480*[DU004];1.e9',
#               'oin_a08->0.55200*[DU004]+1*[DU005];1.e9',
#
#               'num_a02->0.213e14*[DU001];1',
#               'num_a03->0.598e14*[DU001];1',
#               'num_a04->1.192e14*[DU001];1',
#               'num_a05->1.434e14*[DU001]+0.948e13*[DU002];1',
#               'num_a06->2.086e13*[DU002]+0.358e13*[DU003];1',
#               'num_a07->0.255e13*[DU003]+0.592e12*[DU004];1',
#               'num_a08->0.296e12*[DU004]+1.655e11*[DU005];1']
#     '''
#
#     return config

class IcBcMappingItem():
    def __init__(self, mapping_rule_str, output_stream):
        self.mapping_rule_str = mapping_rule_str
        self.output_stream = output_stream  # plugs into the filename template
