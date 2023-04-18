import numpy as np
from emac2wrf.bc_ic_utils import IcBcMappingItem


def get_emac2wrf_mapping():

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
        IcBcMappingItem('so2 -> 1*[SO2];1.e6', 'tr_sulfur'),  # SO2_ave
        IcBcMappingItem('ch4 -> 1*[CH4];1.e6', 'tr_hycarbs'),  # _ave
        IcBcMappingItem('co -> 1*[CO];1.e6', 'tr_alks'),
        IcBcMappingItem('o3 -> 1*[O3];1.e6', 'tr_Ox_HOx'),
        # TODO: I can map H2SO4 too, but the life time is short compared to SO4, I can ignore it
        # VOCs, but currently only longer lived ones
        IcBcMappingItem('hono -> 1*[HONO];1.e6', 'tr_NOx_NOy'),
        IcBcMappingItem('hcho -> 1*[HCHO];1.e6', 'tr_alks'),
        IcBcMappingItem('eth -> 1*[C2H6];1.e6', 'tr_hycarbs'),  # ethane
        IcBcMappingItem('hc3 -> 1*[C3H8];1.e6', 'tr_hycarbs'),  # propane in Volatility Basis Set

        # aerosols, see aerosol_parameters_for_MADE_IC_BC.py on details how to compute the mapping coefficient
        # aerosols in EMAC are in ppv, in WRF in mmr.
        # Compute mmr = vmr * Ma/Mb, where Mb is dry air molar mass = 28.97 g/mol
        # Check EMAC nc for molar masses & densities

        # Sulfate, split into nu + acc
        # TODO: what is mm in SO4mm_?
        IcBcMappingItem('nu0 -> 2.99*[SO4mm_ns]+2.99*[SO4mm_ks];1.e14', 'aer_sulf'),
        IcBcMappingItem('so4ai -> 3.315*[SO4mm_ns]+3.315*[SO4mm_ks];1.e9', 'aer_sulf'),    # 96.06/28.97 = 3.315
        IcBcMappingItem('ac0 -> 3.56*[SO4mm_as]+3.56*[SO4mm_cs];1.e11', 'aer_sulf'),
        IcBcMappingItem('so4aj -> 3.315*[SO4mm_as]+3.315*[SO4mm_cs];1.e9', 'aer_sulf'),
        # Dust, only coarse
        IcBcMappingItem('corn -> 1.65*[DU_as]+1.65*[DU_ai]+1.65*[DU_cs]+1.65*[DU_ci];1.e7', 'aer_du_ss'),
        IcBcMappingItem('soila -> 1.38*[DU_as]+1.38*[DU_ai]+1.38*[DU_cs]+1.38*[DU_ci];1.e9', 'aer_du_ss'),  # 40.08/28.97=1.38
        # Sea salt
        IcBcMappingItem('corn -> 2.01*[SS_ks]+2.01*[SS_as]+2.01[SS_cs];1.e7', 'aer_du_ss'),  # coarse
        IcBcMappingItem('seas -> 2.017*[SS_ks]+2.017*[SS_as]+2.017*[SS_cs];1.e9', 'aer_du_ss'),  # SS molar mass is 58.44; 58.44/28.97=2.017
        # TODO: Cations: Mg, Ca and such (from Dust and Sea Salt)
        # Organics, BC. BC probably should be in the nucleation mode, not acc
        IcBcMappingItem('ac0 -> 3.2*[BC_ks]+3.2*[BC_as]+3.2*[BC_cs]+3.2*[BC_ki];1.e11', 'aer_bc_oc'),  # accumulation
        IcBcMappingItem('ecj -> 0.414*[BC_ks]+0.414*[BC_as]+0.414*[BC_cs]+0.414*[BC_ki];1.e9', 'aer_bc_oc'),  # molar mass 12.01 =>

        # Organics, POA. Check EMAC in nml/oracle.nml.
        # TODO: this rule is slow due to large number of variables. USE nco to preprocess and sum up all POA. Same for SOA by bins
        IcBcMappingItem('ac0 -> 6.41*[fPOA01_ks]+6.41*[fPOA02_ks]+6.41*[fPOA03_ks]+6.41*[fPOA04_ks]+6.41*[fPOA05_ks]+'
                        '6.41*[fPOA01_as]+6.41*[fPOA02_as]+6.41*[fPOA03_as]+6.41*[fPOA04_as]+6.41*[fPOA05_as]+'
                        '6.41*[bbPOA01_ks]+6.41*[bbPOA02_ks]+6.41*[bbPOA03_ks]+6.41*[bbPOA04_ks]+'
                        '6.41*[bbPOA01_as]+6.41*[bbPOA02_as]+6.41*[bbPOA03_as]+6.41*[bbPOA04_as];1.e11', 'aer_OA'),
        IcBcMappingItem('orgpaj -> 8.63*[fPOA01_ks]+8.63*[fPOA02_ks]+8.63*[fPOA03_ks]+8.63*[fPOA04_ks]+8.63*[fPOA05_ks]+'
                        '8.63*[fPOA01_as]+8.63*[fPOA02_as]+8.63*[fPOA03_as]+8.63*[fPOA04_as]+8.63*[fPOA05_as]+'
                        '8.63*[bbPOA01_ks]+8.63*[bbPOA02_ks]+8.63*[bbPOA03_ks]+8.63*[bbPOA04_ks]+'
                        '8.63*[bbPOA01_as]+8.63*[bbPOA02_as]+8.63*[bbPOA03_as]+8.63*[bbPOA04_as];1.e9', 'aer_OA'),  # molar mass 250 =>
        # Organics, SOA, have to assign according to VBS bins (saturation vapor pressure).
        # EMAC vars: *sv semivolatile, *iv intermidiate, *v volatile
        # WRF has 4 bins from 1-10**3 ug/m^3 @ 300K. EMAC bins are in oracle.nml
        # SOA, VBS bin 1, asoa1j. An example only
        # IcBcMappingItem('ac0 -> 6.41*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e11', 'aer_OA'),
        # IcBcMappingItem('orgpaj -> 8.63*[fSOAsv01_ks]+6.41*[bbSOAsv01_ks]+6.41*[fSOAiv01_ks]+6.41*[bbSOAiv01_ks]+6.41*[SOAv0101_ks];1.e9', 'aer_OA'),  # molar mass 250 =>

        # Nitrate
        IcBcMappingItem('nu0 -> 3.56*[NO3m_ns]+3.56*[NO3m_ks];1.e14', 'aer_nitr'),
        IcBcMappingItem('no3ai -> 2.14*[NO3m_ns]+2.14*[NO3m_ks];1.e9', 'aer_nitr'),  # 62.01/28.97
        IcBcMappingItem('ac0 -> 4.24*[NO3m_as]+4.24*[NO3m_cs];1.e11', 'aer_nitr'),
        IcBcMappingItem('no3aj -> 2.14*[NO3m_as]+2.14*[NO3m_cs];1.e9', 'aer_nitr'),
        # Ammonium
        IcBcMappingItem('nu0 -> 3.11*[NH4p_ns]+3.11*[NH4p_ks];1.e14', 'aer_nitr'),
        IcBcMappingItem('nh4ai -> 0.62*[NH4p_ns]+0.62*[NH4p_ks];1.e9', 'aer_nitr'),  # 18.05/28.97
        IcBcMappingItem('ac0 -> 3.70*[NH4p_as]+3.70*[NH4p_cs];1.e11', 'aer_nitr'),
        IcBcMappingItem('nh4aj -> 0.62*[NH4p_as]+0.62*[NH4p_cs];1.e9', 'aer_nitr'),
        ]
    return mappings

