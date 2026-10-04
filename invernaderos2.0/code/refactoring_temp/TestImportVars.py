#%%

from sympy import symbols

mt, mg, m, C, s, W, mg_CO2, J, g, mol_CH2O = symbols('mt mg m C s W mg_CO2 J g mol_CH2O')

mt, mg, m, C, s, W, mg_CO2, J, Pa, kg_water, kg, K, ppm, m_cover, kg_air = symbols('mt mg m C s W mg_CO2 J Pa kg_water kg K ppm m_cover kg_air')  # Symbolic use of base phisical units

mt, mg, m, C, s, W, mg_CO2, J, Pa, kg_water, kg, K, ppm, kmol, kg_air, kg_vapour, mxn, eur = symbols('mt mg m C s W mg_CO2 J Pa kg_water kg K ppm kmol kg_air kg_vapour mxn eur')  # Symbolic use of base phisical units

kW, hour = symbols('kW hour')
truss,fruit, plant = symbols('truss,fruit plant')
mm = symbols('mm')
m = symbols('m')
Euros = symbols('Euros')
ok = 'OK'

"""
from .parameters_dir import PARAMS_DIR
from .parameters_env import day2seconds
from .parameters_dt import DT
from .parameters_ddpg import CONTROLS as config_controls
"""

days = "temp_days" #PARAMS_DIR['days']
dt = 'ModuleClimate' #DT['ModuleClimate']
nrec = 60*60*24 #int(day2seconds(days)/dt)
mt = symbols('mt') #Minutos
MODEL_NOISE = False
mt = symbols('mt') #Minutos


def Struct( typ, varid, prn, desc, rec=nrec, ok=ok, units=1, val=1):
    return { 'typ':typ, 'varid':varid, 'prn':prn, 'desc':desc, 'ok':ok, 'units':units, 'val':val}


################## Constants ##################
OTHER_CONSTANTS = {     
    ################## other constants ################## 
    'etagas':      Struct(typ='Cnts', varid='etagas', prn=r'$\eta_{gas}$',
                    desc="Energy efficiency of natural gas", units= kW * hour * m**-3 , val=9.794, ok='Data hacakaton'),  
    'qh2o':      Struct(typ='Cnts', varid='qh2o', prn=r'$\q_{H_{2}O}$',
                    desc="Costo del agua", units= mxn *m**-3, val=3.5, ok='Data hacakaton'), 
    'qgas':    Struct(typ='Cnts', varid='qgas', prn=r'$q_{gas}$',
                    desc="Cost of natural gas", units= eur * kW**-1 * hour**-1, val=0.09, ok='Data hacakaton'),      
    'q_co2_ext': Struct(typ='Cnts', varid='q_co2_ext', prn=r'$\q_{CO_2}_{ext}$',
                    desc="", units=eur * kg**-1, val=0.3, ok=ok), #val=0.175    # Costo del gas de la fuente externa lo tomamos al precio de la tesis 
    'cost_elect': Struct(typ='Cnts', varid='cost_elect', prn=r'$\cost_{elect}$',
                    desc="Cost of the electricity", units=eur * kW**-1 * hour**-1, val=0.3, ok='Data hacakaton'),    #val=0.68                     
    'T_cal':     Struct(typ='Cnts', varid='T_cal', prn=r'$T_{cal}$',
                    desc="Missing", units=1, val=65, ok='falta descripción y unidades'),          # Temperatura máxima de la caldera  
    'sigma':     Struct(typ='Cnts', varid='sigma', prn=r'$\sigma$',
                    desc="Stefan-Boltzmann constant", units=W * m**-2 * K**-4, val=5.670e-8, ok=ok), # Constante de Stefan-Boltzmann (W m−2)
    'etadrain':  Struct(typ='Cnts', varid='etadrain', prn=r'$\eta_{drain}$',
                    desc="Missing", units=1, val=30, ok='falta descripción y unidades'),
    'HEAT_PIPE' :Struct(typ='Cnts', varid='HEAT_PIPE', prn=r'$HEAT_PIPE$',
                    desc="Temperatura del tubo de calentamiento", units=1, val=90, ok=ok),
    'n_pipes':Struct(typ='Cnts', varid='n_pipes', prn=r'$N_{pipes}$',
                    desc="Numero de tuberias de calentamiento", units=1, val=1, ok='Multiplica h_4 y r_6'),
    'RH':Struct( typ='State', varid='RH', prn=r'$RH$',\
           desc="Relative humidity percentage in the greenhouse air", \
           units=1,rec=nrec, val=50),
    'VPD':Struct( typ='State', varid='VPD', prn=r'$VPD$',\
           desc="Presion de vapor de saturacion", \
           units=1,rec=nrec, val=70),
    'lambda':Struct(typ='Cnts', varid='lambda', prn=r'$\lambda$',\
           desc="Constante de la ecuacion para U7 ", \
           units=1,rec=nrec, val=0.0001),
    'U7_c':Struct(typ='State', varid='U7_c', prn=r'$U_{7c}',\
           desc="Punto final para el control U7 ", \
           units=1,rec=nrec, val=0.1),
    'U6_c':Struct(typ='State', varid='U6_c', prn=r'$U_{6c}',\
           desc="Punto final para el control U6 ", \
           units=1,rec=nrec, val=0.1),
    'LampEndtime':Struct(typ='State', varid='LampEndtime', prn=r'LampEndtime',\
           desc="Hora de apagado de luces ", \
           units=1, rec=nrec, val=20 ), # Unidades horas
    'hour':Struct(typ='State', varid='hour', prn=r'hour',\
           desc="Hora del dia", \
           units=1,rec=nrec, val=0),
    'dia':Struct(typ='State', varid='dia', prn=r'dia',\
           desc="Fecha", \
           units=1,rec=nrec, val=0),
    'Potential_value':Struct(typ='State', varid='Potential_value', prn=r'Potential_value',\
           desc="Value in Euros per sq meter of potted tomato plants currently in greenhouse", \
           units= Euros / m**2,rec=nrec, val=0),
    'PARsum':Struct(typ='State', varid='PARsum', prn=r'PARSum',\
           desc="Suma del PAR que recibe el canopy por dia (moles * m-2 * dia-1)", \
           units=1,rec=nrec, val=0),
    'month':Struct(typ='State', varid='month', prn=r'month',\
           desc="Fecha", \
           units=1,rec=nrec, val=0),
    'year':Struct(typ='State', varid='year', prn=r'year',\
           desc="Fecha", \
           units=1,rec=nrec, val=0)
}


from numpy import zeros


class test:
    def __init__( self, a, b, Vars):
        self.a = a 
        self.b = b
        self.Vars = Vars
        for key in self.Vars.keys():    
            self.__dict__[key] = self.Vars[key]['val']
    
    def Calculo(self):
        ### Use the variables in a regular fashion
        return self.HEAT_PIPE * self.sigma
    
    def Info( self, key):
        ### Update the value
        self.Vars[key]['val'] = self.__dict__[key]
        return self.Vars[key] #return the info


if __name__ == "__main__":
    
    ts = test( a=2, b=3, Vars=OTHER_CONSTANTS)
    print(ts.Calculo())
    print(ts.n_pipes)
    ts.U7_c = 30
    print(ts.Info('U7_c'))

        
        
        
        
# %%
