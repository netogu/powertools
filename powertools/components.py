from dataclasses import dataclass, field
import numpy as np
# import pandas as pd


# ***********************************************************************************************
#                    PASSIVE DEVICES
# ***********************************************************************************************

@dataclass
class Passive:
    #Operating Point
    volt:float = field(default=0)
    Irms:float = field(default=0)
    Ploss:float = field(default=0)

    #Electrical Characteristics
    R:float = field(default=0)
    C:float = field(default=0)

    def _calc_std_value(self,value, eseries):
        frac,decade = np.modf(np.log10(value))
        return round(10**(round(frac*eseries)/eseries),2)*10**(decade)


@dataclass
class Resistor(Passive):
    name:str = field(default='R')

    @property
    def R_e96(self):
        """Calculate closest 1% comercial value [E96]"""
        return self._calc_std_value(self.R, eseries=96.0)

    @property
    def R_e48(self):
        """Calculate closest 2% comercial value [E48]"""
        return self._calc_std_value(self.R, 48.0)

    @property
    def R_e24(self):
        """Calculate closest 5% comercial value [E24]"""
        return self._calc_std_value(self.R, 24.0)

@dataclass
class Capacitor(Passive):
    name:str = field(default='C')

    @property
    def C_e24(self):
        """Calculate closest 5% comercial value [E24]"""
        return self._calc_std_value(self.C, 24.0)

    @property
    def C_e12(self):
        """Calculate closest 10% comercial value [E12]"""
        return self._calc_std_value(self.C, 12.0)

    @property
    def C_e6(self):
        """Calculate closest 20% comercial value [E6]"""
        return self._calc_std_value(self.C, 6.0)
        

@dataclass
class Inductor:
    name:str = field(default='Generic')
    #Operating Point
    von:float = field(default=0,repr=False)
    voff:float = field(default=0,repr=False)
    vsec:float = field(default=0,repr=False)
    vusec:float = field(default=0,repr=True)
    Iavg:float = field(default=0,repr=True)
    Ipk:float = field(default=0,repr=False)
    Ipkpk:float = field(default=0,repr=True)
    Ivalley:float = field(default=0,repr=False)
    Irms:float = field(default=0,repr=True)
    Ploss:float = field(default_factory=lambda:{'core':0.0,'ohm':0.0})

    #Electrical Characteristics
    L:float = field(default=0,repr=True)
    dcr:float = field(default=0,repr=False)
    Isat:float = field(default=0,repr=False)
    Irms_20C:float = field(default=0,repr=False)
    Irms_40C:float = field(default=0,repr=False)
    
    #Mechanical Characteristics
    xdim:float = field(default=0,repr=False)
    ydim:float = field(default=0,repr=False)
    weight:float = field(default=0,repr=False)

    #Misc
    cost:float = field(default=0,repr=False)

    @property
    def Irf(self):
        return self.Ipkpk/self.Iavg

    # def load(self, df:pd.DataFrame, name:str):
    #     self.name = name
    #     self.L = df.loc[name]['L']
    #     self.dcr = df.loc[name]['dcr']
    #     self.Isat = df.loc[name]['Isat']
    #     self.Irms_20C = df.loc[name]['Irms_20C']
    #     self.Irms_40C = df.loc[name]['Irms_40C']
    #     self.xdim = df.loc[name]['xdim']
    #     self.ydim = df.loc[name]['ydim']
    #     self.zdim = df.loc[name]['zdim']
    #     self.weight = df.loc[name]['weight']
    #     self.cost = df.loc[name]['cost']


# ***********************************************************************************************
#                    POWER DEVICES 
# ***********************************************************************************************

@dataclass
class GateDriver:

    name:str = field(default='Generic')
    #Operating Point
    vdrv:float = field(default=0.0)
    Igdrv:float = field(default=0.0)
    Ploss:float = field(default=0.0)
    #Electrical Characteristics
    rpull_up:float = field(default=0.0, repr=False)
    rpull_down:float = field(default=0.0, repr=False)
    Ipk_source:float = field(default=0.0, repr=False)
    Ipk_sink:float = field(default=0.0, repr=False)
    deadtime:float = field(default=0.0, repr=False)
    #Mechanical
    xdim:float = field(default=0.0, repr=False)
    ydim:float = field(default=0.0, repr=False)
    weight:float = field(default=0.0, repr=False)
    #Misc
    cost:float = field(default=0.0, repr=False)

@dataclass
class Mosfet:
    
    name:str = field(default='Generic')
    #Operating Point
    vds:float = field(default=0.0,repr=False)
    Idrms:float = field(default=0.0,repr=True)
    Idpk:float = field(default=0.0,repr=True)
    Tamb:float = field(default=25.0,repr=False)
    Tj:float = field(default=25.0,repr=True)
    Ploss:dict = field(default_factory=lambda: {'ohm':0.0,'sw':0.0,'coss':0.0,'dt':0.0,'qrr':0.0})
    def qsw(self):
    @property
    def Ploss_total(self):
        return sum(self.Ploss.values())

    def connect_gatedrive(self,gatedrive:GateDriver):
        self.gatedrv = gatedrive

    def get_switch_times(self):
        self.Igon  = (self.gatedrv.vdrv - self.vmiller) / (self.rg + self.gatedrv.rpull_up)
        self.Igoff = abs((0.0 - self.vmiller) / (self.rg + self.gatedrv.rpull_down))

        if self.Igon > self.gatedrv.Ipk_source:
            self.Igon= self.gatedrv.Ipk_source
        if self.Igoff > self.gatedrv.Ipk_sink:
            self.Igoff= self.gatedrv.Ipk_sink
                    
        self.tr = self.qsw/self.Igon
        self.tf = self.qsw/self.Igoff
        
    #Electrical Characteristics
    rdson:float = field(default=0,repr=False)
    coss:float = field(default=0,repr=False)
    vfwd:float = field(default=0,repr=False)
    vmiller:float = field(default=0,repr=False)
    qgs:float = field(default=0,repr=False)
    qgd:float = field(default=0,repr=False)
    qtot:float = field(default=0,repr=False)
    qrr:float = field(default=0,repr=False)
    trr:float = field(default=0,repr=False)
    rg:float = field(default=0,repr=False)
    vgsth:float = field(default=0,repr=False)
    
    #Mechanical Characteristics
    xdim:float = field(default=0,repr=False)
    ydim:float = field(default=0,repr=False)
    weight:float = field(default=0,repr=False)
    #Thermal Characteristics
    rthjc:float = field(default=0,repr=False)
    #Misc Characteristics
    cost:float = field(default=0,repr=False)

    # def load(self, df:pd.DataFrame, name:str):
    #     self.name = name
    #     self.rdson = df.loc[name]['rdson']
    #     self.coss = df.loc[name]['coss']
    #     self.vfwd = df.loc[name]['coss']
    #     self.vmiller = df.loc[name]['vmiller']
    #     self.qgs = df.loc[name]['qgs']
    #     self.qgd = df.loc[name]['qgd']
    #     self.qtot = df.loc[name]['qtot']
    #     self.qrr = df.loc[name]['qrr']
    #     self.trr = df.loc[name]['trr']
    #     self.rg = df.loc[name]['rg']
    #     self.vgsth = df.loc[name]['vgsth']
    #     self.rthjc = df.loc[name]['rthjc']
    #     self.xdim = df.loc[name]['xdim']
    #     self.ydim = df.loc[name]['ydim']
    #     self.zdim = df.loc[name]['zdim']
    #     self.cost = df.loc[name]['cost']
