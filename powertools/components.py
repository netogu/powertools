from dataclasses import dataclass, field
import numpy as np
# import pandas as pd


# ***********************************************************************************************
#                    PASSIVE DEVICES
# ***********************************************************************************************
class Passive:
    #Operating Point
    volt:float = 0.0
    irms:float = 0.0
    ploss:float = 0.0

    #Electrical Characteristics
    R:float = 0.0
    C:float = 0.0

    def _calc_std_value(self,value, eseries):
        frac,decade = np.modf(np.log10(value))
        return round(10**(round(frac*eseries)/eseries),2)*10**(decade)


class Resistor(Passive):
    name:str = 'R'

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

class Capacitor(Passive):
    name:str = 'C'
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
        

class Inductor:
    name:str = 'L'
    #Operating Point
    von = 0
    voff = 0
    vsec:float = 0.0
    vusec:float = 0.0
    iavg = 0.0
    ipkpk = 0.0
    ivalley = 0.0
    irms = 0.0
    ploss = 0.0

    #Electrical Characteristics
    value:float = 0.0
    dcr = 0.0
    ilim_sat = 0.0
    irms_20C = 0.0
    irms_40C = 0.0
    
    #Mechanical Characteristics
    xdim:float = 0.0
    ydim:float = 0.0
    weight:float = 0.0

    #Misc
    cost:float = 0.0

    @property
    def irf(self):
        return self.Ipkpk/self.Iavg

    def __repr__(self):

        s  = f"L = {self.value/1e-6:2.2f}uH\n"
        s += f"iLrms = {self.irms:2.2f}Arms\n"
        s += f"iLpeak = {self.ipk:2.2f}A\n"
        s += f"iLvalley = {self.ivalley:2.2f}A\n"
        return s

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

class GateDriver:

    name:str = 'Gate Driver'
    #Operating Point
    vdrv:float = 0.0
    igdrv:float = 0.0
    ploss:float = 0.0
    #Electrical Characteristics
    rpull_up:float = 0.0
    rpull_down:float = 0.0
    ipk_source:float = 0.0
    ipk_sink:float = 0.0
    deadtime:float = 0.0
    #Mechanical
    xdim:float = 0.0
    ydim:float = 0.0
    weight:float = 0.0
    #Misc
    cost:float = 0.0

class Mosfet:
    
    name:str = 'Ideal FET'
    #Operating Point
    vds:float = 0.0
    idrms:float = 0.0
    idpk:float = 0.0
    tamb:float = 0.0
    tj:float = 0.0
    ploss:dict = 0.0

    @property
    def ploss_total(self):
        return sum(self.ploss.values())

    def connect_gatedrive(self,gatedrive:GateDriver):
        self.gatedrv = gatedrive

    def get_switch_times(self):
        self.igon  = (self.gatedrv.vdrv - self.vmiller) / (self.rg + self.gatedrv.rpull_up)
        self.igoff = abs((0.0 - self.vmiller) / (self.rg + self.gatedrv.rpull_down))

        if self.igon > self.gatedrv.ipk_source:
            self.igon= self.gatedrv.ipk_source
        if self.igoff > self.gatedrv.ipk_sink:
            self.igoff= self.gatedrv.ipk_sink
                    
        self.tr = self.qsw/self.igon
        self.tf = self.qsw/self.igoff
        
    #Electrical Characteristics
    rdson:float = 0.0
    coss:float = 0.0
    vfwd:float = 0.0
    vmiller:float = 0.0
    qgs:float = 0.0
    qgd:float = 0.0
    qtot:float = 0.0
    qrr:float = 0.0
    trr:float = 0.0
    rg:float = 0.0
    vgsth:float = 0.0
    
    #Mechanical Characteristics
    xdim:float = 0.0
    ydim:float = 0.0
    weight:float = 0.0
    #Thermal Characteristics
    rthjc:float = 0.0
    #Misc Characteristics
    cost:float = 0.0

    def __repr__(self):

        s  = f"Name = {self.name}\n"
        s += f"Vds = {self.vds:2.2f}V\n"
        s += f"idrms = {self.idrms:2.2f}Arms\n"
        s += f"idpeak = {self.idpk:2.2f}A\n"
        return s

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
