# Needs pylab
import numpy as np
from dataclasses import dataclass, field
# import pandas as pd
from components import Resistor, Capacitor, Inductor, Mosfet, GateDriver


# ***********************************************************************************************
#                    POWER CONVERTERS
# ***********************************************************************************************
@dataclass
class Requirements:

    vin_min:float = 0.0
    vin_max:float = 0.0
    eff:float = 1.0
    fsw:float = 0
    fc:float = 0
    vin_ripple:float = 0
    vout_ss_ripple:float = 0
    vout_ac_ripple:float = 0
    il_rfactor:float = 0



# TODO: create factory class that generates DCDC converter based on requirements

# ***********************************************************************************************
#                    SWITCHMODE DCDC  POWER TOPOLOGIES
# ***********************************************************************************************

class SyncBuck:
    """Synchronous Buck DCDC Converter Class"""

    
    def __init__(self,  vin:float,
                        vout:float,
                        pout:float,
                        requirements:Requirements,
                        inductor:Inductor=Inductor(),
                        outcap:Capacitor=Capacitor(),
                        incap:Capacitor=Capacitor(),
                        active_sw:Mosfet=Mosfet(),
                        passive_sw:Mosfet=Mosfet(),
                        hs_gatedrive:GateDriver=GateDriver(),
                        ls_gatedrive:GateDriver=GateDriver(),
                        n_active_sw:int = 1,
                        n_passive_sw:int = 1,
                        Nphases:int = 1):


        self.Lout = inductor
        self.Cin = incap
        self.Cout = outcap
        self.Qhs = active_sw
        self.Qls= passive_sw 
        self.Gdhs = hs_gatedrive
        self.Gdls = ls_gatedrive
        self.nQhs = n_active_sw
        self.nQls = n_passive_sw
        self.nPhases = Nphases

        #Connect Gate drivers
        #TODO: change gatedrive class to connect to a FET instead
        self.Qhs.connect_gatedrive(self.Gdhs)
        self.Qls.connect_gatedrive(self.Gdls)

        self.req = requirements
        self.eff = self.req.eff if self.req.eff > 0.0 else 1.0
        self.vin = vin if vin > 0.0 else 0.0
        self.vout = vout if vout > 0.0 else 0.0
        self.pout = pout if pout > 0.0 else 0.0





        self._calc_losses = False

    @property
    def pin(self):
        return self.pout/self.eff

    @property
    def iout(self):
        return self.pout/self.vout

    @property
    def iin(self):
        return self.pin/self.vin


    @property
    def dc(self):
        return self.vout/(self.vin + 1e-22) / self.eff

    @property
    def Icin_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.nPhases)
        n = self.nPhases
        return np.sqrt((D-m/n)*((1+m)/n - D))

    @property
    def Icout_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.nPhases)
        n = self.nPhases
        return n/(D*(1-D))*((D-m/n)*((1+m)/n - D))


    def update(self, calc_losses = True):

        
        # ------- Calculate Phase Current ------------ 
        self.iout_phase = self.iout / self.nPhases

        # ------- Calculate Inductor Current ------------ 

        self.Lout.Iavg = self.iout_phase
        self.Lout.von = self.vin - self.vout
        self.Lout.vsec = self.Lout.von * self.dc / self.req.fsw
        self.Lout.vusec = self.Lout.vsec / 1e-6

        if self.Lout.L == 0.0:
            self.Lout.Ipkpk = self.Lout.Iavg * self.req.il_rfactor
            self.Lout.L = self.Lout.vsec / self.Lout.Ipkpk
        else:
            self.Lout.Ipkpk = self.Lout.vsec / self.Lout.L

        self.Lout.Ipk = self.Lout.Iavg + self.Lout.Ipkpk/2
        self.Lout.Ivalley = self.Lout.Iavg - self.Lout.Ipkpk/2
        self.Lout.Irms = np.sqrt((self.Lout.Ipk**2 + self.Lout.Ivalley**2 + self.Lout.Ipk*self.Lout.Ivalley)/3)

        # ------- Calculate Output Capacitance ------------ 
        self.Cout.Q =  self.Lout.Ipkpk / (8 * self.req.fsw)
        self.Cout.Css = self.Cout.Q / self.req.vout_ss_ripple
        load_step = self.iout * 0.5
        self.Cout.Cac = load_step / (2*np.pi*self.req.fc*self.req.vout_ac_ripple)
        self.Cout.Irms = self.Icout_rms_ratio * self.iout
        
        self.Cout.C = max(self.Cout.Cac,self.Cout.Css) * self.nPhases

        # ------- Calculate FET Currents ------------ 

        # Active FET
        self.Qhs.vds = self.vin
        self.Qhs.Idpk = self.Lout.Ipk/self.nQhs
        self.Qhs.Idrms = self.Lout.Iavg*np.sqrt(self.dc)/self.nQhs
        self.Qhs.Idpk_total = self.Qhs.Idpk * self.nQhs
        self.Qhs.Idrms_total = self.Qhs.Idrms * self.nQhs
        
        # Sync|Passive FET
        self.Qls.vds = self.vin
        self.Qls.Idpk = self.Lout.Ipk/self.nQls
        self.Qls.Idrms = self.Lout.Iavg*np.sqrt(1-self.dc)/self.nQls
        self.Qls.Idpk_total = self.Qls.Idpk * self.nQls
        self.Qls.Idrms_total = self.Qls.Idrms * self.nQls

        # ------- Calculate Input Capacitance ------------ 

        self.Cin.Irms = self.Icin_rms_ratio * self.iout
        self.Cin.Q = self.iin * (1-self.dc)/self.req.fsw
        self.Cin.C = self.Cin.Q / self.req.vin_ripple

        if calc_losses == True:
            self._calc_losses = True

            # --------- Loss Calculation ------------ 
            # Gate Driver Loss
            self.Gdhs.Ploss = self.Qhs.qtot * self.nQhs * self.Gdhs.vdrv * self.req.fsw
            self.Gdls.Ploss = self.Qls.qtot * self.nQls * self.Gdls.vdrv * self.req.fsw

            # Active Switch Loss (per Device)
            self.Qhs.get_switch_times()
            self.Qhs.Ploss['ohm'] = self.Qhs.rdson * self.Qhs.Idrms**2
            self.Qhs.Ploss['sw'] = (self.Lout.Iavg/self.nQhs)*self.Qhs.vds*(self.Qhs.tr + self.Qhs.tf)*self.op.fsw
            self.Qhs.Ploss['dt'] = self.Qhs.vfwd * self.Qhs.Idrms * self.Gdhs.deadtime * 2 * self.req.fsw
            self.Qhs.Ploss['coss'] = 1.0 * self.Qhs.coss * self.Qhs.vds * self.req.fsw

            # Sync. Rectifier Loss (per Device)
            self.Qls.Ploss['ohm'] = self.Qls.rdson * self.Qls.Idrms**2
            self.Qls.Ploss['dt'] = self.Qls.vfwd * self.Qls.Idrms * self.Gdls.deadtime * 2 * self.op.fsw

            # Inductor Loss
            self.Lout.Ploss['ohm'] = self.Lout.dcr * self.Lout.Irms**2

    def __repr__(self):
        string = "Sync Buck:"


