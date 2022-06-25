# Needs pylab
import numpy as np
from dataclasses import dataclass, field
# import pandas as pd
from .components import Resistor, Capacitor, Inductor, Mosfet, GateDriver


# ***********************************************************************************************
#                    POWER CONVERTERS
# ***********************************************************************************************

@dataclass
class PowerParams:
    #Operating conditions
    vin:float = field(default=0)
    vout:float = field(default=0)
    pout:float = field(default=0)
    eff:float = field(default=1.0)
    fsw:float = field(default=0)
    #Requirements
    fc:float = field(default=0,repr=False)
    vin_ripple:float = field(default=0,repr=False)
    vout_ss_ripple:float = field(default=0,repr=False)
    vout_ac_ripple:float = field(default=0,repr=False)
    Il_rfactor:float = field(default=0,repr=False)

    @property
    def pin(self):
        return self.pout/self.eff
    @property
    def iout(self):
        return self.pout/self.vout
    @property
    def iin(self):
        return self.pin/self.vin



# ***********************************************************************************************
#                    SWITCHMODE DCDC  POWER TOPOLOGIES
# ***********************************************************************************************

class SyncBuck:
    """Synchronous Buck DCDC Converter Class"""

    def __init__(self,  op:PowerParams, 
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


        self.op = op
        self.Lout = inductor
        self.Cin = incap
        self.Cout = outcap
        self.SWhs = active_sw
        self.SWls = passive_sw
        self.Gdls = hs_gatedrive
        self.Gdhs = ls_gatedrive
        self.n_swhs = n_active_sw
        self.n_swls = n_passive_sw
        self.n_phases = Nphases

        #Connect Gate drivers
        self.SWhs.connect_gatedrive(self.Gdhs)
        self.SWls.connect_gatedrive(self.Gdls)

        self._calc_losses = False


    @property
    def dc(self):
        return self.op.vout/(self.op.vin + 1e-22) / self.op.eff

    @property
    def Icin_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.n_phases)
        n = self.n_phases
        return np.sqrt((D-m/n)*((1+m)/n - D))

    @property
    def Icout_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.n_phases)
        n = self.n_phases
        return n/(D*(1-D))*((D-m/n)*((1+m)/n - D))


    def update(self, calc_losses = True):

        
        # ------- Calculate Phase Current ------------ 
        self.op.iout_phase = self.op.iout / self.n_phases

        # ------- Calculate Inductor Current ------------ 

        self.Lout.Iavg = self.op.iout_phase
        self.Lout.von = self.op.vin - self.op.vout
        self.Lout.vsec = self.Lout.von * self.dc / self.op.fsw
        self.Lout.vusec = self.Lout.vsec / 1e-6

        if self.Lout.L == 0.0:
            self.Lout.Ipkpk = self.Lout.Iavg * self.op.Il_rfactor
            self.Lout.L = self.Lout.vsec / self.Lout.Ipkpk
        else:
            self.Lout.Ipkpk = self.Lout.vsec / self.Lout.L

        self.Lout.Ipk = self.Lout.Iavg + self.Lout.Ipkpk/2
        self.Lout.Ivalley = self.Lout.Iavg - self.Lout.Ipkpk/2
        self.Lout.Irms = np.sqrt((self.Lout.Ipk**2 + self.Lout.Ivalley**2 + self.Lout.Ipk*self.Lout.Ivalley)/3)

        # ------- Calculate Output Capacitance ------------ 
        self.Cout.Q =  self.Lout.Ipkpk / (8 * self.op.fsw)
        self.Cout.Css = self.Cout.Q / self.op.vout_ss_ripple
        load_step = self.op.iout * 0.5
        self.Cout.Cac = load_step / (2*np.pi*self.op.fc*self.op.vout_ac_ripple)
        self.Cout.Irms = self.Icout_rms_ratio * self.op.iout
        
        self.Cout.C = max(self.Cout.Cac,self.Cout.Css) * self.n_phases

        # ------- Calculate FET Currents ------------ 

        # Active FET
        self.SWhs.vds = self.op.vin
        self.SWhs.Idpk = self.Lout.Ipk/self.n_swhs
        self.SWhs.Idrms = self.Lout.Iavg*np.sqrt(self.dc)/self.n_swhs
        self.SWhs.Idpk_total = self.SWhs.Idpk * self.n_swhs
        self.SWhs.Idrms_total = self.SWhs.Idrms * self.n_swhs
        
        # Sync|Passive FET
        self.SWls.vds = self.op.vin
        self.SWls.Idpk = self.Lout.Ipk/self.n_swls
        self.SWls.Idrms = self.Lout.Iavg*np.sqrt(1-self.dc)/self.n_swls
        self.SWls.Idpk_total = self.SWls.Idpk * self.n_swls
        self.SWls.Idrms_total = self.SWls.Idrms * self.n_swls

        # ------- Calculate Input Capacitance ------------ 

        self.Cin.Irms = self.Icin_rms_ratio * self.op.iout
        self.Cin.Q = self.op.iin * (1-self.dc)/self.op.fsw
        self.Cin.C = self.Cin.Q / self.op.vin_ripple

        if calc_losses == True:
            self._calc_losses = True

            # --------- Loss Calculation ------------ 
            # Gate Driver Loss
            self.Gdhs.Ploss = self.SWhs.qtot * self.n_swhs * self.Gdhs.vdrv * self.op.fsw
            self.Gdls.Ploss = self.SWls.qtot * self.n_swls * self.Gdls.vdrv * self.op.fsw

            # Active Switch Loss (per Device)
            self.SWhs.get_switch_times()
            self.SWhs.Ploss['ohm'] = self.SWhs.rdson * self.SWhs.Idrms**2
            self.SWhs.Ploss['sw'] = (self.Lout.Iavg/self.n_swhs)*self.SWhs.vds*(self.SWhs.tr + self.SWhs.tf)*self.op.fsw
            self.SWhs.Ploss['dt'] = self.SWhs.vfwd * self.SWhs.Idrms * self.Gdhs.deadtime * 2 * self.op.fsw
            self.SWhs.Ploss['coss'] = 1.0 * self.SWhs.coss * self.SWhs.vds * self.op.fsw

            # Sync. Rectifier Loss (per Device)
            self.SWls.Ploss['ohm'] = self.SWls.rdson * self.SWls.Idrms**2
            self.SWls.Ploss['dt'] = self.SWls.vfwd * self.SWls.Idrms * self.Gdls.deadtime * 2 * self.op.fsw

            # Inductor Loss
            self.Lout.Ploss['ohm'] = self.Lout.dcr * self.Lout.Irms**2

    def __repr__(self):
        if 
        string = "Sync Buck:"


