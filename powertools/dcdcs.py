import numpy as np
# import pandas as pd
from .components import Resistor, Capacitor, Inductor, Mosfet, GateDriver


# ***********************************************************************************************
#                    POWER CONVERTERS
# ***********************************************************************************************

class PowerParams:

    #Operating conditions
    vin:float = 0.0
    vout:float = 0.0
    pout:float = 0.0
    eff:float = 0.0
    fsw:float = 0.0

    #Requirements
    fc:float = 0.0
    vin_ripple:float = 0.0
    vout_ss_ripple:float = 0.0
    vout_ac_ripple:float = 0.0
    il_rfactor:float = 0.0

    @property
    def pin(self):
        return self.pout/self.eff
    @property
    def iout(self):
        return self.pout/self.vout
    @property
    def iin(self):
        return self.pin/self.vin

    def __repr__(self):
        s = f"Vin = {self.vin:2.2f}V\n"
        s += f"Vout = {self.vout:2.2f}V\n"
        s += f"Pout = {self.pout:2.2f}W\n"
        s += f"Iout = {self.iout:2.2f}A\n"
        s += f"Fsw = {self.fsw/1e3:2.2f}kHz\n"

        return s



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
    def icin_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.n_phases)
        n = self.n_phases
        return np.sqrt((D-m/n)*((1+m)/n - D))

    @property
    def icout_rms_ratio(self):
        D = self.dc
        m = np.floor(self.dc*self.n_phases)
        n = self.n_phases
        return n/(D*(1-D))*((D-m/n)*((1+m)/n - D))


    def update(self, calc_losses = True):

        # ------- Calculate Phase Current ------------ 
        self.op.iout_phase = self.op.iout / self.n_phases

        # ------- Calculate Inductor Current ------------ 

        self.Lout.iavg = self.op.iout_phase
        self.Lout.von = self.op.vin - self.op.vout
        self.Lout.vsec = self.Lout.von * self.dc / self.op.fsw
        self.Lout.vusec = self.Lout.vsec / 1e-6

        if self.Lout.value == 0.0:
            self.Lout.ipkpk = self.Lout.iavg * self.op.il_rfactor
            self.Lout.value = self.Lout.vsec / self.Lout.ipkpk
        else:
            self.Lout.ipkpk = self.Lout.vsec / self.Lout.value

        self.Lout.ipk = self.Lout.iavg + self.Lout.ipkpk/2
        self.Lout.ivalley = self.Lout.iavg - self.Lout.ipkpk/2
        self.Lout.irms = np.sqrt((self.Lout.ipk**2 + self.Lout.ivalley**2 + self.Lout.ipk*self.Lout.ivalley)/3)

        # ------- Calculate Output Capacitance ------------ 
        self.Cout.Q =  self.Lout.ipkpk / (8 * self.op.fsw)
        self.Cout.Css = self.Cout.Q / self.op.vout_ss_ripple
        load_step = self.op.iout * 0.5
        self.Cout.Cac = load_step / (2*np.pi*self.op.fc*self.op.vout_ac_ripple)
        self.Cout.irms = self.icout_rms_ratio * self.op.iout
        
        self.Cout.C = max(self.Cout.Cac,self.Cout.Css) * self.n_phases

        # ------- Calculate FET Currents ------------ 

        # Active FET
        self.SWhs.vds = self.op.vin
        self.SWhs.idpk = self.Lout.ipk/self.n_swhs
        self.SWhs.idrms = self.Lout.iavg*np.sqrt(self.dc)/self.n_swhs
        self.SWhs.idpk_total = self.SWhs.idpk * self.n_swhs
        self.SWhs.idrms_total = self.SWhs.idrms * self.n_swhs
        
        # Sync|Passive FET
        self.SWls.vds = self.op.vin
        self.SWls.idpk = self.Lout.ipk/self.n_swls
        self.SWls.idrms = self.Lout.iavg*np.sqrt(1-self.dc)/self.n_swls
        self.SWls.idpk_total = self.SWls.idpk * self.n_swls
        self.SWls.idrms_total = self.SWls.idrms * self.n_swls

        # ------- Calculate Input Capacitance ------------ 

        self.Cin.irms = self.icin_rms_ratio * self.op.iout
        self.Cin.Q = self.op.iin * (1-self.dc)/self.op.fsw
        self.Cin.C = self.Cin.Q / self.op.vin_ripple

        if calc_losses == True:
            self._calc_losses = True

            # --------- Loss Calculation ------------ 
            # Gate Driver Loss
            self.Gdhs.ploss = self.SWhs.qtot * self.n_swhs * self.Gdhs.vdrv * self.op.fsw
            self.Gdls.ploss = self.SWls.qtot * self.n_swls * self.Gdls.vdrv * self.op.fsw

            # Active Switch Loss (per Device)
            self.SWhs.get_switch_times()
            self.SWhs.ploss['ohm'] = self.SWhs.rdson * self.SWhs.idrms**2
            self.SWhs.ploss['sw'] = (self.Lout.iavg/self.n_swhs)*self.SWhs.vds*(self.SWhs.tr + self.SWhs.tf)*self.op.fsw
            self.SWhs.ploss['dt'] = self.SWhs.vfwd * self.SWhs.idrms * self.Gdhs.deadtime * 2 * self.op.fsw
            self.SWhs.ploss['coss'] = 1.0 * self.SWhs.coss * self.SWhs.vds * self.op.fsw

            # Sync. Rectifier Loss (per Device)
            self.SWls.ploss['ohm'] = self.SWls.rdson * self.SWls.idrms**2
            self.SWls.ploss['dt'] = self.SWls.vfwd * self.SWls.idrms * self.Gdls.deadtime * 2 * self.op.fsw

            # Inductor Loss
            self.Lout.ploss['ohm'] = self.Lout.dcr * self.Lout.irms**2

    def __repr__(self):
        s = "Sync Buck:\n"
        s += str(self.op)
        s += "----------------------------\n"
        s += str(self.Lout)
        s += "----------------------------\n"
        s += "HS FET:\n"
        s += str(self.SWhs)
        s += "----------------------------\n"
        s += "LS FET:\n"
        s += str(self.SWls)
        s += "----------------------------\n"
        s += "Output Capacitance:\n"
        s += f"Cout(ss) = {self.Cout.Css/1e-6:2.2f}uF\n"
        s += f"Cout(ac) = {self.Cout.Cac/1e-6:2.2f}uF\n"
        s += f"Cout:irms = {self.Cout.irms:2.2f}A\n"
        s += "----------------------------\n"
        s += "Input Capacitance:\n"
        s += f"Cin(ss) = {self.Cin.C/1e-6:2.2f}uF\n"
        s += f"Cin:irms = {self.Cin.irms:2.2f}A\n"
        return s

