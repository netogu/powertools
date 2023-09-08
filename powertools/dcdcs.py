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
                        fet_hs:Mosfet=Mosfet(),
                        fet_ls:Mosfet=Mosfet(),
                        gd_hs:GateDriver=GateDriver(),
                        gd_ls:GateDriver=GateDriver(),
                        n_fet_hs:int = 1,
                        n_fet_ls:int = 1,
                        n_phases:int = 1):


        self.op = op
        self.Lout = inductor
        self.Cin = incap
        self.Cout = outcap
        self.FetH = fet_hs
        self.FetL = fet_ls
        self.GdH = gd_hs
        self.GdL = gd_ls
        self.n_fet_hs = n_active_sw
        self.n_fet_ls = n_passive_sw
        self.n_phases = n_phases
        

        #Connect Gate drivers
        self.FetH.connect_gatedrive(self.GdH)
        self.FetL.connect_gatedrive(self.GdL)

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
        self.FetH.vds = self.op.vin
        self.FetH.idpk = self.Lout.ipk/self.n_fet_hs
        self.FetH.idrms = self.Lout.iavg*np.sqrt(self.dc)/self.n_fet_hs
        self.FetH.idpk_total = self.FetH.idpk * self.n_fet_hs
        self.FetH.idrms_total = self.FetH.idrms * self.n_fet_hs
        
        # Sync|Passive FET
        self.FetL.vds = self.op.vin
        self.FetL.idpk = self.Lout.ipk/self.n_fet_ls
        self.FetL.idrms = self.Lout.iavg*np.sqrt(1-self.dc)/self.n_fet_ls
        self.FetL.idpk_total = self.FetL.idpk * self.n_fet_ls
        self.FetL.idrms_total = self.FetL.idrms * self.n_fet_ls

        # ------- Calculate Input Capacitance ------------ 

        self.Cin.irms = self.icin_rms_ratio * self.op.iout
        self.Cin.Q = self.op.iin * (1-self.dc)/self.op.fsw
        self.Cin.C = self.Cin.Q / self.op.vin_ripple

        if calc_losses == True:
            self._calc_losses = True

            # --------- Loss Calculation ------------ 
            # Gate Driver Loss
            self.GdH.ploss = self.FetH.qtot * self.n_fet_hs * self.GdH.vdrv * self.op.fsw
            self.GdL.ploss = self.FetL.qtot * self.n_fet_ls * self.GdL.vdrv * self.op.fsw

            # Active Switch Loss (per Device)
            self.FetH.get_switch_times()
            self.FetH.ploss['ohm'] = self.FetH.rdson * self.FetH.idrms**2
            self.FetH.ploss['sw'] = (self.Lout.iavg/self.n_fet_hs)*self.FetH.vds*(self.FetH.tr + self.FetH.tf)*self.op.fsw
            self.FetH.ploss['dt'] = self.FetH.vfwd * self.FetH.idrms * self.GdH.deadtime * 2 * self.op.fsw
            self.FetH.ploss['coss'] = 0.5 * self.FetH.coss * self.FetH.vds**2 * self.op.fsw

            # Sync. Rectifier Loss (per Device)
            self.FetL.ploss['ohm'] = self.FetL.rdson * self.FetL.idrms**2
            self.FetL.ploss['dt'] = self.FetL.vfwd * self.FetL.idrms * self.GdL.deadtime * 2 * self.op.fsw

            # Inductor Loss
            self.Lout.ploss['ohm'] = self.Lout.dcr * self.Lout.irms**2

    def __repr__(self):
        s = "Sync Buck:\n"
        s += str(self.op)
        s += "----------------------------\n"
        s += "Inductor:\n"
        s += str(self.Lout)
        s += "----------------------------\n"
        s += "HS FET:\n"
        s += str(self.FetH)
        s += "----------------------------\n"
        s += "LS FET:\n"
        s += str(self.FetL)
        s += "----------------------------\n"
        s += "Output Capacitance:\n"
        s += f"Cout(ss) = {self.Cout.Css/1e-6:2.2f}uF\n"
        s += f"Cout(ac) = {self.Cout.Cac/1e-6:2.2f}uF\n"
        s += f"Cout:irms = {self.Cout.irms:2.2f}A\n"
        s += "----------------------------\n"
        s += "Input Capacitance:\n"
        s += f"Cin(ss) = {self.Cin.C/1e-6:2.2f}uF\n"
        s += f"Cin:irms = {self.Cin.irms:2.2f}A\n"
        if self._calc_losses == True:
            s += "----------------------------\n"
            s += "Gate Driver Losses:\n"
            s += f"[HS]: {self.GdH.ploss:2.2f}W\n"
            s += f"[LS]: {self.GdL.ploss:2.2f}W\n"
            s += "----------------------------\n"
            s += "HS FET Losses:\n"
            s += f"[Ohm]: {self.FetH.ploss['ohm']:2.2f}W\n"
            s += f"[Coss]: {self.FetH.ploss['coss']:2.2f}W\n"
            s += f"[SW]: {self.FetH.ploss['sw']:2.2f}W\n"
            s += f"[DT]: {self.FetH.ploss['dt']:2.2f}W\n"
            s += f"[Total]: {self.FetH.ploss_total:2.2f}W\n"
            s += "----------------------------\n"
            s += "LS FET Losses:\n"
            s += f"[Ohm]: {self.FetL.ploss['ohm']:2.2f}W\n"
            s += f"[Coss]: {self.FetL.ploss['coss']:2.2f}W\n"
            s += f"[SW]: {self.FetL.ploss['sw']:2.2f}W\n"
            s += f"[DT]: {self.FetL.ploss['dt']:2.2f}W\n"
            s += f"[Total]: {self.FetL.ploss_total:2.2f}W\n"
            s += "----------------------------\n"
            s += "Inductor Losses:\n"
            s += f"[Ohm]: {self.Lout.ploss['ohm']:2.2f}W\n"
        return s

