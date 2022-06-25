#%%
import numpy as np
from scipy.optimize import fsolve


class r_divider(object):

    def __init__(self, Rtop, Rbot):
        self.Rtop = Rtop
        self.Rbot = Rbot
        self.calc_equivalent_r()
        self.calc_div_ratio()

    def calc_equivalent_r(self):
        self.Req = self.Rtop*self.Rbot/(self.Rtop + self.Rbot)

    def calc_div_ratio(self):
        self.dratio = self.Req / self.Rtop


#%%
class r_feedback(r_divider):

    def __init__(self, Rtop, Rbot, Cf, vref):
        super().__init__(Rtop, Rbot)
        self.Cf = Cf
        self.vref = vref
        self.calc_cutoff()
        self.calc_vout()
        self.calc_ibias()

    def calc_cutoff(self):
        self.fc = 1/(2*np.pi*self.Req*self.Cf+1e-20)

    def calc_vout(self):
        self.vout = self.vref / (self.dratio + 1e-20)
    
    def calc_cf(self, fc_hz):
        self.Cf = 1/(2*np.pi*self.Req*fc_hz+1e-20)

    def calc_ibias(self):
        self.rtot = self.Rtop + self.Rbot
        self.ibias = self.vout/(self.rtot + 1e-20)

    def evaluate(self):
        self.calc_equivalent_r()
        self.calc_div_ratio()
        self.calc_cutoff()
        self.calc_vout()
        self.calc_ibias()
       

    def __str__(self):

        string = "Rtop = {:2.2f}kOhm\n".format(self.Rtop/1e3)
        string = string + "Rbot = {:2.2f}kOhm\n".format(self.Rbot/1e3)
        string = string + "Cf = {:2.2f}nF\n".format(self.Cf/1e-9)
        string = string + "Vref = {:2.2f}V\n".format(self.vref)
        string = string + "Vset = {:2.2f}V\n".format(self.vout)
        string = string + "Fc = {:2.2f}Hz\n".format(self.fc)
        string = string + "Req = {:2.2f}kOhm\n".format(self.Req/1e3)
        string = string + "Rtot = {:2.2f}kOhm\n".format(self.rtot/1e3)
        string = string + "Ibias = {:2.2f}uA\n".format(self.ibias/1e-6)
        
        
        return string 
    






#%%
if __name__ == "__main__":
    
    RTOP = 100e3
    RBOT = 10e3
    CF = 2.2e-9
    VREF = 2.49
    VOUT = 12.0

    def func_solve_rtop(rtop):
        fb = r_feedback(rtop, RBOT, CF, VREF)
        return VOUT - fb.vout

    def func_solve_rbot(rbot):
        fb = r_feedback(RTOP, rbot, CF, VREF)
        return VOUT - fb.vout

    sol = fsolve(func_solve_rtop, RTOP)

    sol_fb = r_feedback(cval.e96_std_value(sol[0]), RBOT, CF, VREF)
    print(sol_fb)
