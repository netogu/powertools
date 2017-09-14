# Power Design Library - Components and Calculations

import numpy as num


class FET(object):

    '''
    -----------------------------------------
    The object represents an N-Channel Mosfet.
    -----------------------------------------
    '''

    def __init__(self, name='No Name'):
        self.Name = name

# MOSFET CHARACTERISTICS
    Vds_max = 0
    Vds = 0
    Ids_rms = 0
    Ids_peak = 0
    Rdson = 0

# LOSSES
    Pdis = 0
    Psw_loss = 0
    Pohm_loss = 0
    Pcoss_loss = 0
    Pqrr_loss = 0
    Pdt_loss = 0
    Pgate_drive = 0

#   GATE CHARACTERISTICS
    Vgs_drive = 0
    Vgs_miller = 0
    Vgs_th = 0

    Rg_internal = 0
    Rg_external = 0

    Qg = 0
    Qgd = 0
    Qgs = 0
    Qsw = 0
    Qrr = 0
    Coss = 0
    Ciss = 0

# BODY DIODE
    Qrr = 0
    Qrr_It = 0
    Qrr_norm = 0
    Vfwd_bodydiode = 0

# THERMAL CHARACTERISTICS
    t_rise = 0
    J_temp = 0


    t_deadtime = 0

    Supplier_link = ''

# FET RELATED FUNCTIONS

    def calc_ohm_loss(self):
        """Calculates the conductive loss of the MOSFET based on the RMS current"""

        self.Pohm_loss = self.Rdson * self.Ids_rms**2

    def calc_sw_trise(self):
        """Calculates the switching losses based on the buck topology"""

        ig_miller = ((self.Vgs_drive - self.Vgs_miller) / (self.Rg_external + self.Rg_internal))
        self.t_rise = self.Qsw/ig_miller
        print('trise = {:4.4f}ns'.format(self.t_rise/1e-9))

    def calc_jtemp(self):
        """Calculates the junction temperature"""
        self.Pdis = self.Psw_loss + self.Pohm_loss
        self.J_temp = self.Pdis * self.Rjc

    def calc_symetric_sw_loss(self,fsw):
        """Calculates switching loss assuming t_rise = t_fall"""
        self.Psw_loss = (self.Vds * self.Ids_peak)/2 * 2 * self.t_rise * fsw

    def calc_coss_loss(self, fsw):
        """Calculates the Coss loss"""
        self.Pcoss_loss = self.Coss * self.Vds**2*fsw/2

    def calc_qrr_loss(self, fsw):
        """Calculates the body diode reverse recovery loss"""
        self.Qrr_norm = self.Ids_peak/self.Qrr_It*self.Qrr
        self.Pqrr_loss = self.Qrr_norm * self.Vds * fsw

    def calc_total_losses(self):
        """Calculates the the total loss in the device"""
        self.Pdis = self.Pcoss_loss + self.Pohm_loss + self.Psw_loss + self.Pdt_loss + self.Pqrr_loss

    def calc_deadtime_loss(self, i_bd, fsw):
        """Calculates the deadtime loss caused by the body diode"""
        self.Pdt_loss = 2 * self.t_deadtime * fsw * self.Vfwd_bodydiode * i_bd

    def calc_bottom_fet_sw_loss_buck(self, iout, rdr, fsw):
        """Calculates the switching losses on the lower switch based on the buck topology"""
        k2r = num.log(self.Vgs_drive/(self.Vgs_drive - self.Vgs_miller)) \
            - num.log(self.Vgs_drive/(self.Vgs_drive - self.Vgs_th))

        k3r = num.log(self.Vgs_drive/(self.Vgs_drive - 0.9*self.Vgs_drive)) \
            - num.log(self.Vgs_drive/(self.Vgs_drive - self.Vgs_miller))

        k2f = num.log(self.Vgs_miller/self.Vgs_th)
        k3f = num.log(0.9*self.Vgs_drive/self.Vgs_miller)

        t2r = k2r*(rdr + self.Rg_internal + self.Rg_external) * self.Ciss
        t3r = k3r*(rdr + self.Rg_internal + self.Rg_external) * self.Ciss

        t2f = k2f*(rdr + self.Rg_external + self.Rg_internal) * self.Ciss
        t3f = k3f*(rdr + self.Rg_external + self.Rg_internal) * self.Ciss

        t2 = t2r + t2f
        t3 = t3r + t3f

        return (t2*self.Vfwd_bodydiode + t3*(self.Vfwd_bodydiode+iout*1.1*self.Rdson)/2)*iout*fsw

    def print_losses(self):

        print('---------------------------')
        print('FET = {}'.format(self.Name))
        print('Conduction Loss = {:4.3f} W '.format(self.Pohm_loss))
        print('Switching Loss  = {:4.3f} W '.format(self.Psw_loss))
        print('Coss Loss = {:4.3f} W '.format(self.Pcoss_loss))
        print('Deadtime Loss = {:4.3f} W '.format(self.Pdt_loss))
        print('Qrr Loss = {:4.3f} W '.format(self.Pqrr_loss))

        print('Total Loss = {:4.3f} W '.format(self.Pdis))
        print('---------------------------\n')


class Inductor(object):

    '''
    -----------------------------------------
    The object represents a Power Inductor.
    -----------------------------------------
    '''

    def __init__(self, name='No Name'):
        self.Name = name

    # General Properties
    DCR = 0
    DCR_temp = 0
    L = 0
    Ilsat = 0
    Ilrms_max = 0
    Il_rms = 0
    Il_avg = 0
    Il_peak = 0
    Il_valley = 0
    Il_pp = 0
    Volt_sec = 0
    Volt_usec = 0
    Rth = 0
    temp_rise = 0

    # Loss Variables
    Pdis = 0
    Pdc_loss = 0
    Pac_loss = 0
    Pcore_loss = 0

    # Magnetic Material Properties
    B_max = 0
    Bpk_gauss = 0
    Bpk = 0
    Bsat = 0
    kf = 0
    kb = 0
    K0 = 0
    K1 = 0
    f_eff = 0  # correction for steinmetz eq.

    def calc_effective_freq(self, f0, duty_cycle):
        """Calculates frequency correctio to be used in the sine Steinmetz equation"""
        self.f_eff = f0/(2*num.pi*(duty_cycle-duty_cycle**2))

    def calc_dcr_temp(self,temp_amb):
        """Calculates the effective DCR at temperature temp_amb"""
        self.DCR_temp = self.DCR*(274.5+temp_amb)/259.5

    def calc_volt_sec(self, VL_a,VL_b,duty_cycle,fsw):
        """Calculates the Volt-sec at the inductor"""
        self.Volt_usec = (VL_a-VL_b)*duty_cycle/fsw/1e-6
        self.Volt_sec = (VL_a-VL_b)*duty_cycle/fsw


    def calc_bpk_gauss(self, etckt,et100):
        """Calculates the flux density peak in gauss"""
        self.Bpk_gauss = etckt/et100*100

    def calc_core_loss(self,f0):
        """Calculates the inductor core losses"""
        self.Pcore_loss = self.K0*self.f_eff**(self.kf-1)*self.Bpk_gauss**self.kb*f0*10**(-14)

    def calc_ac_loss(self,f0):
        """Calculates the inductor AC losses in the copper"""
        self.Pac_loss = self.K1*self.Il_pp**2*num.sqrt(f0)*self.DCR_temp

    def calc_dc_loss(self):
        """Calculates the inductor DC losses in the copper"""
        self.Pdc_loss = self.Il_rms**2 * self.DCR_temp

    def calc_inductor_currents(self):
        """Calculates the current waveform based on a square wave voltage input"""
        self.Il_pp = self.Volt_sec/self.L
        self.Il_peak = self.Il_avg + self.Il_pp/2
        self.Il_valley = self.Il_avg - self.Il_pp/2
        self.Il_rms = num.sqrt((self.Il_peak**2+self.Il_valley**2+self.Il_peak*self.Il_valley)/3.0)


class Capacitor(object):
    def __init__(self, c, esr, rmsrate):
        self.C = c
        self.ESR = esr
        self.Ic_rms_rating = rmsrate
    Ic_peak = 0
    Ic_rms = 0
    Pesr_loss = 0
    Pdis = 0

    def calc_cu_loss(self):
        self.Pesr_loss = self.Ic_rms * self.ESR


class Diode(object):

    Ploss = 0
    If_rms = 0
    If_peak = 0
    If_avg = 0
    Vrrm = 0
    Qrr = 0
    Coss = 0
    Vfwd = 0

    def calc_ploss(self):
        self.Ploss = self.Vf * self.Id_rms


class Source(object):

    Vmax = 0
    Vtyp = 0
    Vmin = 0
    Pdis = 0


class Load(object):
    def __init__(self, i=0, v=0):

        self.V = v
        self.I_tdp = i
        self.I_edp = 0
        self.P_tdp = self.I_tdp * self.V
        self.P_edp = self.I_edp * self.V
        self.R = self.V / (self.I_tdp + 1e-12)


class Buck_DCDC(object):
    def __init__(self, vin=Source(), vout=Load(), il_rf=0, n_phases=0, fsw=0):

        # General Buck DC/DC Characteristics (No components)
        self.N_Phases = n_phases

        # Input / Output Voltage Characteristics
        self.Vin = num.array([vin.Vmin, vin.Vtyp, vin.Vmax])
        self.Vout = vout.V
        self.Vout_ripple_dc = 0
        self.Vout_ripple_ac = 0
        self.Vout_pp = 0
        self.Vin_ripple_dc = 0
        self.Vin_ripple_ac = 0
        self.Vin_pp = 0

        # Input / Output Current Characteristics
        self.Iout_tdp = vout.I_tdp
        self.Iout_edp = vout.I_edp
        self.Iin_tdp = 0
        self.Iin_edp = 0
        self.Iin_rms = 0
        self.load_step = 0
        self.Il_rf = il_rf
        self.Iorip_norm = 0
        self.Icin_norm = 0

        # Input / Output Power Characteristics
        self.P_tdp = vout.P_tdp
        self.P_edp = vout.P_edp
        self.P_tdp = self.Iout_tdp * self.Vout
        self.P_edp = self.Iout_edp * self.Vout

        # DC/DC Switching and Dynamic Characteristics
        self.Fsw = fsw
        self.Tsw = 1/(self.Fsw+1e-12)
        self.ton = 0
        self.toff = 0
        self.fc = 0
        self.t_vmin = 0
        self.t_vmax = 0
        self.DC = 0
        self.Ton_req = 0

        # Output Capacitor Characteristics
        self.Cout_bulk = 0
        self.Cout_cer = 0
        self.Cout_bulk_esr  = 0
        self.Cout_cer_esr = 0
        self.Cout_esr_total = 0
        self.Cout_total = 0
        self.Cout_n_bulk = 0
        self.Cout_n_cer = 0
        self.Vout_ac_cap = 0
        self.Vout_ac_esr = 0

        # Input Capacitor Characteristics

        self.Cin_total = 0
        self.Icin_rms = 0

        # Gate Driver Characteristics
        self.Rdr = 0

    def calculate_required_inductance(self, l=Inductor(), vin=0):
        """ Calculate how much is required to meet the given ripple factor"""

        self.DC = self.Vout / (vin + 1e-22)
        l.Il_avg = self.Iout_tdp / self.N_Phases
        l.Il_pp = l.Il_avg * self.Il_rf
        l.L = (vin - self.Vout)/(l.Il_pp * self.Fsw) * self.DC
        self.Ton_req = 1/self.Fsw * self.DC

        print('---------------------------')
        print('Calculated Inductance L = {:.2f}uH'.format(l.L / 1e-6))
        print('diL = {:4.2f}A'.format(l.Il_pp))
        print('<iL> = {:4.2f}A'.format(l.Il_avg))
        print('Duty Cycle = {:.02f}% '.format(self.DC * 100))
        print('t_on = {:0.2f}ns'.format(self.Ton_req / 1e-9))
        print('---------------------------')

    def calculate_currents_sync_rec_tdp(self, l=Inductor(), qtop=FET(), qbot=FET(), vin=0):
        """ Calculate the inductor and switch currents (TDP)"""
        self.DC = self.Vout / (vin + 1e-12)
        l.calc_volt_sec(vin, self.Vout, self.DC, self.Fsw)
        l.Il_pp = (vin - self.Vout) * self.DC  / (l.L * self.Fsw)
        l.Il_avg = self.Iout_tdp / self.N_Phases
        l.Il_peak = l.Il_avg + l.Il_pp / 2.0
        l.Il_valley = l.Il_avg - l.Il_pp / 2.0
        l.Il_rms = num.sqrt((l.Il_peak ** 2 + l.Il_valley ** 2 + l.Il_peak * l.Il_valley) / 3.0)

        qtop.Ids_peak = l.Il_peak
        qtop.Ids_rms = l.Il_avg * num.sqrt(self.DC)
        qtop.Vds = vin

        qbot.Ids_peak = l.Il_peak
        qbot.Ids_rms = l.Il_avg * num.sqrt(1 - self.DC)
        qbot.Vds = vin

        # PRINT RESULTS
        print('++++++++++++ TDP +++++++++++++++')
        print('Vin = {:.2f}'.format(vin))
        print('Duty Cycle = {:.2f}%'.format(self.DC * 100))
        print('No. of phases = {}'.format(self.N_Phases))
        print('Output Power = {:4.2f}W'.format(self.P_tdp))
        print('Fsw = {:.2f}kHz'.format(self.Fsw / 1e3))
        print('Chosen Inductance = {:.2f}uH'.format(l.L / 1e-6))
        print('---------------------------')
        print('Inductor rms current:  iLrms = {:4.3f} A rms'.format(l.Il_rms))
        print('Inductor peak current: iLpeak = {:4.3f} A '.format(l.Il_peak))
        print('Inductor valley current: iLvalley = {:4.3f} A '.format(l.Il_valley))
        print('Inductor peak-peak current: iLpp = {:4.3f} A '.format(l.Il_pp))
        print('---------------------------\n')

        print('Top MOSFET Vds = {:4.2f}V'.format(qtop.Vds))
        print('Bottom MOSFET Vds = {:4.2f}V'.format(qbot.Vds))
        print('---------------------------')
        print('Top MOSFET rms current:  ids rms  = {:4.3f} A rms'.format(qtop.Ids_rms))
        print('Top MOSFET peak current: ids peak = {:4.3f} A '.format(qtop.Ids_peak))
        print('Bottom MOSFET rms current:  ids rms  = {:4.3f} A rms'.format(qbot.Ids_rms))
        print('Bottom MOSFET peak current: ids peak = {:4.3f} A '.format(qbot.Ids_peak))
        print('---------------------------\n')

    def calculate_currents_passive_rec_tdp(self, l=Inductor(), qtop=FET(), diode=Diode(), vin=0):
        """ Calculate the inductor and switch currents (TDP)"""

        self.DC = self.Vout / (vin + 1e-12)
        l.calc_volt_sec(vin, self.Vout, self.DC, self.Fsw)
        l.Il_pp = (vin - self.Vout) * self.DC / (l.L * self.Fsw)
        l.Il_avg = self.Iout_tdp / self.N_Phases
        l.Il_peak = l.Il_avg + l.Il_pp / 2
        l.Il_valley = l.Il_avg - l.Il_pp / 2
        l.Il_rms = num.sqrt((l.Il_peak ** 2 + l.Il_valley ** 2 + l.Il_peak * l.Il_valley) / 3.0)

        qtop.Ids_peak = l.Il_peak
        qtop.Ids_rms = l.Il_avg * num.sqrt(self.DC)
        qtop.Vds = vin

        diode.If_rms = l.Il_avg * num.sqrt(1 - self.DC)
        diode.If_avg = l.Il_avg * (1-self.DC)
        diode.If_peak = l.Il_peak
        diode.Vrrm = vin


        # PRINT RESULTS
        print('++++++++++++ TDP +++++++++++++++')
        print('Vin = {:.2f}'.format(vin))
        print('Duty Cycle = {:.2f}%'.format(self.DC * 100))
        print('No. of phases = {}'.format(self.N_Phases))
        print('Output Power = {:4.2f}W'.format(self.P_tdp))
        print('Fsw = {:.2f}kHz'.format(self.Fsw/1e3))
        print('Chosen Inductance = {:.2f}uH'.format(l.L/1e-6))
        print('---------------------------')
        print('Inductor rms current:  iLrms = {:4.3f} A rms'.format(l.Il_rms))
        print('Inductor peak current: iLpeak = {:4.3f} A '.format(l.Il_peak))
        print('Inductor valley current: iLvalley = {:4.3f} A '.format(l.Il_valley))
        print('Inductor peak-peak current: iLpp = {:4.3f} A '.format(l.Il_pp))
        print('---------------------------\n')

        print('Top MOSFET Vds = {:4.2f}V'.format(qtop.Vds))
        print('Diode Blocking Voltage = {:4.2f}V'.format(diode.Vrrm))
        print('---------------------------')
        print('Top MOSFET rms current:  ids rms  = {:4.3f} A rms'.format(qtop.Ids_rms))
        print('Top MOSFET peak current: ids peak = {:4.3f} A '.format(qtop.Ids_peak))
        print('Diode avg current:  If avg  = {:4.3f} A rms'.format(diode.If_avg))
        print('Diode peak current: If peak = {:4.3f} A '.format(diode.If_peak))
        print('---------------------------\n')

    def calculate_currents_edp(self, l=Inductor(), qtop=FET(), qbot=FET(), vin=0):

        duty = self.Vout / (vin + 1e-12)
        l.calc_volt_sec(vin, self.Vout, duty, self.Fsw)
        l.Il_pp = (vin - self.Vout) * duty / (l.L * self.Fsw)
        l.Il_avg = self.Iout_edp / self.N_Phases
        l.Il_peak = l.Il_avg + l.Il_pp / 2
        l.Il_valley = l.Il_avg - l.Il_pp / 2
        l.Il_rms = num.sqrt((l.Il_peak ** 2 + l.Il_valley ** 2 + l.Il_peak * l.Il_valley) / 3.0)

        qtop.Ids_peak = l.Il_peak
        qtop.Ids_rms = l.Il_avg * num.sqrt(duty)
        qtop.Vds = vin

        qbot.Ids_peak = l.Il_peak
        qbot.Ids_rms = l.Il_avg * num.sqrt(1 - duty)
        qbot.Vds = vin

        # PRINT RESULTS
        print('+++++++++++++ EDP ++++++++++++++')
        print('Vin = {:.2f}'.format(vin))
        print('Duty Cycle = {:.2f}'.format(duty))
        print('No. of phases = {}'.format(self.N_Phases))
        print('Output Power = {:4.2f}W'.format(self.P_edp))
        print('Chosen Inductance = {:.2f}uH'.format(l.L / 1e-6))
        print('---------------------------')
        print('Inductor rms current:  iLrms = {:4.3f} A rms'.format(l.Il_rms))
        print('Inductor peak current: iLpeak = {:4.3f} A '.format(l.Il_peak))
        print('Inductor valley current: iLvalley = {:4.3f} A '.format(l.Il_valley))
        print('Inductor peak-peak current: iLpp = {:4.3f} A '.format(l.Il_pp))
        print('---------------------------\n')

        print('Top MOSFET Vds = {:4.2f}V'.format(qtop.Vds))
        print('Bottom MOSFET Vds = {:4.2f}V'.format(qbot.Vds))
        print('---------------------------')
        print('Top MOSFET rms current:  ids rms  = {:4.3f} A rms'.format(qtop.Ids_rms))
        print('Top MOSFET peak current: ids peak = {:4.3f} A \n'.format(qtop.Ids_peak))
        print('Bottom MOSFET rms current:  ids rms  = {:4.3f} A rms'.format(qbot.Ids_rms))
        print('Bottom MOSFET peak current: ids peak = {:4.3f} A '.format(qbot.Ids_peak))
        print('---------------------------\n')

    def calculate_losses_sync_rec_tdp(self, qtop=FET(), qbot=FET(), l=Inductor()):
        """ Calculate the DC and AC FET losses"""

        # ------------   Active Switch ----------------------

        # Active Switch Conduction Loss
        qtop.calc_ohm_loss()

        # Active Switch Switching Loss
        qtop.calc_sw_trise()

        # gdriver_Isource = 1.8
        # gdriver_Isink = 3.5
        # qtop.t_rise = Qtop.Qsw / gdriver_Isource
        # qtop.t_fall = Qtop.Qsw / gdriver_Isink

        # Active Switch switching Loss (based on triangle)
        qtop.calc_symetric_sw_loss(self.Fsw)

        # Active Switch Coss Loss
        qtop.Pcoss_loss = (qtop.Coss + qbot.Coss/2.0) * qtop.Vds ** 2 * self.Fsw / 2

        # Top MOSFET Qrr losses
        # Normalize Qrr for peak current
        qbot.Qrr_n = l.Il_valley / qbot.Qrr_Ipeak * qbot.Qrr
        # If ta = tb
        qtop.Pqrr_loss = 4 / 6.0 * qbot.Qrr_n * qtop.Vds * self.Fsw

        # Sum Total MOSFET Losses
        qtop.calc_total_losses()
        # Qtop.Pdis = Qtop.PQrr_loss+Qtop.Pcoss_loss+Qtop.Pohm_loss+Qtop.Psw_loss

        # ------------   Sync Switch ----------------------

        # Conduction Loss
        qbot.calc_ohm_loss()

        # SW losses (very low)
        qbot.Psw_loss = qbot.calc_bottom_fet_sw_loss_buck(self.Iout_tdp, self.Rdr, self.Fsw)

        # Deadtime Loss
        qbot.calc_deadtime_loss(l.Il_avg, self.Fsw)

        # Coss Loss
        qbot.Pcoss_loss = (qbot.Coss/2.0) * qtop.Vds ** 2 * self.Fsw / 2

        # Sum all losses
        qbot.calc_total_losses()

    def calculate_losses_passive_rec_tdp(self, qtop=FET(), diode=Diode()):
        """ Calculate the DC and AC FET losses"""

        # ------------   Active Switch ----------------------

        # Active Switch Conduction Loss
        qtop.calc_ohm_loss()

        # Active Switch Switching Loss
        qtop.calc_sw_trise()

        # gdriver_Isource = 1.8
        # gdriver_Isink = 3.5
        # qtop.t_rise = Qtop.Qsw / gdriver_Isource
        # qtop.t_fall = Qtop.Qsw / gdriver_Isink

        # Active Switch switching Loss (based on triangle)
        qtop.calc_symetric_sw_loss(self.Fsw)

        # Active Switch Coss Loss
        qtop.calc_coss_loss(self.Fsw)

        # Top MOSFET Qrr losses
        # Normalize Qrr for peak current
        # qbot.Qrr_n = l.Il_valley / qbot.Qrr_Ipeak * diode.Qrr
        # If ta = tb
        # qtop.Pqrr_loss = 4 / 6.0 * qbot.Qrr_n * qtop.Vds * self.Fsw

        # Sum Total MOSFET Losses
        qtop.calc_total_losses()
        # Qtop.Pdis = Qtop.PQrr_loss+Qtop.Pcoss_loss+Qtop.Pohm_loss+Qtop.Psw_loss

        # ------------   Diode ----------------------

        diode.Ploss = diode.Vfwd * diode.If_avg

    def calculate_equivalent_output_cap_esr(self):

        """ Calculate equivalent Capacitance and ESR for output filter"""

        self.Cout_total = self.Cout_n_bulk * self.Cout_bulk + self.Cout_n_cer * self.Cout_cer
        bulk_eq_esr = self.Cout_bulk_esr / self.Cout_n_bulk
        cer_eq_esr = self.Cout_cer_esr / self.Cout_n_cer
        self.Cout_esr_total = bulk_eq_esr * cer_eq_esr / ( bulk_eq_esr + cer_eq_esr)

    def calculate_output_capacitance(self, l=Inductor()):

        f0 = self.Fsw / num.pi * num.sqrt(2 * self.Vout_ripple_dc / (self.Vout * (1 - self.DC)))
        cout_ss = 1 / ((2 * num.pi * f0)**2 * l.L)
        cout_ac = self.load_step / (2 * num.pi * self.fc * self.Vout_ripple_ac)

        print('---------------------------------------------------------')
        print('LC cutoff fo = {:.2f}kHz to get {:.2f}mV of output ripple'.format(f0 / 1e3, self.Vout_ripple_dc / 1e-3))
        print('Capacitance to meet ss ripple = {:.2f}uF'.format(cout_ss/1e-6))
        print('Capacitance to meet AC ripple = {:.2f}uF with a load step of {:.2f}A'
              .format(cout_ac / 1e-6, self.load_step))
        print('----------------------------------------------------------\n')

    def calculate_eff_output_current_ripple_factor(self):

        # Multiphase Ripple Cancellation
        m = num.floor(self.DC * self.N_Phases)
        numerator = self.N_Phases * ((self.DC - m / self.N_Phases) * ((1 + m) / self.N_Phases - self.DC))
        denominator = ((1 - self.DC) * self.DC)
        self.Iorip_norm =  numerator / denominator

    def calculate_eff_input_rms_current_factor(self):

        # Input RMS / Load Current

        m = num.floor(self.DC * self.N_Phases)

        self.Icin_norm = num.sqrt((self.DC - m/self.N_Phases) * ((1+m)/self.N_Phases - self.DC))

    def calculate_output_voltage_ripple(self, l=Inductor()):

        self.ton = self.DC / self.Fsw
        self.toff = (1 - self.DC) / self.Fsw

        r_esr = self.Cout_esr_total
        real_ripple = l.Il_pp  * r_esr * (1 - (self.t_vmax / self.toff + self.t_vmin / self.ton))
        reactive_ripple = l.Il_pp  / (2 * self.Cout_total) * (self.t_vmax + self.t_vmin -
                                                             (self.t_vmax**2 / self.toff + self.t_vmin**2/self.ton))

        self.Vout_pp = real_ripple + reactive_ripple

    def calculate_input_capacitance(self,l=Inductor()):

        q_cin = self.DC * self.Iout_edp * (1 - self.DC) * self.Tsw
        self.Cin_total = q_cin / self.Vin_ripple_dc

        print('---------------------------------------------------------')
        print('Capacitance to meet Input Voltage ripple = {:.2f}uF'.format(self.Cin_total / 1e-6))
        print('----------------------------------------------------------\n')

    def calculate_cin_rms_current(self, l=Inductor()):

        self.Icin_rms = self.Iout_edp * num.sqrt(self.DC * (1-self.DC) +
                        1.0/12.0 * (self.Vout / (l.L * self.Fsw *self.Iout_edp))**2 * (1 - self.DC)**2 * self.DC)

















            # def plot_loss_piechart(losses,colors,labels,title_str):
#     """Displays all losses in a PieChart"""
#
#     losses_norm = [ v/min(losses) for v in losses]
#
#     def make_autopct(values):
#         def my_autopct(pct):
#             total = sum(values)
#             val = pct*total/100.0
#             return '{p:.3f}%\n({v:.3f}W)'.format(p=pct,v=val)
#         return my_autopct
#
#     figure(figsize=(6,6))
#     plot=pie(losses_norm, colors=colors, labels=labels,shadow='True',autopct=make_autopct(losses))
#     plot=axis('equal')
#     plot=title(title_str)