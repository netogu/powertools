import powertools.components as comp
from powertools.dcdcs import PowerParams, SyncBuck


op = PowerParams()
op.vin = 48.0
op.vout = 12.0
op.pout = 100.0
op.fsw = 400e3
op.il_rfactor = 0.4
op.vin_ripple = op.vin * 0.01
op.vout_ss_ripple = 0.1
op.vout_ac_ripple = op.vout * 0.03
op.fc = op.fsw / 15.0
op.eff = 1.0

L = comp.Inductor()
L.name = '784325065'
L.value = 6.0e-6
L.dcr = 12.0e-3

GD1 = comp.GateDriver()
GD1.name = 'LM5148-HO'
GD1.vdrv = 5.0
GD1.rpull_down = 1.0
GD1.rpull_up = 1.0
GD1.ipk_source = 2.2
GD1.ipk_sink = 3.2
GD1.deadtime = 40e-9

GD2 = comp.GateDriver()
GD2.name = 'LM5148-LO'
GD2.vdrv = 5.0
GD2.rpull_down = 1.0
GD2.rpull_up = 1.0
GD2.ipk_source = 2.2
GD2.ipk_sink = 3.2
GD2.deadtime = 40e-9

FET1 = comp.Mosfet()
FET1.name = 'PSMN041-80YL'
FET1.rdson = 41e-3*1.5
FET1.coss = 100e-12
FET1.vgsth = 2.1
FET1.vmiller = 2.4
FET1.rg = 1.2
FET1.qgs = 2.5e-9
FET1.qgd = 4.3e-9
FET1.qsw = 0.8e-9 + 2.5e-9
FET1.qtot = 11.9e-9
FET1.vfwd = 1.2
FET1.qrr = 22e-9
FET1.trr = 22e-9

FET2 = comp.Mosfet()
FET2.name = 'PSMN014-80YL'
FET2.rdson = 15e-3*1.5
FET2.coss = 283e-12
FET2.vgsth = 2.1
FET2.vmiller = 2.8
FET2.qgs = 8.1e-9
FET2.qgd = 8.7e-9
FET2.qsw = 8.1e-9/2 + 8.7e-9
FET2.qtot = 29e-9
FET2.rg = 1.2
FET2.vfwd = 1.2
FET2.qrr = 25e-9
FET2.trr = 29e-9

# FET2 = comp.Mosfet()
# FET2.name = 'ISC0602NLS'
# FET2.rdson = 9.5e-3*1.5
# FET2.coss = 290e-12
# FET2.vmiller = 3.0
# FET2.qgs = 4.3e-9
# FET2.qgd = 3.9e-9
# FET2.qsw = 5.7e-9
# FET2.qtot = 14e-9
# FET2.rg = 1.2
# FET2.vgsth = 2.3
# FET2.vfwd = 1.0
# FET2.qrr = 20e-9
# FET2.trr = 28e-9

buck1 = SyncBuck(op=op, active_sw=FET1, passive_sw=FET2, hs_gatedrive=GD1, ls_gatedrive=GD2, inductor=L)
buck1.update(calc_losses=True)
print(buck1)


