import powertools.components as comp
from powertools.dcdcs import PowerParams, SyncBuck


op = PowerParams()
op.vin = 48.0
op.vout = 12.0
op.pout = 120.0
op.fsw = 400e3
op.il_rfactor = 0.4
op.vin_ripple = op.vin * 0.01
op.vout_ss_ripple = 0.1
op.vout_ac_ripple = op.vout * 0.03
op.fc = op.fsw / 15.0
op.eff = 1.0


buck1 = SyncBuck(op=op)
buck1.update(calc_losses=False)
print(buck1)


