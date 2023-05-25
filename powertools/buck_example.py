from dcdcs import Requirements, SyncBuck 

req = Requirements(vin_min=30.0,
                   vin_max=60.0,
                   fsw=200e3,
                   il_rfactor=0.4,
                   vin_ripple = 30.0*0.01,
                   vout_ss_ripple = 0.1,
                   vout_ac_ripple = 0.3,
                   fc = 200e3/5.0)

buck1 = SyncBuck(vin=48.0, vout=12.0, pout= 25.0, requirements=req)
buck1.update(calc_losses=False)


