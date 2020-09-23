import subprocess
import fileinput

temp=5e7
temp=temp-.1e7
while temp<=5.1e7:
 temp=temp+.1e7
 Temp="{:.1e}".format(temp) 
 print(Temp)
 
 mdot=2e-9
 mdot=mdot-.1e-9
 while mdot<=2.1e-9:
  mdot=mdot+.1e-9
  Mdot="{:.1e}".format(mdot) 

  with open("../I_Files/I_Acc_KS.dat",'w') as g:
   g.write("ACCRETION RATE CONTROL:"+'\n')
   g.write("2"+'\n')
   g.write(Mdot+'\n')
   g.write("5e3"+'\n')
   g.write("1e30"+'\n')
   g.write("12.5"+'\n')
   g.write("3."+'\n')
   g.write("1e-4"+'\n')
   g.write("0.2"+'\n')
   g.write("0.7"+'\n')
   g.write('\n')
  g.close()
  with open("I.dat",'w') as f:
   f.write("PRINT OUT CONTROL:"+'\n')
   f.write("1	PSCREEN"+'\n')
   f.write("0.      DEBUG"+'\n')
   f.write("0       ISTEP DEBUG"+'\n')
   f.write("1	PTEFF"+'\n')
   f.write("0	PTEMP"+'\n')
   f.write("0	PSTAR"+'\n')
   f.write("1	IDUMP1"+'\n')
   f.write("111	IDUMP2"+'\n')
   f.write("421	IDUMP3"+'\n')
   f.write("8.5e3	TEMPMIN  "+'\n')
   f.write(Temp+"	TEMPINI"+'\n')
   f.write("ELECTRON SPECIFIC HEAT in OUTER CRUST CONTROL:------------"+'\n')
   f.write("0       ICVEL_NODEG"+'\n')
   f.write("LATTICE SPECIFIC HEAT CONTROL:------------"+'\n')
   f.write("3       ICVION"+'\n')
   f.write("EFFECTIVE MASS CONTROL:-----------------------------------"+'\n')
   f.write("5	EMNCO"+'\n')
   f.write("5	EMNCR"+'\n')
   f.write("3	EMP"+'\n')
   f.write("INITIAL ROTATIONAL PERIOD:--------------------------------"+'\n')
   f.write("0.1	p0: initial spin period"+'\n')
   f.write("TPRINT: times at which you select the profiles---------"+'\n')
   f.write("\n")
  f.close()

  with open("DCH_HZ08.in",'w') as g:
   g.write("  'NEW'"+'\n')
   g.write("BASIC MODEL FILES:"+'\n')
   g.write("  '../EOS/Crust/Crust_EOS_BSk20.dat'"+'\n')
   g.write("  '../EOS/BSk20_EOS_Acc_Fe.dat'"+'\n')
   g.write("  '../TOV/Profile/Prof_BSk20_Acc_Fe_1.0.dat'"+'\n')
   g.write("OTHER MODEL FILES:'"+'\n')
   g.write("  '../I_Files/I_Struct_1.4e14-4.0e11-1e8_fine.dat'"+'\n')
   g.write("  '../I_Files/I_Bound_Acc.dat'"+'\n')
   g.write("  '../I_Files/I_Pairing_0-0-0.dat'"+'\n')
   g.write("  '../I_Files/I_Neutrino_1.dat'"+'\n')
   g.write("  '../I_Files/I_Conduct_42.dat'"+'\n')
   g.write("  '../I_Files/I_Heat_pycno_bsk20.dat'"+'\n')
   g.write("  '../I_Files/I_Bfield_0.dat'"+'\n')
   g.write("  '../I_Files/I_Acc_KS.dat'"+'\n')
   g.write("OUTPUT FILES:"+'\n')
   g.write("  'I.dat'"+'\n')
   g.write("  'Teff_"+Temp+"_"+Mdot+"_0_Acc_5e3_1e8_1.0.dat'"+'\n')
   g.write("  'whatever'"+'\n')
   g.write("  'whatever'"+'\n')
   g.write("  'Stuff'"+'\n')
  g.close()
  

  subprocess.call('../NSCool.out', shell=True)

