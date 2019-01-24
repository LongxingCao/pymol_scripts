import cmd

def alignsg(sele):
   allsg = cmd.get_model("(("+sele+") and (name VSG))").atom
   allsf = cmd.get_model("(("+sele+") and (name S5,S6,S7,S8))").atom
   while allsg:
      closest = [9e9,None,None]
      for sga in allsg:
         for sfa in allsf:
            sg = Vec(sga.coord)
            sf = Vec(sfa.coord)
            if (sg-sf).length() < closest[0]:
               closest = [(sg-sf).length(),sga,sfa]
      sga,sfa = closest[1:]
      allsg.remove(sga)
      allsf.remove(sfa)
      if closest[0] > 10.0: break
      cmd.do("alter_state 1, (("+sele+") and (resi %s and name %s))"%(sga.resi,sga.name)+",x=%f"%sfa.coord[0])
      cmd.do("alter_state 1, (("+sele+") and (resi %s and name %s))"%(sga.resi,sga.name)+",y=%f"%sfa.coord[1])
      cmd.do("alter_state 1, (("+sele+") and (resi %s and name %s))"%(sga.resi,sga.name)+",z=%f"%sfa.coord[2])           

for i in range(1,10):
   alignsg("wire%i"%i)
   