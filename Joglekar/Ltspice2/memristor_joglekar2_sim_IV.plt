[Transient Analysis]
{
   Npanes: 2
   Active Pane: 1
   {
      traces: 1 {336592898,0,"Ix(U1:TE)"}
      Parametric: "V(n001)"
      X: (' ',2,-1.2,0.02,1.2)
      Y[0]: ('m',1,-0.001,0.0002,0.0012)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Amps: ('m',0,0,1,-0.001,0.0002,0.0012)
      Log: 0 0 0
      NeyeDiagPeriods: 0
      TeyeEnd: 3
   },
   {
      traces: 1 {268959747,0,"V(nc_01)"}
      Parametric: "V(n001)"
      X: (' ',2,-1.2,0.02,1.2)
      Y[0]: (' ',2,0.21,0.07,1.05)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,2,0.21,0.07,1.05)
      Log: 0 0 0
      NeyeDiagPeriods: 0
      TeyeEnd: 3
   }
}
