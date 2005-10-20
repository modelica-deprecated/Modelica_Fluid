model TestWaterPump "Test case for WaterPump" 
  extends Modelica.Icons.Example;
  
annotation (
  Diagram,
  experiment(StopTime=10, Tolerance=1e-006),
  Documentation(info="<HTML>
<p>This model tests the <tt>Pump</tt> model with the check valve option active.
<p>The sink pressure is varied sinusoidally with a period of 10 s, so as to operate the pump in all the possible working conditions, including stopped flow.
<p>
Simulation Interval = [0...10] sec <br> 
Integration Algorithm = DASSL <br>
Algorithm Tolerance = 1e-6 
<p><b>Revision history:</b></p>
<ul>
<li><i>5 Feb 2004</i>
    by <a href=\"mailto:francesco.schiavo@polimi.it\">Francesco
Schiavo</a>:<br>
       First release.</li>
</ul>
</HTML>"));
  Sources.FixedAmbient_pTX Source(redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=3e5) 
  annotation (extent=[-100,20; -80,40]);
  Sources.PrescribedAmbient_pTX SinkP1(p=3e5, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
  annotation (extent=[36,26; 16,46]);
  Components.Pump Pump1(
    rho0=1000,
    OpPoints=true,
    pin_start=1e5,
    pout_start=4e5,
    hstart=1e5,
    P_cons={800,1800,2000},
    head_nom={60,30,0},
    q_nom={0,0.001,0.0015},
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    redeclare package SatMedium = Modelica.Media.Water.StandardWater,
    ComputeNPSHa=true,
    CheckValve=true,
    M=0.1, 
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Components.PumpCharacteristics.LinearFlowCharacteristic 
        (q_nom=[0,1])) 
                      annotation (extent=[-66,14; -32,46]);
  Modelica.Blocks.Sources.Constant Constant1 
  annotation (extent=[-60,60; -40,80]);
  Modelica.Blocks.Sources.Sine Sine1(
    freqHz=0.1,
    startTime=0,
    offset=5e5,
    amplitude=4e5) 
                 annotation (extent=[80,60; 60,80]);
  Components.ValveLinear Valve(Kv=1e-5, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
  annotation (extent=[-16,26; 2,46]);
equation 
  connect(Constant1.y, Valve.opening)     annotation (points=[-39,70; -7,70; -7,
        44],        style(color=74, rgbcolor={0,0,127}));
  connect(Valve.port_b, SinkP1.port)     annotation (points=[2.9,36; 15,36],
                       style(color=69, rgbcolor={0,127,255}));
  connect(Valve.port_a, Pump1.outlet)     annotation (points=[-16.9,36; -26,36; 
        -26,41.84; -38.8,41.84], style(color=69, rgbcolor={0,127,255}));
  connect(Pump1.inlet, Source.port) annotation (points=[-62.6,33.52; -62.6,30; 
        -79,30],               style(color=69, rgbcolor={0,127,255}));
  connect(Sine1.y, SinkP1.p_in) annotation (points=[59,70; 50,70; 50,42; 38,42], 
      style(color=74, rgbcolor={0,0,127}));
end TestWaterPump;
