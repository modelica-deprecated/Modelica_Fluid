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
  Sources.SourceP Source(redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p0=3e5) 
  annotation (extent=[-80,20; -60,40]);
  Sources.SourceP SinkP1(p0=3e5, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
  annotation (extent=[60,22; 40,42]);
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
    M=0.1)            annotation (extent=[-42,20; -22,40]);
  Modelica.Blocks.Sources.Constant Constant1 
  annotation (extent=[-72, 64; -52, 84]);
  Modelica.Blocks.Sources.Sine Sine1(
    freqHz=0.1,
    startTime=0,
    offset=5e5,
    amplitude=4e5) 
                 annotation (extent=[92,72; 72,92]);
  Components.ValveLinear Valve(Kv=1e-5, redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
  annotation (extent=[4,22; 24,42]);
equation 
  connect(SinkP1.in_p, Sine1.y) annotation (points=[56,38.4; 56,82; 71,82],
      style(color=74, rgbcolor={0,0,127}));
  connect(Constant1.y, Valve.opening)     annotation (points=[-51,74; 14,
        74; 14,40], style(color=74, rgbcolor={0,0,127}));
  connect(Valve.port_b, SinkP1.port)     annotation (points=[25,32; 40,32],
                       style(color=69, rgbcolor={0,127,255}));
  connect(Valve.port_a, Pump1.outlet)     annotation (points=[3,32; -10,32; -10,
        37.4; -26,37.4],         style(color=69, rgbcolor={0,127,255}));
  connect(Pump1.inlet, Source.port) annotation (points=[-40,32.2; -50,32.2; -50,
        30; -60,30],           style(color=69, rgbcolor={0,127,255}));
end TestWaterPump;
