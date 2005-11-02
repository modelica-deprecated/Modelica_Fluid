package DrumBoiler 
  "Drum boiler example, see Franke, Rode, Krueger: On-line Optimization of Drum Boiler Startup, 3rd International Modelica Conference, Linkoping, 2003" 
  
  model DrumBoilerSimulation "Simulate start-up of DrumBoiler" 
    extends Modelica.Icons.Example;
    DrumBoiler drumBoiler            annotation (extent=[-20, -40; 40, 20]);
    Modelica.Blocks.Sources.TimeTable q_F_Tab(table=[0, 0; 3600, 400; 7210,
          400]) annotation (extent=[-80,0; -60,20]);
    Modelica.Blocks.Sources.TimeTable Y_Valve_Tab(table=[0, 1; 3600, 1; 7210,
           1]) annotation (extent=[-80, -40; -60, -20]);
    annotation (
      Diagram,
      experiment(StopTime=7200),
      Documentation(info="<HTML>
<p>
Apply a ramp to fuel input and hold outlet valve open.
Simulate for 7200 seconds.
</p>
</HTML>"));
  equation 
    connect(q_F_Tab.y, drumBoiler.q_F)       annotation (points=[-59,10; -40,10; 
          -40,-31; -21.35,-31],       style(rgbcolor={0,0,127}));
    connect(Y_Valve_Tab.y, drumBoiler.Y_Valve)       annotation (points=[-59,-30; 
          -44,-30; -44,-37; -21.35,-37],         style(
        rgbcolor={0,0,127},
        fillColor=7,
        fillPattern=1));
  end DrumBoilerSimulation;
  
  model DrumBoiler 
    "Complete drum boiler model, including evaporator and supplementary components" 
    import Modelica_Fluid;
    
    import Modelica.SIunits.Conversions.*;
    
    Modelica_Fluid.Components.Evaporator evaporator(
      m_D=300e3,
      cp_D=500,
      V_t=100,
      p_start=1e5,
      V_l_start=67,
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      initOption=Modelica_Fluid.Types.InitTypes.InitialValues) 
                          annotation (extent=[-46,-29; -26,-9]);
    annotation (
      uses(Modelica_Fluid(version="0.72")),
      Diagram,
      Coordsys(
        extent=[-100, -100; 100, 100],
        grid=[1, 1],
        component=[20, 20]),
      Icon(
        Rectangle(extent=[-100, 100; 100, -100], style(fillColor=7)),
        Text(
          extent=[-151, 165; 138, 102],
          style(fillColor=7, fillPattern=1),
          string="%name"),
        Text(
          extent=[-79, 67; 67, 21],
          string="drum",
          style(
            color=0,
            fillColor=7,
            fillPattern=1)),
        Text(
          extent=[-90, -14; 88, -64],
          string="boiler",
          style(
            color=0,
            fillColor=7,
            fillPattern=1))));
    Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
      annotation (extent=[-46,-63; -26,-43],   rotation=90);
    Modelica.Blocks.Interfaces.RealInput q_F(redeclare type SignalType = Real (
           unit="MW")) "Thermal power to the evaporator" 
      annotation (extent=[-109,-65; -100,-75]);
    Modelica.Blocks.Interfaces.RealInput Y_Valve(redeclare type SignalType = 
          Real (unit="1")) 
      annotation (extent=[-109,-95; -100,-85]);
    Modelica_Fluid.Sources.FixedAmbient sink(p=from_bar(0.5), redeclare package
        Medium = Modelica.Media.Water.StandardWater) 
      annotation (extent=[77,-29; 97,-9],   rotation=180);
    Modelica_Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium 
        = Modelica.Media.Water.StandardWater) 
      annotation (extent=[40,-29; 20,-9],  rotation=180);
    Modelica_Fluid.Sensors.Temperature temperature(redeclare package Medium = 
          Modelica.Media.Water.StandardWater) 
      annotation (extent=[-10,-9; 10,-29]);
    Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
          Modelica.Media.Water.StandardWater) 
      annotation (extent=[10,14; 30,34]);
    Modelica.Blocks.Continuous.PI controller(T=120, k=10) 
      annotation (extent=[-51,23; -65,37]);
    Modelica_Fluid.Sources.PrescribedMassFlowRate_hX pump(
                                             h=5e5, redeclare package Medium = 
          Modelica.Media.Water.StandardWater) 
      annotation (extent=[-80,-30; -60,-10]);
    Modelica.Blocks.Math.Feedback feedback 
      annotation (extent=[-26,20; -46,40]);
    Modelica.Blocks.Sources.Constant levelSetPoint(k=67) 
      annotation (extent=[-43,50; -30,63]);
    Modelica.Blocks.Interfaces.RealOutput T_S(redeclare type SignalType = 
          Real (unit="degC")) 
      annotation (extent=[100,56; 108,64]);
    Modelica.Blocks.Interfaces.RealOutput p_S(redeclare type SignalType = 
          Real (unit="bar")) 
      annotation (extent=[100,20; 108,28]);
    Modelica.Blocks.Interfaces.RealOutput qm_S(redeclare type SignalType = 
          Modelica.SIunits.MassFlowRate) 
      annotation (extent=[100,-4; 108,4],    rotation=0);
    Modelica.Blocks.Interfaces.RealOutput V_l(redeclare type SignalType = 
          Modelica.SIunits.Volume) 
      annotation (extent=[100,88; 108,96]);
  public 
    Modelica.Blocks.Math.Gain MW2W(k=1e6) 
      annotation (extent=[-95,-75.5; -85,-64.5]);
    Modelica.Blocks.Math.Gain Pa2bar(k=1e-5) annotation (extent=[37,19; 47,29]);
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin K2degC 
      annotation (extent=[38,55; 48,65]);
    Modelica.Blocks.Nonlinear.Limiter limiter(uMin=0, uMax=500) 
      annotation (extent=[-85,23; -71,37], rotation=180);
    Modelica_Fluid.Components.ValveLinear SteamValve(redeclare package Medium 
        = Modelica.Media.Water.StandardWater, Kv=2e-5) 
      annotation (extent=[53,-10; 66,-26]);
  equation 
    connect(furnace.port, evaporator.heatPort) 
      annotation (points=[-36,-43; -36,-30],   style(color=42));
    connect(controller.u,feedback.y) 
      annotation (points=[-49.6,30; -45,30], style(rgbcolor={0,0,127}));
    connect(feedback.u2,      evaporator.V) 
      annotation (points=[-36,22; -36,-9; -32,-8],
                                            style(rgbcolor={0,0,127}));
    connect(levelSetPoint.y,feedback.u1)             annotation (points=[
          -29.35,56.5; -22,56.5; -22,30; -28,30],
                                           style(rgbcolor={0,0,127}));
    connect(massFlowRate.m_flow, qm_S) 
      annotation (points=[30,-8; 30,0; 104,0],        style(rgbcolor={0,0,127}));
    connect(evaporator.V, V_l) 
      annotation (points=[-32,-8; -32,11; -16,11; -16,92; 104,92],
                                                     style(rgbcolor={0,0,127}));
    connect(MW2W.y,furnace.Q_flow)       annotation (points=[-84.5,-70; -36,-70;
          -36,-63],          style(rgbcolor={0,0,127}));
    connect(pressure.p, Pa2bar.u) 
      annotation (points=[31,24; 36,24], style(color=74, rgbcolor={0,0,127}));
    connect(Pa2bar.y, p_S) 
      annotation (points=[47.5,24; 104,24],style(color=74, rgbcolor={0,0,127}));
    connect(q_F, MW2W.u) annotation (points=[-104.5,-70; -96,-70],
                                                                 style(color=74,
          rgbcolor={0,0,127}));
    connect(K2degC.Celsius, T_S) annotation (points=[48.5,60; 104,60],style(
          color=74, rgbcolor={0,0,127}));
    connect(controller.y, limiter.u) annotation (points=[-65.7,30; -69.6,30],
        style(color=74, rgbcolor={0,0,127}));
    connect(limiter.y, pump.m_flow_in) annotation (points=[-85.7,30; -90,30;
          -90,-14; -79.3,-14], style(color=74, rgbcolor={0,0,127}));
    connect(temperature.port_b, massFlowRate.port_a) annotation (points=[11,-19;
          19,-19],                      style(color=69, rgbcolor={0,127,255}));
    connect(temperature.T, K2degC.Kelvin) annotation (points=[0,-8; 0,60; 37,60],
               style(color=74, rgbcolor={0,0,127}));
    connect(pressure.port, massFlowRate.port_a) annotation (points=[20,13; 20,
          -20; 19,-19],          style(color=69, rgbcolor={0,127,255}));
    connect(evaporator.steam, temperature.port_a) annotation (points=[-25,-19;
          -11,-19], style(color=69, rgbcolor={0,127,255}));
    connect(pump.port, evaporator.feedwater) annotation (points=[-59,-20; -47,
          -19], style(color=69, rgbcolor={0,127,255}));
    connect(massFlowRate.port_b, SteamValve.port_a) annotation (points=[41,-19;
          47,-19; 47,-18; 52.35,-18], style(color=69, rgbcolor={0,127,255}));
    connect(SteamValve.port_b, sink.port) annotation (points=[66.65,-18; 71,-18;
          71,-19; 76,-19], style(color=69, rgbcolor={0,127,255}));
    connect(SteamValve.opening, Y_Valve) annotation (points=[59.5,-24.4; 59.5,
          -90; -104.5,-90], style(color=74, rgbcolor={0,0,127}));
  end DrumBoiler;
end DrumBoiler;
