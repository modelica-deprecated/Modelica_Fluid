model TestSimpleGenericOrifice 
  extends Modelica.Icons.Example;
  parameter Real D_a = 0.1;
  parameter Real D_b = 0.2;
  parameter Real A_rel = (D_a/D_b)^2;
  parameter Real zeta =  (1 - A_rel)^2;
  
  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
    extends Modelica.Media.Interfaces.PartialMedium "Medium in all components" 
                                                      annotation (
    choicesAllMatching =                                                                            true);
  
  annotation (Diagram,
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Coordsys(extent=[-100,-100; 100,100]),
    Documentation(info="<html>
</html>"));
  Components.Sources.PrescribedAmbient_pTX ambient_a(redeclare package Medium 
      = Medium) 
    annotation (extent=[-40,40; -20,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
    annotation (extent=[-80,40; -60,60]);
  
  Components.Sources.FixedAmbient_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,40; 40,60]);
  Components.Sources.FixedAmbient_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,10; 40,30]);
  Modelica_Fluid.Components.PressureLosses.SimpleGenericOrifice expansion1(
    redeclare package Medium = Medium,
    zeta=zeta,
    diameter=D_a) 
             annotation (extent=[0,40; 20,60]);
  Modelica_Fluid.Components.PressureLosses.SimpleGenericOrifice expansion2(
    redeclare package Medium = Medium,
    from_dp=false,
    zeta=zeta,
    diameter=D_a) annotation (extent=[0,10; 20,30]);
  
  Modelica_Fluid.Components.PressureLosses.SuddenExpansion expansion3(
    redeclare package Medium = Medium,
    use_Re=false,
    D_a=D_a,
    D_b=D_b) annotation (extent=[0,-20; 20,0]);
  Components.Sources.FixedAmbient_pTX ambient_p3(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,-20; 40,0]);
  Modelica_Fluid.Components.PressureLosses.SuddenExpansion expansion4(
    redeclare package Medium = Medium,
    use_Re=false,
    D_a=D_a,
    D_b=D_b,
    from_dp=false) 
             annotation (extent=[0,-60; 20,-40]);
  Components.Sources.FixedAmbient_pTX ambient_p4(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,-60; 40,-40]);
  inner Components.Ambient ambient annotation (extent=[-82,-90; -62,-70]);
equation 
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (points=[-59,50; -52,50; -52,56;
        -42,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient_a.port, expansion1.port_a) 
                                           annotation (points=[-20,50; 0,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion1.port_b, ambient_p1.port) annotation (points=[20,50; 40,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion2.port_b, ambient_p2.port) annotation (points=[20,20; 40,20],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion2.port_a, ambient_a.port) 
                                           annotation (points=[0,20; -10,20;
        -10,50; -20,50], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(ambient_a.port, expansion3.port_a) 
                                           annotation (points=[-20,50; -10,50;
        -10,-10; 0,-10], style(color=69, rgbcolor={0,127,255}));
  connect(expansion3.port_b, ambient_p3.port) annotation (points=[20,-10; 40,
        -10], style(color=69, rgbcolor={0,127,255}));
  connect(ambient_a.port, expansion4.port_a) 
                                           annotation (points=[-20,50; -10,50;
        -10,-50; 0,-50], style(color=69, rgbcolor={0,127,255}));
  connect(expansion4.port_b, ambient_p4.port) annotation (points=[20,-50; 40,
        -50], style(color=69, rgbcolor={0,127,255}));
end TestSimpleGenericOrifice;
